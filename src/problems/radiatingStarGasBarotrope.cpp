#include "octotigerII/problems/radiatingStarGasBarotrope.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace octotigerII::problems {
namespace {
Real value(auto q) { return units::value(q); }
constexpr Real maximumDensityRatio = 4;
// Eight-point Gauss-Legendre quadrature, paired about each cell midpoint.
constexpr std::array<Real,4> nodes{0.18343464249564980494,0.52553240991632898582,
	0.79666647741362673959,0.96028985649753623168};
constexpr std::array<Real,4> weights{0.36268378337836198297,0.31370664587788728734,
	0.22238103445337447054,0.10122853629037625915};
}

RadiatingStarGasBarotrope::RadiatingStarGasBarotrope(Parameters parameters)
	: parameters_(parameters), eos_(parameters.eos) {
	if (!std::isfinite(parameters.cutoffDensityFraction)
		|| !(parameters.cutoffDensityFraction>0 && parameters.cutoffDensityFraction<1))
		throw std::invalid_argument("Gas barotrope requires 0 < cutoff density fraction < 1");
	gasConstant_=value(constants::boltzmann/constants::atomicMassUnit)/parameters.eos.meanMolecularWeight;
	hScale_=gasConstant_*value(eos_.centralTemperature());
	auto const cut=eos_.atDensity(cutoffDensity());
	Real const gamma=1+1/parameters.eos.index;
	envelopeGamma_=1+(gamma-cut.beta)/(4-3*cut.beta);
	cutoffPressure_=value(cut.gasPressure);
	cutoffH_=envelopeGamma_/(envelopeGamma_-1)*cutoffPressure_/value(cut.density)/hScale_;
	if (!(cutoffH_>0) || !std::isfinite(cutoffH_))
		throw std::invalid_argument("Gas barotrope cutoff is outside the representable range");

	Real const lower=std::log(parameters.cutoffDensityFraction), upper=std::log(maximumDensityRatio);
	// The small log-density spacing makes interpolation (not quadrature) the
	// limiting error. Exact endpoint slopes retain the matched envelope derivative.
	int const cells=std::max(1,static_cast<int>(std::ceil((upper-lower)/0.002)));
	table_.resize(cells+1);
	long double accumulated=cutoffH_;
	for (int i=0;i<=cells;++i) {
		Real const x=lower+(upper-lower)*Real(i)/cells;
		if (i) {
			Real const midpoint=(table_[i-1].x+x)/2, half=(x-table_[i-1].x)/2;
			long double integral=0;
			for (std::size_t k=0;k<nodes.size();++k)
				integral+=weights[k]*(coreDerivative(midpoint-half*nodes[k])+coreDerivative(midpoint+half*nodes[k]));
			accumulated+=half*integral;
		}
		table_[i]={x,static_cast<Real>(accumulated),coreDerivative(x)};
		if (i) {
			// This sufficient monotonicity condition also rejects nonrepresentable
			// increments. No extrapolation or nonmonotone inverse is accepted.
			Real const secant=(table_[i].h-table_[i-1].h)/(x-table_[i-1].x);
			Real const a=table_[i-1].derivative/secant,b=table_[i].derivative/secant;
			if (!(secant>0 && a>0 && b>0 && a*a+b*b<=9))
				throw std::runtime_error("Gas barotrope interpolation lost monotonicity");
		}
	}
}

Real RadiatingStarGasBarotrope::coreDerivative(Real x) const {
	auto const q=eos_.atDensity(parameters_.eos.centralDensity*std::exp(x));
	Real const slope=(1+1/parameters_.eos.index-q.beta)/(4-3*q.beta);
	// d h_g / d log rho = d p_g / d rho = R T (1+d log T/d log rho).
	return gasConstant_*value(q.temperature)*(1+slope)/hScale_;
}

RadiatingStarGasBarotrope::Interpolation RadiatingStarGasBarotrope::interpolate(std::size_t i, Real x) const {
	auto const& a=table_[i];auto const& b=table_[i+1];
	Real const width=b.x-a.x,t=(x-a.x)/width,change=b.h-a.h;
	Real const start=width*a.derivative,end=width*b.derivative;
	Real const quadratic=3*change-2*start-end,cubic=start+end-2*change;
	return {a.h+t*(start+t*(quadratic+t*cubic)),(start+t*(2*quadratic+3*t*cubic))/width};
}

RadiatingStarGasBarotrope::State RadiatingStarGasBarotrope::coreState(units::Density density, Real h, Real derivative) const {
	auto const q=eos_.atDensity(density);
	return {density,q.temperature,q.gasPressure,units::VelocitySquared::from_value(h*hScale_),
		units::Quantity<-5,1,2>::from_value(value(density)/(hScale_*derivative))};
}

RadiatingStarGasBarotrope::State RadiatingStarGasBarotrope::atDensity(units::Density density) const {
	if (!units::finite(density) || density<units::Density{})
		throw std::invalid_argument("Gas barotrope density must be finite and nonnegative");
	if (density>maximumDensity()) throw std::out_of_range("Gas barotrope density exceeds 4 rho_c");
	if (density==units::Density{}) return {};
	if (density==maximumDensity()) return coreState(density,table_.back().h,table_.back().derivative);
	if (density<=cutoffDensity()) {
		Real const ratio=Real(density/cutoffDensity()),power=std::pow(ratio,envelopeGamma_-1);
		Real const h=cutoffH_*hScale_*power,pressure=cutoffPressure_*ratio*power;
		Real const temperature=(envelopeGamma_-1)*h/(envelopeGamma_*gasConstant_);
		return {density,units::Temperature::from_value(temperature),
			units::Pressure::from_value(pressure),units::VelocitySquared::from_value(h),
			units::Quantity<-5,1,2>::from_value(value(density)/((envelopeGamma_-1)*h))};
	}
	Real const x=std::log(Real(density/parameters_.eos.centralDensity));
	auto const next=std::upper_bound(table_.begin(),table_.end(),x,[](Real key,Node const& node){return key<node.x;});
	std::size_t const i=std::min(std::size_t(next-table_.begin()-1),table_.size()-2);
	auto const q=interpolate(i,x);
	return coreState(density,q.h,q.derivative);
}

RadiatingStarGasBarotrope::State RadiatingStarGasBarotrope::atIntegralH(units::VelocitySquared integralH) const {
	if (!units::finite(integralH)) throw std::invalid_argument("Gas barotrope integral must be finite");
	if (integralH<=units::VelocitySquared{}) return {};
	if (integralH>maximumIntegralH()) throw std::out_of_range("Gas barotrope integral exceeds its 4 rho_c table");
	if (integralH==maximumIntegralH()) return atDensity(maximumDensity());
	Real const h=value(integralH)/hScale_;
	if (h<=cutoffH_) {
		Real const density=value(cutoffDensity())*std::pow(h/cutoffH_,1/(envelopeGamma_-1));
		Real const temperature=(envelopeGamma_-1)*value(integralH)/(envelopeGamma_*gasConstant_);
		return {units::Density::from_value(density),units::Temperature::from_value(temperature),
			units::Pressure::from_value(density*gasConstant_*temperature),integralH,
			units::Quantity<-5,1,2>::from_value(density/((envelopeGamma_-1)*value(integralH)))};
	}
	auto const next=std::upper_bound(table_.begin(),table_.end(),h,[](Real key,Node const& node){return key<node.h;});
	std::size_t const i=std::min(std::size_t(next-table_.begin()-1),table_.size()-2);
	Real lower=table_[i].x,upper=table_[i+1].x;
	Real x=lower+(upper-lower)*(h-table_[i].h)/(table_[i+1].h-table_[i].h);
	for (int iteration=0;iteration<16;++iteration) {
		auto const q=interpolate(i,x);
		Real const residual=q.h-h;
		if (std::abs(residual)<=2*std::numeric_limits<Real>::epsilon()*h) break;
		if (residual>0) upper=x; else lower=x;
		Real const candidate=x-residual/q.derivative;
		x=candidate>lower && candidate<upper ? candidate : (lower+upper)/2;
	}
	auto const q=interpolate(i,x);
	auto result=coreState(parameters_.eos.centralDensity*std::exp(x),h,q.derivative);
	result.integralH=integralH;
	return result;
}

units::Density RadiatingStarGasBarotrope::cutoffDensity() const {return parameters_.cutoffDensityFraction*parameters_.eos.centralDensity;}
units::Density RadiatingStarGasBarotrope::maximumDensity() const {return maximumDensityRatio*parameters_.eos.centralDensity;}
units::VelocitySquared RadiatingStarGasBarotrope::cutoffIntegralH() const {return units::VelocitySquared::from_value(cutoffH_*hScale_);}
units::VelocitySquared RadiatingStarGasBarotrope::centralIntegralH() const {return atDensity(parameters_.eos.centralDensity).integralH;}
units::VelocitySquared RadiatingStarGasBarotrope::maximumIntegralH() const {return units::VelocitySquared::from_value(table_.back().h*hScale_);}

} // namespace octotigerII::problems
