#include "octotigerII/problems/whiteDwarfStructure.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace octotigerII::problems {
namespace {
constexpr Real pi=std::numbers::pi_v<Real>;
Real value(auto q) { return units::value(q); }

void legendre(Real x,int degree,std::vector<Real>& p,std::vector<Real>* derivative=nullptr) {
	p.resize(degree+1);p[0]=1;
	if(derivative){derivative->resize(degree+1);(*derivative)[0]=0;}
	if(degree==0)return;
	p[1]=x;if(derivative)(*derivative)[1]=1;
	for(int l=2;l<=degree;++l){
		p[l]=((2*l-1)*x*p[l-1]-(l-1)*p[l-2])/l;
		if(derivative)(*derivative)[l]=((2*l-1)*(p[l-1]+x*(*derivative)[l-1])-(l-1)*(*derivative)[l-2])/l;
	}
}

void gaussLegendre(int count,std::vector<Real>& nodes,std::vector<Real>& weights) {
	nodes.resize(count);weights.resize(count);
	std::vector<Real> p,dp;
	for(int k=0;k<count;++k){
		Real x=std::cos(pi*(k+Real(.75))/(count+Real(.5)));
		for(int n=0;n<32;++n){
			legendre(x,count,p,&dp);Real const delta=p[count]/dp[count];x-=delta;
			if(std::abs(delta)<4*std::numeric_limits<Real>::epsilon())break;
		}
		legendre(x,count,p,&dp);
		nodes[k]=(1+x)/2;weights[k]=1/((1-x*x)*dp[count]*dp[count]);
	}
}
}

Real WhiteDwarfStructure::densityRatio(Real h) const {
	if(h<=0)return 0;
	auto const coldH=units::VelocitySquared::from_value(value(centralEnthalpy_)*h/(1+parameters_.thermalPressureFraction));
	return value(cold_.densityFromEnthalpy(coldH))/value(parameters_.centralDensity);
}

WhiteDwarfStructure::WhiteDwarfStructure(Parameters p)
	: parameters_(p),cold_(p.meanMassPerElectron) {
	if(!(p.centralDensity>units::Density{}) || !units::finite(p.centralDensity)
		|| !(p.meanMassPerIon>0) || !std::isfinite(p.meanMassPerIon)
		|| !(p.thermalPressureFraction>0 && p.thermalPressureFraction<.01)
		|| !std::isfinite(p.spinFractionOfSphericalBreakup)
		|| p.spinFractionOfSphericalBreakup<0 || p.spinFractionOfSphericalBreakup>.6
		|| p.radialCells<32 || p.angularPoints<4 || p.maxMultipole<0
		|| p.maxMultipole%2 || p.maxMultipole>=2*p.angularPoints || p.maxMultipole>32
		|| p.maxIterations<1 || !(p.tolerance>0) || !(p.relaxation>0 && p.relaxation<=1))
		throw std::invalid_argument("Invalid rotating white-dwarf reference parameters");
	centralEnthalpy_=(1+p.thermalPressureFraction)*cold_.enthalpy(p.centralDensity);
	scaleLength_=units::sqrt(centralEnthalpy_/(4*pi*constants::G*p.centralDensity));
	// Integrate the spherical barotrope first, fixing the breakup normalization.
	struct SphericalNode { Real r,h,m; };
	std::vector<SphericalNode> spherical;
	Real const dr=.002;
	spherical.push_back({0,1,0});
	SphericalNode y{dr,1-dr*dr/6,dr*dr*dr/3};
	spherical.push_back(y);
	auto rhs=[&](Real r,Real h,Real m){return std::array<Real,2>{-m/(r*r),r*r*densityRatio(h)};};
	for(int i=0;i<100000 && y.h>0;++i){
		auto const a=rhs(y.r,y.h,y.m);
		auto const b=rhs(y.r+dr/2,y.h+dr*a[0]/2,y.m+dr*a[1]/2);
		auto const c=rhs(y.r+dr/2,y.h+dr*b[0]/2,y.m+dr*b[1]/2);
		auto const d=rhs(y.r+dr,y.h+dr*c[0],y.m+dr*c[1]);
		y={y.r+dr,y.h+dr*(a[0]+2*b[0]+2*c[0]+d[0])/6,
			y.m+dr*(a[1]+2*b[1]+2*c[1]+d[1])/6};
		spherical.push_back(y);
	}
	if(y.h>0 || spherical.size()<3)throw std::runtime_error("White-dwarf spherical reference has no closed surface");
	auto const& previous=spherical[spherical.size()-2];
	Real const fraction=previous.h/(previous.h-y.h);
	sphericalRadius_=previous.r+fraction*(y.r-previous.r);
	sphericalMass_=previous.m+fraction*(y.m-previous.m);
	omegaSquared_=p.spinFractionOfSphericalBreakup*p.spinFractionOfSphericalBreakup
		*sphericalMass_/std::pow(sphericalRadius_,3);
	outerRadius_=1.6*sphericalRadius_;step_=outerRadius_/p.radialCells;
	gaussLegendre(p.angularPoints,mu_,angularWeights_);
	int const angles=p.angularPoints,count=p.radialCells+1,modes=p.maxMultipole/2+1;
	legendre_.resize(angles*modes);
	std::vector<Real> polynomials;
	for(int k=0;k<angles;++k){
		legendre(mu_[k],p.maxMultipole,polynomials);
		for(int l=0;l<modes;++l)legendre_[k*modes+l]=polynomials[2*l];
	}
	std::vector<Real> density(count*angles),next(density.size());
	for(int j=0;j<count;++j){
		Real const r=j*step_;
		Real ratio=0;
		if(r<sphericalRadius_){
			auto const it=std::lower_bound(spherical.begin(),spherical.end(),r,
				[](SphericalNode const& node,Real radius){return node.r<radius;});
			if(it==spherical.begin())ratio=1;
			else {
				auto const& a=*(it-1);auto const& b=*it;
				ratio=densityRatio(a.h+(b.h-a.h)*(r-a.r)/(b.r-a.r));
			}
		}
		for(int k=0;k<angles;++k)density[j*angles+k]=ratio;
	}
	bool converged=false;
	for(int iteration=0;iteration<p.maxIterations;++iteration){
		poisson(density);centralPotential_=coefficients_[0];
		Real error=0;
		for(int k=0;k<angles;++k){
			bool outside=false;
			for(int j=0;j<count;++j){
				Real phi=0;
				for(int l=0;l<modes;++l)phi+=coefficients_[j*modes+l]*legendre_[k*modes+l];
				Real const r=j*step_,h=1+centralPotential_-phi+
					Real(.5)*omegaSquared_*r*r*(1-mu_[k]*mu_[k]);
				if(h<=0)outside=true;
				Real const predicted=outside?0:densityRatio(h);
				error=std::max(error,std::abs(predicted-density[j*angles+k]));
				next[j*angles+k]=density[j*angles+k]+p.relaxation*(predicted-density[j*angles+k]);
				if(j==count-1 && !outside)
					throw std::runtime_error("Rotating white dwarf exceeds the SCF domain");
			}
		}
		density.swap(next);
		diagnostics_.iterations=iteration+1;diagnostics_.densityResidual=error;
		if(error<p.tolerance){converged=true;break;}
	}
	if(!converged)throw std::runtime_error("Rotating white-dwarf SCF did not converge");
	poisson(density);centralPotential_=coefficients_[0];
	equatorialRadius_=surface(0);polarRadius_=surface(1);
	Real potentialEnergy=0,rotationalEnergy=0,pressureIntegral=0;
	for(int j=0;j<count;++j){
		Real const r=j*step_,radialWeight=(j==0||j==count-1?.5:1)*step_*r*r;
		for(int k=0;k<angles;++k){
			Real phi=0;
			for(int l=0;l<modes;++l)phi+=coefficients_[j*modes+l]*legendre_[k*modes+l];
			Real const ratio=density[j*angles+k],rotation=omegaSquared_*r*r*(1-mu_[k]*mu_[k]);
			if(ratio>1e-12){
				auto const h=cold_.enthalpy(p.centralDensity*ratio)*(1+p.thermalPressureFraction);
				diagnostics_.bernoulliResidual=std::max(diagnostics_.bernoulliResidual,
					std::abs(value(h/centralEnthalpy_)+phi-rotation/2-1-centralPotential_));
			}
			Real const weight=radialWeight*angularWeights_[k];
			potentialEnergy+=weight*ratio*phi/2;
			rotationalEnergy+=weight*ratio*rotation/2;
			if(ratio>0)pressureIntegral+=weight*value((1+p.thermalPressureFraction)
				*cold_.pressure(p.centralDensity*ratio)/(p.centralDensity*centralEnthalpy_));
		}
	}
	diagnostics_.virialResidual=(2*rotationalEnergy+potentialEnergy+3*pressureIntegral)/std::abs(potentialEnergy);
}

void WhiteDwarfStructure::poisson(std::vector<Real> const& density) {
	int const count=parameters_.radialCells+1,angles=parameters_.angularPoints,modes=parameters_.maxMultipole/2+1;
	coefficients_.assign(count*modes,0);coefficientDerivatives_.assign(count*modes,0);exteriorMoments_.resize(modes);
	std::vector<Real> q(count),interior(count),exterior(count);
	for(int mode=0;mode<modes;++mode){
		int const l=2*mode;
		for(int j=0;j<count;++j){
			q[j]=0;
			for(int k=0;k<angles;++k)q[j]+=angularWeights_[k]*density[j*angles+k]*legendre_[k*modes+mode];
		}
		if(l>0)q[0]=0;
		interior[0]=0;
		for(int j=1;j<count;++j){
			Real const a=(j-1)*step_,b=j*step_;
			interior[j]=interior[j-1]+step_*(q[j-1]*std::pow(a,l+2)+q[j]*std::pow(b,l+2))/2;
		}
		exterior[count-1]=0;
		for(int j=count-2;j>=0;--j){
			Real const a=j*step_,b=(j+1)*step_;
			exterior[j]=exterior[j+1]+step_*((j==0?0:q[j]*std::pow(a,1-l))+
				q[j+1]*std::pow(b,1-l))/2;
		}
		exteriorMoments_[mode]=interior.back();
		if(l==0){coefficients_[0]=-exterior[0];massIntegral_=interior.back();}
		for(int j=1;j<count;++j){
			Real const r=j*step_,inner=interior[j]/std::pow(r,l+1),outer=std::pow(r,l)*exterior[j];
			coefficients_[j*modes+mode]=-inner-outer;
			coefficientDerivatives_[j*modes+mode]=((l+1)*inner-l*outer)/r;
		}
	}
}

WhiteDwarfStructure::Potential WhiteDwarfStructure::potential(Real r,Real mu) const {
	int const modes=parameters_.maxMultipole/2+1;
	std::array<Real,33> p{},dp{};
	p[0]=1;p[1]=mu;dp[1]=1;
	for(int l=2;l<=parameters_.maxMultipole;++l){
		p[l]=((2*l-1)*mu*p[l-1]-(l-1)*p[l-2])/l;
		dp[l]=((2*l-1)*(p[l-1]+mu*dp[l-1])-(l-1)*dp[l-2])/l;
	}
	Potential result;
	int const cell=r>=outerRadius_?0:std::clamp(int(r/step_),0,parameters_.radialCells-1);
	Real const t=r>=outerRadius_?0:(r-cell*step_)/step_,t2=t*t,t3=t2*t;
	for(int mode=0;mode<modes;++mode){
		int const l=2*mode;
		Real y,dy;
		if(r>=outerRadius_){
			y=-exteriorMoments_[mode]/std::pow(r,l+1);dy=-(l+1)*y/r;
		}else{
			Real const a=coefficients_[cell*modes+mode],b=coefficients_[(cell+1)*modes+mode];
			Real const da=coefficientDerivatives_[cell*modes+mode],db=coefficientDerivatives_[(cell+1)*modes+mode];
			y=(2*t3-3*t2+1)*a+(t3-2*t2+t)*step_*da+(-2*t3+3*t2)*b+(t3-t2)*step_*db;
			dy=(6*t2-6*t)*(a-b)/step_+(3*t2-4*t+1)*da+(3*t2-2*t)*db;
		}
		result.value+=y*p[l];result.radialDerivative+=dy*p[l];result.muDerivative+=y*dp[l];
	}
	return result;
}

Real WhiteDwarfStructure::surface(Real mu) const {
	auto h=[&](Real r){return 1+centralPotential_-potential(r,mu).value+
		Real(.5)*omegaSquared_*r*r*(1-mu*mu);};
	for(int j=1;j<=parameters_.radialCells;++j)if(h(j*step_)<=0){
		Real lower=(j-1)*step_,upper=j*step_;
		for(int iteration=0;iteration<64;++iteration){
			Real const mid=(lower+upper)/2;
			if(h(mid)>0)lower=mid;else upper=mid;
		}
		return (lower+upper)/2;
	}
	throw std::runtime_error("Rotating white dwarf has no closed surface in the SCF domain");
}

WhiteDwarfStructure::State WhiteDwarfStructure::sample(units::Length cylindricalRadius,units::Length z) const {
	if(!units::finite(cylindricalRadius)||!units::finite(z)||cylindricalRadius<units::Length{})
		throw std::invalid_argument("Invalid white-dwarf sample position");
	Real const scale=value(scaleLength_),R=value(cylindricalRadius)/scale,Z=value(z)/scale;
	Real const r=std::hypot(R,Z),mu=r==0?0:Z/r,sintheta=r==0?0:R/r;
	if(!std::isfinite(r))throw std::invalid_argument("White-dwarf sample radius is not finite");
	auto const phi=potential(r,mu);
	Real const normalizedH=r>equatorialRadius_?0:
		std::max(Real(0),1+centralPotential_-phi.value+Real(.5)*omegaSquared_*R*R);
	State result;
	result.density=parameters_.centralDensity*densityRatio(normalizedH);
	if(result.density>units::Density{}){
		auto const coldPressure=cold_.pressure(result.density);
		result.ionPressure=parameters_.thermalPressureFraction*coldPressure;
		result.pressure=coldPressure+result.ionPressure;
		result.temperature=result.ionPressure*(parameters_.meanMassPerIon*constants::atomicMassUnit)/
			(result.density*constants::boltzmann);
		result.radiationEnergy=units::EnergyDensity::from_value(value(constants::radiation)*
			std::pow(value(result.temperature),4));
	}
	result.potential=phi.value*centralEnthalpy_;
	Real const radial=r==0?0:phi.radialDerivative*sintheta-phi.muDerivative*mu*sintheta/r;
	Real const vertical=r==0?0:phi.radialDerivative*mu+phi.muDerivative*(1-mu*mu)/r;
	result.potentialGradient={units::Acceleration::from_value(value(centralEnthalpy_)/scale*radial),
		units::Acceleration::from_value(value(centralEnthalpy_)/scale*vertical)};
	result.effectivePotentialGradient={units::Acceleration::from_value(value(centralEnthalpy_)/scale*
		(radial-omegaSquared_*R)),result.potentialGradient[1]};
	return result;
}

std::array<units::EnergyFlux,2> WhiteDwarfStructure::diffusionFlux(State const& state,Real opacityCgs) const {
	if(!(opacityCgs>0)||!std::isfinite(opacityCgs))
		throw std::invalid_argument("White-dwarf diffusion opacity must be positive and finite");
	std::array<units::EnergyFlux,2> result{};
	if(state.density==units::Density{})return result;
	auto const coldPressure=cold_.pressure(state.density);
	auto const derivative=cold_.pressureDerivative(state.density);
	Real const logarithmicSlope=value(state.density*derivative/coldPressure)-1;
	Real const dEdH=4*value(state.radiationEnergy)/
		((1+parameters_.thermalPressureFraction)*value(derivative))*logarithmicSlope;
	Real const factor=value(constants::c)/(3*opacityCgs*value(state.density));
	for(int axis=0;axis<2;++axis)
		result[axis]=units::EnergyFlux::from_value(factor*dEdH*value(state.effectivePotentialGradient[axis]));
	return result;
}

units::Mass WhiteDwarfStructure::mass() const {
	return units::Mass::from_value(4*pi*value(parameters_.centralDensity)*std::pow(value(scaleLength_),3)*massIntegral_);
}
units::InverseTime WhiteDwarfStructure::angularVelocity() const {
	return units::InverseTime::from_value(std::sqrt(omegaSquared_*4*pi*value(constants::G)*value(parameters_.centralDensity)));
}
units::Velocity WhiteDwarfStructure::centralSoundSpeed() const {
	auto const rho=parameters_.centralDensity;
	auto const coldSoundSquared=cold_.pressureDerivative(rho);
	auto const ionPressure=parameters_.thermalPressureFraction*cold_.pressure(rho);
	return units::sqrt(coldSoundSquared+Real(5)/3*ionPressure/rho);
}
} // namespace octotigerII::problems
