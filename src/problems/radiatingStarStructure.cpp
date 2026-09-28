#include "octotigerII/problems/radiatingStarStructure.hpp"
#include "octotigerII/problems/laneEmden.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace octotigerII::problems {
namespace {
constexpr Real pi = std::numbers::pi_v<Real>;
Real value(auto quantity) { return units::value(quantity); }
bool positive(Real x) { return std::isfinite(x) && x > 0; }

void legendre(Real x, int degree, std::vector<Real>& p, std::vector<Real>* derivative = nullptr) {
	p.resize(degree + 1); p[0] = 1;
	if (derivative) { derivative->resize(degree + 1); (*derivative)[0] = 0; }
	if (degree == 0) return;
	p[1] = x;
	if (derivative) (*derivative)[1] = 1;
	for (int l = 2; l <= degree; ++l) {
		p[l] = ((2*l-1)*x*p[l-1] - (l-1)*p[l-2])/l;
		if (derivative) (*derivative)[l] = ((2*l-1)*(p[l-1]+x*(*derivative)[l-1]) - (l-1)*(*derivative)[l-2])/l;
	}
}

void gaussLegendre(int count, std::vector<Real>& nodes, std::vector<Real>& weights) {
	nodes.resize(count); weights.resize(count);
	std::vector<Real> p, dp;
	for (int k = 0; k < count; ++k) {
		Real x = std::cos(pi*(k+0.75)/(count+0.5));
		for (int iteration = 0; iteration < 32; ++iteration) {
			legendre(x, count, p, &dp);
			Real const change = p[count]/dp[count];
			x -= change;
			if (std::abs(change) < 4*std::numeric_limits<Real>::epsilon()) break;
		}
		legendre(x, count, p, &dp);
		nodes[k] = (1+x)/2;
		weights[k] = 1/((1-x*x)*dp[count]*dp[count]);
	}
}
}

RadiatingStarEos::RadiatingStarEos(Parameters parameters) : parameters_(parameters) {
	if (!(parameters.index > 3 && parameters.index < 5) || !std::isfinite(parameters.index)
		|| !positive(value(parameters.centralDensity)) || !(parameters.centralBeta > 0 && parameters.centralBeta < 1)
		|| !positive(parameters.meanMolecularWeight)) throw std::invalid_argument("Radiating-star EOS requires 3<n<5, positive density/mu, and 0<beta<1");
	gasConstant_ = value(constants::boltzmann/constants::atomicMassUnit)/parameters.meanMolecularWeight;
	centralTemperature_ = std::cbrt(3*gasConstant_*value(parameters.centralDensity)*(1-parameters.centralBeta)
		/(value(constants::radiation)*parameters.centralBeta));
	centralPressure_ = value(parameters.centralDensity)*gasConstant_*centralTemperature_/parameters.centralBeta;
	centralH_ = (parameters.index+1)*centralPressure_/value(parameters.centralDensity);
	scaleLength_ = std::sqrt(centralH_/(4*pi*value(constants::G)*value(parameters.centralDensity)));
	if (!positive(centralTemperature_) || !positive(centralPressure_) || !positive(centralH_) || !positive(scaleLength_))
		throw std::invalid_argument("Radiating-star EOS scaling is outside the representable range");
}

RadiatingStarEos::State RadiatingStarEos::atDensity(units::Density density) const {
	Real const rho = value(density);
	if (!std::isfinite(rho) || rho < 0) throw std::invalid_argument("Radiating-star density must be finite and nonnegative");
	if (rho == 0) return {};
	Real const ratio = rho/value(parameters_.centralDensity), n = parameters_.index, gamma = 1+1/n;
	Real const pressure = centralPressure_*std::pow(ratio, gamma);
	if (!positive(pressure)) throw std::overflow_error("Radiating-star pressure is outside the representable range");
	// Bracket the unique temperature root by its pure gas and pure radiation
	// limits, then solve a dimensionless polynomial with coefficients <=1.
	Real const upper = std::min(pressure/(rho*gasConstant_), std::pow(3*pressure/value(constants::radiation), Real(0.25)));
	Real const gas = rho*gasConstant_*upper/pressure;
	Real const radiation = value(constants::radiation)*std::pow(upper, 4)/(3*pressure);
	Real lower = 0, higher = 1, t = 0.8;
	for (int iteration = 0; iteration < 64; ++iteration) {
		Real const f = gas*t+radiation*std::pow(t,4)-1;
		if (std::abs(f) < 4*std::numeric_limits<Real>::epsilon()) break;
		if (f > 0) higher = t; else lower = t;
		Real next = t-f/(gas+4*radiation*t*t*t);
		if (!(next > lower && next < higher)) next = (lower+higher)/2;
		t = next;
	}
	Real const temperature = t*upper, beta = rho*gasConstant_*temperature/pressure;
	Real const h = centralH_*std::pow(ratio,1/n), slope = (gamma-beta)/(4-3*beta);
	Real const b = 4*(1-beta)*(gamma-beta)/(gamma*(4-3*beta));
	Real const db = -4*(3*beta*beta-8*beta+gamma+4)/(gamma*std::pow(4-3*beta,2));
	Real const betaH = (n-3)*beta*(1-beta)/(h*(4-3*beta));
	Real const radiationEntropy = 4*value(constants::radiation)*std::pow(temperature,3)/(3*rho);
	Real const centralRadiationEntropy = 4*gasConstant_*(1-parameters_.centralBeta)/parameters_.centralBeta;
	State state;
	state.density = density;
	state.pressure = units::Pressure::from_value(pressure);
	state.gasPressure = units::Pressure::from_value(rho*gasConstant_*temperature);
	state.radiationEnergy = units::EnergyDensity::from_value(value(constants::radiation)*std::pow(temperature,4));
	state.temperature = units::Temperature::from_value(temperature);
	state.integralH = units::VelocitySquared::from_value(h);
	state.beta = beta; state.radiationPressureDerivative = b; state.radiationPressureSecondDerivative = db*betaH;
	state.gamma1 = beta+std::pow(4-3*beta,2)/(12-10.5*beta);
	state.nabla = slope/gamma; state.nablaAd = (8-6*beta)/(32-24*beta-3*beta*beta);
	state.entropyDifference = gasConstant_*(1.5*std::log(temperature/centralTemperature_)-std::log(ratio))
		+ radiationEntropy-centralRadiationEntropy;
	state.entropyDerivativeLogDensity = gasConstant_*(1.5*slope-1)+radiationEntropy*(3*slope-1);
	return state;
}

RadiatingStarEos::State RadiatingStarEos::atIntegralH(units::VelocitySquared h) const {
	if (!units::finite(h) || h < units::VelocitySquared{}) throw std::invalid_argument("Radiating-star H must be finite and nonnegative");
	return atDensity(parameters_.centralDensity*std::pow(value(h)/centralH_, parameters_.index));
}
units::Pressure RadiatingStarEos::centralPressure() const { return units::Pressure::from_value(centralPressure_); }
units::Temperature RadiatingStarEos::centralTemperature() const { return units::Temperature::from_value(centralTemperature_); }
units::VelocitySquared RadiatingStarEos::centralIntegralH() const { return units::VelocitySquared::from_value(centralH_); }
units::Length RadiatingStarEos::scaleLength() const { return units::Length::from_value(scaleLength_); }

RadiatingStarStructure::RadiatingStarStructure(Parameters parameters) : parameters_(parameters), eos_(parameters.eos) {
	if (!std::isfinite(parameters.spinFractionOfSphericalBreakup) || parameters.spinFractionOfSphericalBreakup < 0
		|| parameters.spinFractionOfSphericalBreakup > 0.6 || parameters.radialCells < 32 || parameters.angularPoints < 4
		|| parameters.maxMultipole < 0 || parameters.maxMultipole%2 || parameters.maxMultipole >= 2*parameters.angularPoints
		|| parameters.maxMultipole > 32 || parameters.maxIterations < 1 || !positive(parameters.tolerance)
		|| !(parameters.relaxation > 0 && parameters.relaxation <= 1))
		throw std::invalid_argument("Invalid radiating-star SCF resolution, spin, or convergence controls");
	LaneEmden spherical(parameters.eos.index);
	sphericalRadius_ = spherical.surface();
	omegaSquared_ = std::pow(parameters.spinFractionOfSphericalBreakup,2)*spherical.surfaceMass()/std::pow(sphericalRadius_,3);
	outerRadius_ = 1.6*sphericalRadius_; step_ = outerRadius_/parameters.radialCells;
	gaussLegendre(parameters.angularPoints, mu_, angularWeights_);
	int const angles = parameters.angularPoints, count = parameters.radialCells+1, modes = parameters.maxMultipole/2+1;
	legendre_.resize(angles*modes);
	std::vector<Real> p;
	for (int k = 0; k < angles; ++k) {
		legendre(mu_[k], parameters.maxMultipole, p);
		for (int l = 0; l < modes; ++l) legendre_[k*modes+l] = p[2*l];
	}
	std::vector<Real> density(count*angles), next(density.size());
	for (int j = 0; j < count; ++j) {
		Real const initial = std::pow(spherical(j*step_).theta, parameters.eos.index);
		for (int k = 0; k < angles; ++k) density[j*angles+k] = initial;
	}
	bool converged = false;
	for (int iteration = 0; iteration < parameters.maxIterations; ++iteration) {
		poisson(density);
		centralPotential_ = coefficients_[0];
		Real error = 0;
		for (int k = 0; k < angles; ++k) {
			bool outside = false;
			for (int j = 0; j < count; ++j) {
				Real phi = 0;
				for (int l = 0; l < modes; ++l) phi += coefficients_[j*modes+l]*legendre_[k*modes+l];
				Real const r = j*step_, h = 1+centralPotential_-phi+0.5*omegaSquared_*r*r*(1-mu_[k]*mu_[k]);
				if (h <= 0) outside = true;
				Real const predicted = outside ? 0 : std::pow(h, parameters.eos.index);
				error = std::max(error,std::abs(predicted-density[j*angles+k]));
				next[j*angles+k] = density[j*angles+k]+parameters.relaxation*(predicted-density[j*angles+k]);
				if (j == count-1 && !outside) throw std::runtime_error("Rotating SCF star exceeds the reference domain; reduce spin");
			}
		}
		density.swap(next);
		diagnostics_.iterations = iteration+1; diagnostics_.densityResidual = error;
		if (error < parameters.tolerance) { converged = true; break; }
	}
	if (!converged) throw std::runtime_error("Radiating-star SCF did not converge");
	poisson(density); centralPotential_ = coefficients_[0];
	equatorialRadius_ = surface(0); polarRadius_ = surface(1);
	diagnostics_.bernoulliResidual = 0;
	Real potentialEnergy = 0, rotationalEnergy = 0, pressureIntegral = 0;
	for (int j = 0; j < count; ++j) {
		Real const r = j*step_, radialWeight = (j == 0 || j == count-1 ? 0.5 : 1)*step_*r*r;
		for (int k = 0; k < angles; ++k) {
			Real phi = 0;
			for (int l = 0; l < modes; ++l) phi += coefficients_[j*modes+l]*legendre_[k*modes+l];
			Real const rho = density[j*angles+k], rotation = omegaSquared_*r*r*(1-mu_[k]*mu_[k]);
			if (rho > 1e-12) diagnostics_.bernoulliResidual = std::max(diagnostics_.bernoulliResidual,
				std::abs(std::pow(rho,1/parameters.eos.index)+phi-0.5*rotation-1-centralPotential_));
			Real const weight = radialWeight*angularWeights_[k];
			potentialEnergy += weight*rho*phi/2;
			rotationalEnergy += weight*rho*rotation/2;
			pressureIntegral += weight*std::pow(rho,1+1/parameters.eos.index)/(parameters.eos.index+1);
		}
	}
	diagnostics_.virialResidual = (2*rotationalEnergy+potentialEnergy+3*pressureIntegral)/std::abs(potentialEnergy);
}

void RadiatingStarStructure::poisson(std::vector<Real> const& density) {
	int const count = parameters_.radialCells+1, angles = parameters_.angularPoints, modes = parameters_.maxMultipole/2+1;
	coefficients_.assign(count*modes,0); coefficientDerivatives_.assign(count*modes,0); exteriorMoments_.resize(modes);
	std::vector<Real> q(count), interior(count), exterior(count);
	for (int mode = 0; mode < modes; ++mode) {
		int const l = 2*mode;
		for (int j = 0; j < count; ++j) {
			q[j] = 0;
			for (int k = 0; k < angles; ++k) q[j] += angularWeights_[k]*density[j*angles+k]*legendre_[k*modes+mode];
		}
		if (l > 0) q[0] = 0;
		interior[0] = 0;
		for (int j = 1; j < count; ++j) {
			Real const a = (j-1)*step_, b = j*step_;
			interior[j] = interior[j-1]+0.5*step_*(q[j-1]*std::pow(a,l+2)+q[j]*std::pow(b,l+2));
		}
		exterior[count-1] = 0;
		for (int j = count-2; j >= 0; --j) {
			Real const a = j*step_, b = (j+1)*step_;
			exterior[j] = exterior[j+1]+0.5*step_*((j == 0 ? 0 : q[j]*std::pow(a,1-l))+q[j+1]*std::pow(b,1-l));
		}
		exteriorMoments_[mode] = interior.back();
		if (l == 0) { coefficients_[0] = -exterior[0]; massIntegral_ = interior.back(); }
		for (int j = 1; j < count; ++j) {
			Real const r = j*step_, inner = interior[j]/std::pow(r,l+1), outer = std::pow(r,l)*exterior[j];
			coefficients_[j*modes+mode] = -inner-outer;
			coefficientDerivatives_[j*modes+mode] = ((l+1)*inner-l*outer)/r;
		}
	}
}

RadiatingStarStructure::Potential RadiatingStarStructure::potential(Real r, Real mu) const {
	int const modes = parameters_.maxMultipole/2+1;
	// Small fixed arrays keep the repeated cell-query path allocation free.
	std::array<Real,33> p{}, dp{};
	p[0] = 1; p[1] = mu; dp[1] = 1;
	for (int l = 2; l <= parameters_.maxMultipole; ++l) {
		p[l] = ((2*l-1)*mu*p[l-1]-(l-1)*p[l-2])/l;
		dp[l] = ((2*l-1)*(p[l-1]+mu*dp[l-1])-(l-1)*dp[l-2])/l;
	}
	Potential result;
	int const cell = r >= outerRadius_ ? 0 : std::clamp(int(r/step_),0,parameters_.radialCells-1);
	Real const t = r >= outerRadius_ ? 0 : (r-cell*step_)/step_, t2 = t*t, t3 = t2*t;
	for (int mode = 0; mode < modes; ++mode) {
		int const l = 2*mode;
		Real y, dy;
		if (r >= outerRadius_) {
			y = -exteriorMoments_[mode]/std::pow(r,l+1); dy = -(l+1)*y/r;
		} else {
			Real const a = coefficients_[cell*modes+mode], b = coefficients_[(cell+1)*modes+mode];
			Real const da = coefficientDerivatives_[cell*modes+mode], db = coefficientDerivatives_[(cell+1)*modes+mode];
			y = (2*t3-3*t2+1)*a+(t3-2*t2+t)*step_*da+(-2*t3+3*t2)*b+(t3-t2)*step_*db;
			dy = (6*t2-6*t)*(a-b)/step_+(3*t2-4*t+1)*da+(3*t2-2*t)*db;
		}
		result.value += y*p[l]; result.radialDerivative += dy*p[l]; result.muDerivative += y*dp[l];
	}
	return result;
}

Real RadiatingStarStructure::surface(Real mu) const {
	auto h = [&](Real r) { return 1+centralPotential_-potential(r,mu).value+0.5*omegaSquared_*r*r*(1-mu*mu); };
	for (int j = 1; j <= parameters_.radialCells; ++j) if (h(j*step_) <= 0) {
		Real lower = (j-1)*step_, upper = j*step_;
		for (int iteration = 0; iteration < 64; ++iteration) {
			Real const mid = (lower+upper)/2;
			if (h(mid) > 0) lower = mid; else upper = mid;
		}
		return (lower+upper)/2;
	}
	throw std::runtime_error("Rotating SCF has no closed surface in the reference domain");
}

RadiatingStarStructure::State RadiatingStarStructure::sample(units::Length cylindricalRadius, units::Length z) const {
	if (!units::finite(cylindricalRadius) || !units::finite(z) || cylindricalRadius < units::Length{})
		throw std::invalid_argument("Invalid radiating-star sample position");
	Real const scale = value(eos_.scaleLength()), R = value(cylindricalRadius)/scale, Z = value(z)/scale;
	Real const r = std::hypot(R,Z), mu = r == 0 ? 0 : Z/r, sintheta = r == 0 ? 0 : R/r;
	if (!std::isfinite(r)) throw std::invalid_argument("Radiating-star position exceeds the dimensionless range");
	auto const phi = potential(r,mu);
	Real const Hc = value(eos_.centralIntegralH()), h = r > equatorialRadius_ ? 0 : std::max(Real(0),1+centralPotential_-phi.value+0.5*omegaSquared_*R*R);
	State result;
	result.thermodynamics = eos_.atIntegralH(units::VelocitySquared::from_value(Hc*h));
	result.potential = units::VelocitySquared::from_value(Hc*phi.value);
	Real const radial = r == 0 ? 0 : phi.radialDerivative*sintheta-phi.muDerivative*mu*sintheta/r;
	Real const vertical = r == 0 ? 0 : phi.radialDerivative*mu+phi.muDerivative*(1-mu*mu)/r;
	result.potentialGradient = {units::Acceleration::from_value(Hc/scale*radial), units::Acceleration::from_value(Hc/scale*vertical)};
	result.effectivePotentialGradient = {units::Acceleration::from_value(Hc/scale*(radial-omegaSquared_*R)), result.potentialGradient[1]};
	return result;
}

RadiatingStarStructure::Diffusion RadiatingStarStructure::diffusion(State const& state, Real opacityCgs) const {
	if (!positive(opacityCgs)) throw std::invalid_argument("Diffusion-reference opacity must be positive and finite");
	Diffusion result;
	if (state.thermodynamics.density == units::Density{}) return result;
	auto const& q = state.thermodynamics;
	Real const factor = value(constants::c)/opacityCgs, b = q.radiationPressureDerivative;
	Real g2 = 0;
	for (int d = 0; d < 2; ++d) {
		Real const g = value(state.effectivePotentialGradient[d]); g2 += g*g;
		result.flux[d] = units::EnergyFlux::from_value(factor*b*g);
	}
	Real const omega = value(angularVelocity());
	result.heating = factor*(b*(4*pi*value(constants::G)*value(q.density)-2*omega*omega)-q.radiationPressureSecondDerivative*g2);
	result.fluxFactor = factor*b*std::sqrt(g2)/(value(constants::c)*value(q.radiationEnergy));
	return result;
}
units::Mass RadiatingStarStructure::mass() const { return units::Mass::from_value(4*pi*value(parameters_.eos.centralDensity)*std::pow(value(eos_.scaleLength()),3)*massIntegral_); }
units::Length RadiatingStarStructure::equatorialRadius() const { return equatorialRadius_*eos_.scaleLength(); }
units::Length RadiatingStarStructure::polarRadius() const { return polarRadius_*eos_.scaleLength(); }
units::Length RadiatingStarStructure::sphericalRadius() const { return sphericalRadius_*eos_.scaleLength(); }
units::InverseTime RadiatingStarStructure::angularVelocity() const { return units::InverseTime::from_value(std::sqrt(omegaSquared_*4*pi*value(constants::G)*value(parameters_.eos.centralDensity))); }
units::VelocitySquared RadiatingStarStructure::bernoulliConstant() const { return (1+centralPotential_)*eos_.centralIntegralH(); }
} // namespace octotigerII::problems
