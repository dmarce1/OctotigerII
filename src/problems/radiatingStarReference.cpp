// Deferred table-reader prototype; deliberately absent from the build.
#include "octotigerII/problems/radiatingStarReference.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace octotigerII::problems {
namespace {
Real value(auto quantity) { return units::value(quantity); }

RadiatingStarEos::Parameters eosParameters(RadiatingStarReference::Parameters const& parameters,
	RadiatingStarReference::Metadata const& metadata) {
	if (parameters.centralDensity != units::Density::from_value(1)
		|| parameters.meanMolecularWeight != Real(0.6)) {
		throw std::invalid_argument("Radiating-star tables currently require rho_c=1 g/cm^3 and mean molecular weight 0.6");
	}
	return {metadata.index, parameters.centralDensity, metadata.centralBeta, parameters.meanMolecularWeight};
}

RadiatingStarAtmosphere::Parameters atmosphereParameters(RadiatingStarReference::Parameters const& parameters,
	RadiatingStarReference::Metadata const& metadata) {
	RadiatingStarAtmosphere::Parameters result;
	result.eos = eosParameters(parameters, metadata);
	if (!(parameters.atmosphereTolerance > 0 && parameters.atmosphereTolerance <= 1e-6)
		|| !std::isfinite(parameters.atmosphereTolerance)) {
		throw std::invalid_argument("Invalid spherical-background integration tolerance");
	}
	result.opacityLengthScale = metadata.opacity * std::sqrt((metadata.index + 1) / metadata.centralBeta);
	result.cutoffDensityFraction = metadata.cutoffDensityFraction;
	result.tolerance = parameters.atmosphereTolerance;
	result.maximumSteps = 1000000;
	return result;
}

struct CubicValue { Real value{}, derivative{}, secondDerivative{}; };

CubicValue cubic(Real a, Real da, Real b, Real db, Real width, Real t) {
	Real const change = b - a;
	Real const first = width * da;
	Real const second = 3 * change - 2 * first - width * db;
	Real const third = first + width * db - 2 * change;
	return {a + t * (first + t * (second + t * third)),
		(first + t * (2 * second + 3 * t * third)) / width,
		(2 * second + 6 * t * third) / (width * width)};
}

std::size_t interval(std::vector<Real> const& axis, Real coordinate) {
	auto const next = std::upper_bound(axis.begin(), axis.end(), coordinate);
	if (next == axis.begin()) return 0;
	return std::min(std::size_t(next - axis.begin() - 1), axis.size() - 2);
}

bool agrees(Real a, Real b) {
	return std::abs(a - b) <= 64 * std::numeric_limits<Real>::epsilon()
		* std::max({Real(1), std::abs(a), std::abs(b)});
}
}

RadiatingStarReference::Parsed RadiatingStarReference::parse(std::string_view text) {
	std::istringstream input{std::string(text)};
	input.imbue(std::locale::classic());
	std::string magic;
	int version = 0;
	if (!(input >> magic >> version) || magic != "OCTOII_RADIATING_STAR" || version != 1) {
		throw std::invalid_argument("Unsupported radiating-star reference header");
	}
	Parsed result;
	auto& metadata = result.metadata;
	if (!(input >> metadata.index >> metadata.centralBeta >> metadata.opacity >> metadata.cutoffDensityFraction
		>> metadata.spin >> metadata.extent >> metadata.lightSpeed)) {
		throw std::invalid_argument("Incomplete radiating-star reference metadata");
	}
	for (Real item : {metadata.index, metadata.centralBeta, metadata.opacity, metadata.cutoffDensityFraction,
		metadata.spin, metadata.extent, metadata.lightSpeed}) {
		if (!std::isfinite(item)) throw std::invalid_argument("Nonfinite radiating-star reference metadata");
	}
	if (!(metadata.index > 3 && metadata.index < 5 && metadata.centralBeta > 0 && metadata.centralBeta < 1
		&& metadata.opacity > 0 && metadata.cutoffDensityFraction > 0 && metadata.cutoffDensityFraction < 1
		&& metadata.extent > 1 && metadata.lightSpeed > 0)) {
		throw std::invalid_argument("Unphysical radiating-star reference metadata");
	}
	long long coreRows = 0, exteriorRows = 0, angularPoints = 0, fields = 0;
	if (!(input >> coreRows >> exteriorRows >> angularPoints >> fields)
		|| coreRows < 2 || exteriorRows < 2 || angularPoints < 2 || fields != fieldCount
		|| coreRows > 32769 || exteriorRows > 32769 || angularPoints > 2049
		|| (coreRows + exteriorRows) * angularPoints > 1000000) {
		throw std::invalid_argument("Invalid or excessive radiating-star table dimensions");
	}
	auto readAxis = [&](std::vector<Real>& axis, std::size_t count) {
		axis.resize(count);
		for (std::size_t i = 0; i < count; ++i) {
			if (!(input >> axis[i]) || !std::isfinite(axis[i]) || axis[i] < 0
				|| (i && axis[i] <= axis[i - 1])) {
				throw std::invalid_argument("Radiating-star axes must be finite and strictly increasing");
			}
		}
	};
	readAxis(result.mu, angularPoints);
	readAxis(result.core.radius, coreRows);
	readAxis(result.exterior.radius, exteriorRows);
	if (result.mu.front() != 0 || result.mu.back() != 1 || result.core.radius.front() != 0
		|| result.core.radius.back() != result.exterior.radius.front()) {
		throw std::invalid_argument("Radiating-star axes must include the center, equator, pole, and a shared interface");
	}
	auto readGrid = [&](Grid& grid) {
		grid.nodes.resize(grid.radius.size() * result.mu.size());
		for (auto& node : grid.nodes) {
			for (auto& field : node) {
				if (!(input >> field.value >> field.radial >> field.angular >> field.mixed)
					|| !std::isfinite(field.value) || !std::isfinite(field.radial)
					|| !std::isfinite(field.angular) || !std::isfinite(field.mixed)) {
					throw std::invalid_argument("Incomplete or nonfinite radiating-star field data");
				}
			}
		}
	};
	readGrid(result.core);
	readGrid(result.exterior);
	input >> std::ws;
	if (!input.eof()) throw std::invalid_argument("Unexpected trailing radiating-star reference data");
	std::size_t const angles = result.mu.size();
	for (std::size_t k = 0; k < angles; ++k) {
		auto const& inner = result.core.nodes[(result.core.radius.size() - 1) * angles + k];
		auto const& outer = result.exterior.nodes[k];
		for (std::size_t field = 0; field < fieldCount; ++field) {
			if (!agrees(inner[field].value, outer[field].value)) {
				throw std::invalid_argument("Radiating-star core/exterior values disagree at their interface");
			}
		}
	}
	// These are coordinate regularity constraints, not a test of force balance.
	for (auto const* grid : {&result.core, &result.exterior}) {
		for (std::size_t i = 0; i < grid->radius.size(); ++i) {
			auto const& equator = grid->nodes[i * angles];
			for (std::size_t field = 0; field < fieldCount; ++field) {
				bool const odd = field == 3;
				if (!(odd ? agrees(equator[field].value, 0) && agrees(equator[field].radial, 0)
					: agrees(equator[field].angular, 0) && agrees(equator[field].mixed, 0))) {
					throw std::invalid_argument("Radiating-star table violates equatorial parity");
				}
			}
		}
	}
	for (std::size_t k = 0; k < angles; ++k) {
		auto const& center = result.core.nodes[k];
		for (std::size_t field = 0; field < fieldCount; ++field) {
			bool const scalar = field == 0 || field == 1 || field == 4;
			if (!agrees(center[field].angular, 0)
				|| !(scalar ? agrees(center[field].value, result.core.nodes[0][field].value)
					&& agrees(center[field].radial, 0) && agrees(center[field].mixed, 0)
					: agrees(center[field].value, 0))) {
				throw std::invalid_argument("Radiating-star table violates central regularity");
			}
		}
	}
	return result;
}

RadiatingStarReference::RadiatingStarReference(Parameters parameters, std::string_view text)
	: RadiatingStarReference(parameters, parse(text)) {}

RadiatingStarReference::RadiatingStarReference(Parameters parameters, Parsed parsed)
	: parameters_(parameters), metadata_(parsed.metadata),
	  baseline_(atmosphereParameters(parameters, parsed.metadata)),
	  barotrope_({eosParameters(parameters, parsed.metadata), parsed.metadata.cutoffDensityFraction}),
	  mu_(std::move(parsed.mu)), core_(std::move(parsed.core)), exterior_(std::move(parsed.exterior)) {
	hScale_ = value(constants::boltzmann / constants::atomicMassUnit) / parameters.meanMolecularWeight
		* value(eos().centralTemperature());
	pressureScale_ = value(parameters.centralDensity) * hScale_;
	lengthScale_ = value(eos().scaleLength()) / std::sqrt((metadata_.index + 1) / metadata_.centralBeta);
	surfacePotential_ = value(baseline_.sample(baseline_.surfaceRadius()).potential);
	Real const physicalLightSpeed = value(constants::c) / std::sqrt(hScale_);
	if (!agrees(metadata_.lightSpeed, physicalLightSpeed)) {
		throw std::invalid_argument("Radiating-star table light speed does not match its physical scaling");
	}
	Real const interfaceRadius = value(baseline_.transitionRadius()) / lengthScale_;
	if (std::abs(core_.radius.back() / interfaceRadius - 1) > 1e-6) {
		throw std::invalid_argument("Radiating-star table interface does not match its spherical background");
	}
	if (std::abs(exterior_.radius.back() / core_.radius.back() / metadata_.extent - 1) > 1e-12) {
		throw std::invalid_argument("Radiating-star extent metadata disagrees with its radial axis");
	}
}

std::array<RadiatingStarReference::Jet,RadiatingStarReference::fieldCount>
RadiatingStarReference::interpolate(Grid const& grid, Real r, Real mu) const {
	std::size_t const i = interval(grid.radius, r), k = interval(mu_, mu);
	Real const dr = grid.radius[i + 1] - grid.radius[i], dm = mu_[k + 1] - mu_[k];
	Real const tr = (r - grid.radius[i]) / dr, tm = (mu - mu_[k]) / dm;
	std::array<Jet,fieldCount> result{};
	for (std::size_t field = 0; field < fieldCount; ++field) {
		std::array<CubicValue,2> angular{}, radial{};
		for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
			auto const& a = grid.nodes[(i + endpoint) * mu_.size() + k][field];
			auto const& b = grid.nodes[(i + endpoint) * mu_.size() + k + 1][field];
			angular[endpoint] = cubic(a.value, a.angular, b.value, b.angular, dm, tm);
			radial[endpoint] = cubic(a.radial, a.mixed, b.radial, b.mixed, dm, tm);
		}
		auto const value = cubic(angular[0].value, radial[0].value, angular[1].value, radial[1].value, dr, tr);
		auto const derivative = cubic(angular[0].derivative, radial[0].derivative,
			angular[1].derivative, radial[1].derivative, dr, tr);
		result[field] = {value.value, value.derivative, derivative.value, derivative.derivative, value.secondDerivative};
	}
	return result;
}

RadiatingStarReference::State RadiatingStarReference::sample(units::Length cylindricalRadius, units::Length z) const {
	if (!units::finite(cylindricalRadius) || !units::finite(z) || cylindricalRadius < units::Length{}) {
		throw std::invalid_argument("Radiating-star sample coordinates must be finite with R >= 0");
	}
	auto const physicalRadius = units::hypot(cylindricalRadius, z);
	if (physicalRadius > maximumRadius()) throw std::out_of_range("Sample lies outside the radiating-star reference table");
	Real const r = std::min(value(physicalRadius) / lengthScale_, exterior_.radius.back());
	Real const mu = r > 0 ? std::min(std::abs(value(z) / value(physicalRadius)), Real(1)) : 0;
	Real const sine = r > 0 ? value(cylindricalRadius) / value(physicalRadius) : 1;
	Real const symmetry = sine * sine, sign = z < units::Length{} ? -1 : 1;
	bool const inCore = physicalRadius < transitionRadius();
	auto const correction = interpolate(inCore ? core_ : exterior_, r, mu);
	auto const& dh = correction[0];
	auto const& dlogEnergy = correction[1];
	auto const& df = correction[2];
	auto const& q = correction[3];
	auto const& dphi = correction[4];
	auto const& w = correction[5];
	auto const background = baseline_.sample(physicalRadius);
	auto const backgroundGas = barotrope_.atDensity(background.density);
	Real const h0 = background.density > units::Density{} ? value(backgroundGas.integralH) / hScale_
		: (surfacePotential_ - value(background.potential)) / hScale_;
	Real const h0r = background.density > units::Density{}
		? background.densityDerivative / value(backgroundGas.densityDerivative) * lengthScale_ / hScale_
		: -value(background.gravityMagnitude) * lengthScale_ / hScale_;
	auto const h = units::VelocitySquared::from_value((h0 + dh.value) * hScale_);
	// Preserve the exact spherical value when no enthalpy correction is present;
	// an inverse-table round trip must not create fictitious centrifugal support.
	auto const gas = dh.value == 0 && background.density > units::Density{}
		? backgroundGas : barotrope_.atIntegralH(h);
	Real const hRadial = h0r + dh.radial;
	Real const phiRadial = value(background.gravityMagnitude) * lengthScale_ / hScale_ + dphi.radial;
	Real const e0 = value(background.radiationEnergy) / pressureScale_;
	Real const e0r = background.radiationEnergyDerivative * lengthScale_ / pressureScale_;
	Real const f0 = background.fluxFactor;
	Real const f0r = (background.radiationFluxDerivative / value(constants::c)
		- f0 * background.radiationEnergyDerivative) * lengthScale_ / value(background.radiationEnergy);
	Real energy = 0, energyRadial = 0, energyMu = 0, deltaEnergy = 0;
	Real temperatureH = 0;
	if (gas.density > units::Density{}) {
		if (gas.density <= barotrope_.cutoffDensity()) {
			temperatureH = (barotrope_.envelopeGamma() - 1) * value(gas.temperature) / value(gas.integralH);
		} else {
			auto const thermodynamics = eos().atDensity(gas.density);
			Real const slope = (1 + 1 / metadata_.index - thermodynamics.beta) / (4 - 3 * thermodynamics.beta);
			temperatureH = value(gas.temperature) / value(gas.density) * slope * value(gas.densityDerivative);
		}
	}
	if (inCore) {
		if (!(gas.density > units::Density{})) throw std::runtime_error("Radiating-star opaque core contains vacuum");
		Real const logTemperatureRatio = std::log(value(gas.temperature) / value(background.temperature));
		deltaEnergy = e0 * std::expm1(4 * logTemperatureRatio);
		energy = e0 + deltaEnergy;
		Real const energyH = 4 * energy / value(gas.temperature) * temperatureH * hScale_;
		auto const backgroundEos = eos().atDensity(background.density);
		Real const backgroundSlope = (1 + 1 / metadata_.index - backgroundEos.beta) / (4 - 3 * backgroundEos.beta);
		Real const backgroundEnergyH = 4 * e0 / value(background.density)
			* backgroundSlope * value(backgroundGas.densityDerivative) * hScale_;
		energyRadial = e0r + (energyH - backgroundEnergyH) * h0r + energyH * dh.radial;
		energyMu = energyH * dh.angular;
	} else {
		deltaEnergy = e0 * std::expm1(dlogEnergy.value);
		energy = e0 + deltaEnergy;
		energyRadial = e0r * std::exp(dlogEnergy.value) + energy * dlogEnergy.radial;
		energyMu = energy * dlogEnergy.angular;
	}
	Real const radialFactor = f0 + df.value;
	Real const radialFactorRadial = f0r + df.radial;
	Real const radialFlux = energy * radialFactor;
	Real const polarFlux = energy * q.value;
	Real const deltaRadialFlux = deltaEnergy * f0 + energy * df.value;
	Real const opacity = background.opacityCgs * value(parameters_.centralDensity) * lengthScale_;
	Real const fluxFactorSquared = radialFactor * radialFactor + symmetry * q.value * q.value
		+ symmetry * std::pow(w.value / energy, 2);
	if (!(energy > 0) || !std::isfinite(energy) || !std::isfinite(fluxFactorSquared) || fluxFactorSquared > 1) {
		throw std::runtime_error("Radiating-star interpolant left the physical M1 flux cone");
	}
	Real omegaSquared = 0;
	if (gas.density > units::Density{}) {
		if (r > 0) {
			omegaSquared = (dh.radial + dphi.radial - mu * (dh.angular + dphi.angular) / r
				- opacity * (deltaRadialFlux - mu * polarFlux)) / r;
		} else {
			omegaSquared = dh.radialSecond + dphi.radialSecond
				- opacity * (deltaEnergy * f0r + energy * df.radial);
		}
		if (!std::isfinite(omegaSquared) || omegaSquared < 0) {
			throw std::runtime_error("Radiating-star reference requires negative centrifugal support");
		}
	}
	Real heating = 0;
	if (inCore) {
		if (r > 0) {
			heating = energyRadial * radialFactor + energy * radialFactorRadial
				+ (2 * radialFlux + symmetry * (energyMu * q.value + energy * q.angular) - 2 * mu * polarFlux) / r;
		} else {
			heating = energy * (3 * radialFactorRadial + q.mixed);
		}
	}
	auto cylindricalGradient = [&](Real radial, Real angular) {
		if (r == 0) return std::array<Real,2>{0, 0};
		return std::array<Real,2>{sine * (radial - mu * angular / r),
			sign * (mu * radial + symmetry * angular / r)};
	};
	auto const hg = cylindricalGradient(hRadial, dh.angular);
	auto const pg = cylindricalGradient(phiRadial, dphi.angular);
	auto const eg = cylindricalGradient(energyRadial, energyMu);
	Real const fluxScale = value(constants::c) * pressureScale_;
	State result;
	result.density = gas.density;
	result.temperature = gas.temperature;
	result.gasPressure = gas.gasPressure;
	result.integralH = h;
	result.radiationEnergy = units::EnergyDensity::from_value(energy * pressureScale_);
	result.potential = background.potential + units::VelocitySquared::from_value(dphi.value * hScale_);
	result.radialFlux = units::EnergyFlux::from_value(radialFlux * fluxScale);
	result.thetaFlux = units::EnergyFlux::from_value(-sign * sine * polarFlux * fluxScale);
	result.cylindricalFlux = units::EnergyFlux::from_value(sine * (radialFlux - mu * polarFlux) * fluxScale);
	result.verticalFlux = units::EnergyFlux::from_value(sign * (mu * radialFlux + symmetry * polarFlux) * fluxScale);
	result.azimuthalFlux = units::EnergyFlux::from_value(sine * w.value * fluxScale);
	result.angularVelocity = units::InverseTime::from_value(std::copysign(std::sqrt(omegaSquared * hScale_) / lengthScale_, metadata_.spin));
	for (std::size_t d = 0; d < 2; ++d) {
		result.integralHGradient[d] = units::Acceleration::from_value(hg[d] * hScale_ / lengthScale_);
		result.potentialGradient[d] = units::Acceleration::from_value(pg[d] * hScale_ / lengthScale_);
		result.gravity[d] = -result.potentialGradient[d];
		result.densityGradient[d] = gas.densityDerivative * result.integralHGradient[d];
		result.temperatureGradient[d] = units::Quantity<-1,0,0,1>::from_value(temperatureH * value(result.integralHGradient[d]));
		result.gasPressureGradient[d] = gas.density * result.integralHGradient[d];
		result.radiationEnergyGradient[d] = units::Quantity<-2,1,-2>::from_value(eg[d] * pressureScale_ / lengthScale_);
	}
	result.mu = sign * mu;
	result.radialFluxFactor = radialFactor;
	result.polarFluxFactor = sign * q.value;
	result.radialFactorRadialDerivative = units::Quantity<-1,0,0>::from_value(radialFactorRadial / lengthScale_);
	result.polarFactorRadialDerivative = units::Quantity<-1,0,0>::from_value(sign * q.radial / lengthScale_);
	result.radialFactorMuDerivative = sign * df.angular;
	result.polarFactorMuDerivative = q.angular;
	result.azimuthalFluxCoefficient = w.value;
	result.azimuthalCoefficientRadialDerivative = units::Quantity<-1,0,0>::from_value(w.radial / lengthScale_);
	result.azimuthalCoefficientMuDerivative = sign * w.angular;
	result.energyRadialDerivative = units::Quantity<-2,1,-2>::from_value(energyRadial * pressureScale_ / lengthScale_);
	result.energyMuDerivative = units::EnergyDensity::from_value(sign * energyMu * pressureScale_);
	result.opacityCgs = background.opacityCgs;
	result.photonHeatingCgs = heating * fluxScale / lengthScale_;
	return result;
}

units::Length RadiatingStarReference::lengthScale() const { return units::Length::from_value(lengthScale_); }
units::Pressure RadiatingStarReference::pressureScale() const { return units::Pressure::from_value(pressureScale_); }
units::VelocitySquared RadiatingStarReference::velocitySquaredScale() const { return units::VelocitySquared::from_value(hScale_); }
units::Length RadiatingStarReference::maximumRadius() const { return units::Length::from_value(exterior_.radius.back() * lengthScale_); }
units::Length RadiatingStarReference::transitionRadius() const { return units::Length::from_value(core_.radius.back() * lengthScale_); }

} // namespace octotigerII::problems
