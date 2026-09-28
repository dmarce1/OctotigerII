/** @file
 * @brief Spherical full-M1 equilibrium control for the rotating radiative star.
 */
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/radiatingStarAtmosphere.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/synchronization/shared_mutex.hpp>
#endif

namespace octotigerII::radiating_sphere {
namespace {
using Model = problems::RadiatingStarAtmosphere;

std::shared_ptr<Model const> model(Config const& c) {
	Model::Parameters p;
	p.eos = {c.star.polytropicIndex, c.star.centralDensity,
		c.radiatingStar.centralGasFraction, c.hydro.meanMolecularWeight};
	problems::RadiatingStarEos const eos(p.eos);
	p.opacityLengthScale = c.radiation.opacity * units::value(c.star.centralDensity) * units::value(eos.scaleLength());
	p.cutoffDensityFraction = c.radiatingStar.opacityCutoffFraction;
	p.tolerance = c.radiatingStar.structureTolerance;
	p.maximumRadialStep = Real(1) / c.radiatingStar.radialCells;
	std::array<Real,8> const key{p.eos.index, units::value(p.eos.centralDensity), p.eos.centralBeta,
		p.eos.meanMolecularWeight, p.opacityLengthScale, p.cutoffDensityFraction, p.tolerance, p.maximumRadialStep};
#ifdef OCTOTIGERII_WITH_HPX
	static hpx::shared_mutex mutex;
#else
	static std::shared_mutex mutex;
#endif
	static std::map<decltype(key), std::shared_ptr<Model const>> cache;
	{
		std::shared_lock lock(mutex);
		if (auto it=cache.find(key); it!=cache.end()) return it->second;
	}
	auto value=std::make_shared<Model const>(p);
	std::unique_lock lock(mutex);
	// Active evaluators own their model independently of this small reuse cache.
	if (cache.size()>=8) cache.clear();
	return cache.try_emplace(key,std::move(value)).first->second;
}

units::Length radius(Config const& c, mesh::PhysicalCoordinates const& position) {
	units::Length result{};
	for (int d=0; d<ndim; ++d) result=units::hypot(result,position[d]-c.star.center[d]);
	return result;
}

verification::ExactState state(Config const& c, Model const& star, mesh::PhysicalCoordinates const& position,
	Real probeScale=1) {
	auto const r=radius(c,position);
	auto const profile=star.sample(r/probeScale);
	verification::ExactState result;
	result.hydro.density()=std::max(profile.density,c.star.atmosphereFraction*c.star.centralDensity);
	// Continue the actual transparent-envelope polytrope to the density floor;
	// its normalization differs from the central gas-pressure coefficient.
	auto const pressureFloor=star.sample(star.transitionRadius()).gasPressure
		*std::pow(c.star.atmosphereFraction/c.radiatingStar.opacityCutoffFraction,star.diagnostics().envelopeGamma);
	result.hydro.pressure()=std::max(profile.gasPressure,pressureFloor);
	result.radiation.energy()=profile.radiationEnergy;
	result.gravity.potential()=profile.potential;
	if(r>units::Length{}) for(int d=0;d<ndim;++d) {
		Real const direction=Real((position[d]-c.star.center[d])/r);
		result.radiation.radiativeFlux(d)=direction*profile.radiationFlux;
		result.gravity.acceleration(d)=-direction*profile.gravityMagnitude;
	}
	return result;
}
}

ProblemBoundary problemBoundary(Config const& c) {
	auto star=model(c);
	return [c,star](auto const& position, units::Time){ return state(c,*star,position); };
}

ProblemRadiationMaterial problemRadiationMaterial(Config const& c) {
	auto star=model(c);
	return [c,star](auto const& position, units::Time) {
		auto const value=star->sample(radius(c,position));
		return RadiationMaterial{units::Quantity<2,-1,0>::from_value(value.opacityCgs),
			units::Quantity<-1,1,-3>::from_value(value.heatingCgs)};
	};
}

verification::Reference problemReference(Config const& c) {
	verification::Reference result;
	result.name="Spherical full-M1 radiative equilibrium (transparent envelope; numerical gas floor)";
	result.evaluate=radiating_sphere::problemBoundary(c);
	return result;
}

void problemDefaults(Config& c) {
	c.star.centralDensity=units::Density::from_value(1);
	c.star.polytropicIndex=3.5;
	c.star.atmosphereFraction=1e-12;
	c.radiatingStar.rotationFraction=0;
	c.hydro.gamma=Real(5)/3;
	c.hydro.meanMolecularWeight=0.6;
	c.mesh.cells=8;
	c.mesh.level=0;
	c.mesh.boundary=finiteVolume::BoundaryConditions::uniform(finiteVolume::BoundaryCondition::Outflow);
	c.amr.enabled=true;
	c.amr.minLevel=0;
	c.amr.maxLevel=5;
	c.amr.refineDensity=0.001*c.star.centralDensity;
	c.amr.shadowTolerance=0;
	c.radiation.lightSpeedRatio=1;
	// RADIATION in the manifest initializes explicit M1 moments. Do not use the
	// separate "add LTE radiation to a gas-only problem" initialization switch.
	c.radiation.enabled=false;
}

void validateProblem(Config const& c) {
	if (c.radiatingStar.rotationFraction!=0)
		throw std::invalid_argument("radiating-sphere is the nonrotating control; use the rotating reference for nonzero spin");
	if (std::abs(c.hydro.gamma-Real(5)/3)>1e-12 || c.radiation.enabled)
		throw std::invalid_argument("radiating-sphere requires gas gamma=5/3 and explicit M1 moments (radiation.enabled must remain off)");
	if (!(c.star.atmosphereFraction>0 && c.star.atmosphereFraction<1e-3)
		|| !std::isfinite(c.star.atmosphereFraction) || c.radiatingStar.radialCells<64 || c.radiatingStar.radialCells>4096
		|| !(c.radiatingStar.structureTolerance>0 && c.radiatingStar.structureTolerance<=1e-6))
		throw std::invalid_argument("Invalid radiating-star atmosphere fraction, reference resolution, or tolerance");
	if (!c.mesh.boundary.all(finiteVolume::BoundaryCondition::Outflow))
		throw std::invalid_argument("The isolated radiating sphere requires outflow boundaries");
	for (auto center:c.star.center) if(!units::finite(center)) throw std::invalid_argument("Nonfinite stellar center");
	if(c.frame.omega!=units::InverseTime{} && (c.star.center[0]!=units::Length{} || c.star.center[1]!=units::Length{}))
		throw std::invalid_argument("A stationary radiating sphere on a rotating grid must be centered on the grid rotation axis; its prescribed source and opacity use grid coordinates");
	auto star=model(c);
	for (int d=0;d<ndim;++d) if(c.star.center[d]-star->surfaceRadius()<=c.mesh.lower || c.star.center[d]+star->surfaceRadius()>=c.mesh.upper)
		throw std::invalid_argument("The complete radiating star must fit inside the domain");
	if (c.star.atmosphereFraction*c.star.centralDensity<units::Density::from_value(1e-14))
		throw std::invalid_argument("Radiating-star numerical atmosphere falls below the hydro density floor");
}

void initializeProblem(Snapshot& data, Config const& c, bool refinementProbe) {
	auto star=model(c);
	Real const scale=refinementProbe?Real(initialFeatureWidth(star->eos().scaleLength(),data.cellWidth)/star->eos().scaleLength()):1;
	hydro::HydroSystem const gas(c);
	data.layout.forEachInterior([&](auto const& cell,std::size_t i) {
		auto const position=data.layout.cellCenter(data.lower,data.cellWidth,cell);
		auto const value=state(c,*star,position,scale);
		data.hydro.values()[i]=gas.conservedState(value.hydro);
		data.radiation.values()[i]=value.radiation;
	});
}
} // namespace octotigerII::radiating_sphere
