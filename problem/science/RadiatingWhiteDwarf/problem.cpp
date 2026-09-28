/** @file
 * @brief Rotating cold-degenerate WD with thermal ions and full M1 radiation.
 */
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/whiteDwarfStructure.hpp"
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

namespace octotigerII::radiatingWhiteDwarf {
namespace {
using Model=problems::WhiteDwarfStructure;

std::shared_ptr<Model const> model(Config const& config) {
	Model::Parameters p;
	p.centralDensity=config.star.centralDensity;
	p.meanMassPerElectron=config.hydro.meanMassPerElectron;
	p.meanMassPerIon=config.hydro.meanMolecularWeight;
	p.thermalPressureFraction=config.whiteDwarf.thermalPressureFraction;
	p.spinFractionOfSphericalBreakup=config.whiteDwarf.rotationFraction;
	p.radialCells=config.whiteDwarf.radialCells;
	p.angularPoints=config.whiteDwarf.angularPoints;
	p.maxMultipole=config.whiteDwarf.multipoles;
	p.tolerance=config.whiteDwarf.structureTolerance;
	using Key=std::array<Real,9>;
	Key const key{units::value(p.centralDensity),p.meanMassPerElectron,p.meanMassPerIon,
		p.thermalPressureFraction,p.spinFractionOfSphericalBreakup,Real(p.radialCells),
		Real(p.angularPoints),Real(p.maxMultipole),p.tolerance};
#ifdef OCTOTIGERII_WITH_HPX
	static hpx::shared_mutex mutex;
#else
	static std::shared_mutex mutex;
#endif
	static std::map<Key,std::shared_ptr<Model const>> cache;
	{
		std::shared_lock lock(mutex);
		if(auto found=cache.find(key);found!=cache.end())return found->second;
	}
	auto value=std::make_shared<Model const>(p);
	std::unique_lock lock(mutex);
	if(cache.size()>=8)cache.clear();
	return cache.try_emplace(key,std::move(value)).first->second;
}

verification::ExactState state(Config const& config,Model const& star,
	mesh::PhysicalCoordinates position,Real probeScale=1) {
	for(int axis=0;axis<ndim;++axis)position[axis]-=config.star.center[axis];
	auto const R=units::hypot(position[0],position[1]);
	auto const profile=star.sample(R/probeScale,position[2]/probeScale);
	auto const floorDensity=config.star.atmosphereFraction*config.star.centralDensity;
	verification::ExactState result;
	result.hydro.density()=std::max(profile.density,floorDensity);
	if(profile.density>=floorDensity){
		result.hydro.pressure()=profile.pressure;
	}else{
		hydro::WhiteDwarfEos const cold(config.hydro.meanMassPerElectron);
		result.hydro.pressure()=(1+config.whiteDwarf.thermalPressureFraction)*cold.pressure(floorDensity);
	}
	Real const materialFraction=Real(profile.density/result.hydro.density());
	result.hydro.velocity(0)=-materialFraction*star.angularVelocity()*position[1];
	result.hydro.velocity(1)=materialFraction*star.angularVelocity()*position[0];

	auto const coldPressure=hydro::WhiteDwarfEos(config.hydro.meanMassPerElectron).pressure(result.hydro.density());
	auto const ionPressure=result.hydro.pressure()-coldPressure;
	auto const temperature=ionPressure*(config.hydro.meanMolecularWeight*constants::atomicMassUnit)/
		(result.hydro.density()*constants::boltzmann);
	auto const E0=units::EnergyDensity::from_value(units::value(constants::radiation)*
		std::pow(units::value(temperature),4));
	std::array<units::EnergyFlux,ndim> comovingFlux{};
	Real const nx=R>units::Length{}?Real(position[0]/R):0;
	Real const ny=R>units::Length{}?Real(position[1]/R):0;
	if(profile.density>floorDensity){
		auto const diffusion=star.diffusionFlux(profile,config.radiation.opacity+config.radiation.scatteringOpacity);
		comovingFlux[0]=nx*diffusion[0];
		comovingFlux[1]=ny*diffusion[0];
		comovingFlux[2]=diffusion[1];
	}
	units::EnergyFlux norm{};
	for(auto flux:comovingFlux)norm=units::hypot(norm,flux);
	auto const maximumFlux=(1-32*epsilonR)*constants::c*E0;
	if(norm>maximumFlux){
		for(auto& flux:comovingFlux)flux*=Real(maximumFlux/norm);
		norm=maximumFlux;
	}
	Real const f=Real(norm/(constants::c*E0));
	Real const f2=std::min(Real(1),f*f);
	auto const transversePressure=E0*((1-f2)/(1+std::sqrt(4-3*f2)));
	Real beta2=0;
	for(int axis=0;axis<ndim;++axis){
		Real const beta=Real(result.hydro.velocity(axis)/constants::c);
		beta2+=beta*beta;
	}
	if(!(beta2<1))throw std::invalid_argument("White-dwarf rotation must be subluminal");
	Real const boost2=1/(1-beta2),boost=std::sqrt(boost2);
	result.radiation.energy()=boost2*(E0+beta2*transversePressure);
	for(int axis=0;axis<ndim;++axis)
		result.radiation.radiativeFlux(axis)=boost*comovingFlux[axis]+
			boost2*(E0+transversePressure)*result.hydro.velocity(axis);
	result.gravity.potential()=profile.potential;
	result.gravity.acceleration(0)=-nx*profile.potentialGradient[0];
	result.gravity.acceleration(1)=-ny*profile.potentialGradient[0];
	result.gravity.acceleration(2)=-profile.potentialGradient[1];
	return result;
}
}

ProblemBoundary problemBoundary(Config const& config) {
	auto star=model(config);
	return [config,star](auto const& position,units::Time){return state(config,*star,position);};
}

verification::Reference problemReference(Config const& config) {
	verification::Reference result;
	result.name="Initial rotating degenerate white dwarf with LTE M1 radiation";
	result.evaluate=radiatingWhiteDwarf::problemBoundary(config);
	return result;
}

void problemDefaults(Config& config) {
	config.star.centralDensity=units::Density::from_value(1e9);
	config.star.atmosphereFraction=1e-12;
	config.hydro.eos="white-dwarf";
	config.hydro.gamma=Real(5)/3;
	config.hydro.meanMolecularWeight=14;
	config.hydro.meanMassPerElectron=2;
	config.hydro.dualEnergy.enabled=true;
	config.frame.omega={};
	config.mesh.cells=8;
	config.mesh.level=0;
	config.mesh.boundary=finiteVolume::BoundaryConditions::uniform(finiteVolume::BoundaryCondition::Outflow);
	config.amr.enabled=true;
	config.amr.minLevel=0;
	config.amr.maxLevel=4;
	config.amr.refineDensity=.001*config.star.centralDensity;
	config.amr.shadowTolerance=0;
	config.radiation.lightSpeedRatio=1;
	config.radiation.opacity=.2;
	config.radiation.closedBoundary=true;
	config.radiation.enabled=false;
}

void validateProblem(Config const& config) {
	if(config.hydro.eos!="white-dwarf" || std::abs(config.hydro.gamma-Real(5)/3)>1e-12
		|| !config.hydro.dualEnergy.enabled || config.radiation.enabled)
		throw std::invalid_argument("The rotating WD requires cold-electron/ideal-ion EOS, gamma=5/3, dual energy, and explicit M1 moments");
	if(config.radiation.lightSpeedRatio!=1 || !config.radiation.closedBoundary ||
		config.radiation.opacityModel != "constant" ||
		!(config.radiation.opacity + config.radiation.scatteringOpacity > 0))
		throw std::invalid_argument("The rotating WD requires full light speed, positive constant flux opacity, and a closed radiation boundary");
	if(config.frame.omega!=units::InverseTime{} ||
		!config.mesh.boundary.all(finiteVolume::BoundaryCondition::Outflow))
		throw std::invalid_argument("The rotating WD uses an inertial grid and outflow gas boundaries");
	if(!(config.star.atmosphereFraction>0 && config.star.atmosphereFraction<1e-3))
		throw std::invalid_argument("Invalid WD numerical atmosphere fraction");
	for(auto center:config.star.center)if(!units::finite(center))
		throw std::invalid_argument("Nonfinite WD center");
	auto const star=model(config);
	for(int axis=0;axis<ndim;++axis){
		auto const radius=axis==2?star->polarRadius():star->equatorialRadius();
		if(config.star.center[axis]-radius<=config.mesh.lower ||
			config.star.center[axis]+radius>=config.mesh.upper)
			throw std::invalid_argument("The complete rotating WD must fit inside the domain");
	}
	if(config.star.atmosphereFraction*config.star.centralDensity<units::Density::from_value(1e-14))
		throw std::invalid_argument("WD numerical atmosphere falls below the hydro density floor");
}

void initializeProblem(Snapshot& data,Config const& config,bool refinementProbe) {
	auto star=model(config);
	Real const scale=refinementProbe?Real(initialFeatureWidth(star->scaleLength(),data.cellWidth)/star->scaleLength()):1;
	hydro::HydroSystem const gas(config.hydro);
	data.layout.forEachInterior([&](auto const& cell,std::size_t i){
		auto const position=data.layout.cellCenter(data.lower,data.cellWidth,cell);
		auto const value=state(config,*star,position,scale);
		data.hydro.values()[i]=gas.conservedState(value.hydro);
		data.radiation.values()[i]=value.radiation;
	});
}
} // namespace octotigerII::radiatingWhiteDwarf
