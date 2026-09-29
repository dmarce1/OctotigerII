#include "testSupport.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/radiatingStarStructure.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <cmath>
#include <iostream>
#include <numbers>

using namespace octotigerII;
namespace {
Real value(auto q){return units::value(q);}
Config configuration(int level=1){
	auto c=test::parseConfig({"--mesh.cells=4","--mesh.level="+std::to_string(level),
		"--amr.enabled=off","--output.enabled=off","--runtime.stopTime=0"});
	problems::RadiatingStarEos const eos({c.star.polytropicIndex,c.star.centralDensity,c.radiatingStar.centralGasFraction,c.hydro.meanMolecularWeight});
	c.mesh.lower=-8.0*eos.scaleLength();c.mesh.upper=8.0*eos.scaleLength();
	c.amr.minLevel=c.amr.maxLevel=level;c.validate();return c;
}
mesh::PhysicalCoordinates axisPoint(units::Length r){mesh::PhysicalCoordinates x{};x[0]=r;return x;}
}

TEST(RadiatingSphere, InitializesPhysicalMomentsAndFrozenMaterialWithConsistentGasEos){
	auto const c=configuration();EXPECT_EQ(c.hydro.gamma,Real(5)/3);EXPECT_EQ(c.hydro.meanMolecularWeight,.6);
	EXPECT_TRUE(c.hydroEnabled());EXPECT_TRUE(c.radiationEnabled());EXPECT_TRUE(c.gravityEnabled());
	EXPECT_FALSE(c.radiation.enabled);EXPECT_TRUE(problemHasRadiationMaterial(c));
	auto const boundary=problemBoundary(c);auto const material=problemRadiationMaterial(c);
	problems::RadiatingStarEos const eos({c.star.polytropicIndex,c.star.centralDensity,c.radiatingStar.centralGasFraction,c.hydro.meanMolecularWeight});
	for(Real x:{Real(0),Real(.5),Real(2),Real(3),Real(5),Real(7)}){
		auto const position=axisPoint(x*eos.scaleLength());auto const a=boundary(position,{}),b=boundary(position,units::Time::from_value(100));
		test::expectStateNear(a.hydro,b.hydro,0);test::expectStateNear(a.radiation,b.radiation,0);
		auto const opacity=material(position,{}),later=material(position,units::Time::from_value(100));
		EXPECT_EQ(opacity.opacity,later.opacity);EXPECT_EQ(opacity.photonPower,later.photonPower);
		EXPECT_GE(value(opacity.opacity),0);EXPECT_GE(value(opacity.photonPower),0);
		for(int d=0;d<ndim;++d) {EXPECT_EQ(a.hydro.velocity(d),units::Velocity{}); }
		if(x<3.5){
			auto const temperature=a.hydro.pressure()*c.hydro.meanMolecularWeight*constants::atomicMassUnit/(a.hydro.density()*constants::boltzmann);
			EXPECT_NEAR(Real(a.radiation.energy()/(constants::radiation*boost::units::pow<4>(temperature))),1,3e-13);
		}else if(x>4){EXPECT_EQ(value(opacity.opacity),0);EXPECT_EQ(value(opacity.photonPower),0);}
	}
	hydro::HydroSystem const gas(c.hydro);auto const initial=initialSnapshot(c,{0,{}});
	initial.layout.forEachInterior([&](auto const& cell,std::size_t i){
		auto const exact=boundary(initial.layout.cellCenter(initial.lower,initial.cellWidth,cell),{});
		test::expectStateNear(initial.hydro.values()[i],gas.conservedState(exact.hydro));
		test::expectStateNear(initial.radiation.values()[i],exact.radiation);
	});
}

TEST(RadiatingSphere, ProblemCallbacksRecoverMechanicalAndLuminosityBalances){
	auto const c=configuration();auto const boundary=problemBoundary(c);auto const material=problemRadiationMaterial(c);
	problems::RadiatingStarEos const eos({c.star.polytropicIndex,c.star.centralDensity,c.radiatingStar.centralGasFraction,c.hydro.meanMolecularWeight});
	// Differentiate inside a reference interpolation interval; a wider stencil
	// averages the small, resolved variation of its prescribed photon source.
	Real const alpha=value(eos.scaleLength()),h=1e-7*alpha;
	for(Real x:{Real(.3),Real(1),Real(2),Real(3.4),Real(4),Real(5)}){
		Real const r=x*alpha;auto const q=boundary(axisPoint(units::Length::from_value(r)),{});
		auto const a=boundary(axisPoint(units::Length::from_value(r-h)),{}),b=boundary(axisPoint(units::Length::from_value(r+h)),{});
		auto const m=material(axisPoint(units::Length::from_value(r)),{});
		Real const force=value(m.opacity)*value(q.hydro.density())*value(q.radiation.radiativeFlux(0))/value(constants::c);
		Real const gravity=-value(q.hydro.density())*value(q.gravity.acceleration(0));
		Real const derivative=value(b.hydro.pressure()-a.hydro.pressure())/(2*h);
		EXPECT_NEAR((derivative+gravity-force)/std::max({std::abs(derivative),gravity,std::abs(force)}),0,5e-6);
		Real const divFlux=value(b.radiation.radiativeFlux(0)-a.radiation.radiativeFlux(0))/(2*h)+2*value(q.radiation.radiativeFlux(0))/r;
		Real const heat=value(m.photonPower),scale=std::max(std::abs(heat),value(q.radiation.radiativeFlux(0))/r);
		EXPECT_NEAR((divFlux-heat)/scale,0,3e-5) << x;
	}
}

TEST(RadiatingSphere, RotatingCoordinatesKeepThePrescribedMaterialStationary){
	auto c=configuration();
	c.frame.omega=units::InverseTime::from_value(1e-5);
	EXPECT_NO_THROW(c.validate());
	c.star.center[0]=units::Length::from_value(1);
	EXPECT_THROW(c.validate(),std::invalid_argument);
	c.frame.omega=units::InverseTime{};
	EXPECT_NO_THROW(c.validate());
}

TEST(RadiatingSphere, EvolvedDriftDecreasesWithMeshSpacingAndAccountsForSourceEnergy){
	struct Result{Real densityDrift{},rmsMach{},energyResidual{},massResidual{},radiationDrift{};int steps{};};
	units::Time interval{};
	std::array<Result,2> result{};
	for(int resolution=0;resolution<2;++resolution){
		auto const c=configuration(resolution+1);Runtime runtime(c);runtime.solveGravity();
		auto const initial=runtime.snapshots();auto const before=diagnose(initial,c);
		if(resolution==0)interval=8.0*runtime.stableTimestep();
		units::Time elapsed{};
		while(elapsed<interval){auto const dt=std::min(runtime.stableTimestep(),interval-elapsed);runtime.advanceCoupled(dt);elapsed+=dt;++result[resolution].steps;}
		auto const snapshots=runtime.snapshots();auto const after=diagnose(snapshots,c);auto const boundary=runtime.boundaryTransport();
		auto const out=boundary.outward.gasEnergy-boundary.inward.gasEnergy+boundary.outward.potentialEnergy-boundary.inward.potentialEnergy
			+(boundary.outward.radiationEnergy-boundary.inward.radiationEnergy)/c.radiation.lightSpeedRatio;
		auto const source=runtime.radiationSourceEnergy()/c.radiation.lightSpeedRatio;
		auto& r=result[resolution];r.energyResidual=Real((after.rslaTotalEnergy+out-before.rslaTotalEnergy-source)/before.rslaTotalEnergyNorm);
		r.massResidual=Real((after.mass+boundary.outward.mass-boundary.inward.mass-before.mass)/before.mass);
		r.rmsMach=std::sqrt(Real(2.0*after.kineticEnergy/(c.hydro.gamma*(c.hydro.gamma-1)*after.thermalEnergy)));
		r.radiationDrift=Real((after.radiationEnergy-before.radiationEnergy)/before.radiationEnergy);
		long double densityChange=0,densityNorm=0;
		ASSERT_EQ(initial.size(),snapshots.size());
		for(std::size_t b=0;b<initial.size();++b) {for(std::size_t i=0;i<initial[b].hydro.values().size();++i){
			Real const volume=std::pow(value(initial[b].cellWidth),ndim);
			densityChange+=std::abs(value(snapshots[b].hydro.values()[i].density()-initial[b].hydro.values()[i].density()))*volume;
			densityNorm+=value(initial[b].hydro.values()[i].density())*volume;
		}
		}
		r.densityDrift=densityChange/densityNorm;
		std::cout<<"radiating sphere cells="<<(8<<resolution)<<" steps="<<r.steps<<" time="<<value(interval)
			<<" densityL1="<<r.densityDrift<<" rmsMach="<<r.rmsMach<<" radiationDrift="<<r.radiationDrift
			<<" energyBudget="<<r.energyResidual<<" massBudget="<<r.massResidual<<'\n';
		EXPECT_NEAR(r.energyResidual,0,2e-12);EXPECT_NEAR(r.massResidual,0,2e-12);
		EXPECT_GT(source,units::Energy{});EXPECT_GE(after.minimumRadiationEnergy,units::EnergyDensity{});
	}
	EXPECT_LT(result[1].densityDrift,result[0].densityDrift);
	EXPECT_LT(result[1].rmsMach,result[0].rmsMach);
}
