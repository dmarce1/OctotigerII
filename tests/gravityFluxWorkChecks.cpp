#include "testSupport.hpp"
#include "octotigerII/gravity/fluxWork.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_here.hpp>
#endif

using namespace octotigerII;

namespace {
class GravityFluxWork : public ::testing::Test {
protected:
	mesh::MeshLayout layout{2};
	storage::Layout cells{{layout.cellCount(), 1}, 1}, faces{{allFaceCount(2), 1}, 1};
	storage::PartitionSet store{
#ifdef OCTOTIGERII_WITH_HPX
		{hpx::find_here()}
#else
		{0}
#endif
	};
	storage::ColumnFields<gravity::State> field{cells, store, "test.gravity", false, 2};
	storage::Field<units::MassFlux> mass{faces, store, 2, "test.massFlux"};
	Subgrid block;
	units::Time dt = units::Time::from_value(0.1);

	void SetUp() override {
		block.layout = layout;
		block.cellWidth = units::Length::from_value(0.25);
		block.interior = cells.ranges()[0];
		block.massFlux = faces.ranges()[0];
	}

	void setPotential(Real inside, Real neighbor, Real acceleration = 0, unsigned bank = 0) {
		auto handle = field.handle();
		for (auto range : cells.ranges()) {
			auto output = handle.output(range, bank);
			for (std::size_t i = 0; i < range.count; ++i) {
				gravity::State state;
				state.potential() = units::VelocitySquared::from_value(range.offset == 0 ? inside : neighbor);
				state.acceleration(0) = units::Acceleration::from_value(acceleration);
				output.put(i, state);
			}
			handle.commit(range, bank, output);
		}
	}

	void setFlux(mesh::Coordinates const& coordinate, Real own, Real fine, unsigned bank = 0) {
		auto output = mass.handle().output(block.massFlux, bank);
		std::fill_n(output.data(), output.size(), units::MassFlux{});
		output.data()[allFaceIndex(layout, 0, coordinate)] = units::MassFlux::from_value(own);
		mass.handle().commit(block.massFlux, bank, output);
		auto neighbor = mass.handle().output(faces.ranges()[1], bank);
		neighbor.data()[0] = units::MassFlux::from_value(fine);
		mass.handle().commit(faces.ranges()[1], bank, neighbor);
	}
};
}

TEST_F(GravityFluxWork, InternalFaceUsesMeanEndpointPotentialAndSelectedFluxBank) {
	setPotential(0, 0, 0, 0);
	setPotential(0, 0, 0, 1);
	mesh::Coordinates right{};
	right[0] = 1;
	for (unsigned bank = 0; bank < 2; ++bank) {
		auto values = field.handle().output(block.interior, bank);
		auto state = values.at(layout.index(right));
		state.potential() = units::VelocitySquared::from_value(bank ? 4 : 2);
		values.put(layout.index(right), state);
		field.handle().commit(block.interior, bank, values);
	}
	setFlux(right, 3, 0, 1);
	auto const result = gravity::fluxWork(block, {}, field.handle(), 0, field.handle(), 1, mass.handle(), 1, dt);
	for (std::size_t i = 0; i < result.work.size(); ++i)
		EXPECT_NEAR(units::value(result.work[i]), i == 0 || i == layout.index(right) ? -1.8 : 0, 1e-14);
	EXPECT_EQ(result.boundary.inward.potentialEnergy, units::Energy{});
	EXPECT_EQ(result.boundary.outward.potentialEnergy, units::Energy{});
}

TEST_F(GravityFluxWork, CanonicalAndProvisionalFacesChooseFineAndOwnFluxes) {
	// Zero potential in the recipient models a cell outside a masked component;
	// mass crossing its boundary must still receive the component's work.
	setPotential(0, 6);
	for (int sign : {-1, 1}) {
		mesh::Coordinates cell{}, coordinate{};
		cell[0] = sign > 0 ? 1 : 0;
		coordinate[0] = sign > 0 ? 2 : 0;
		setFlux(coordinate, sign * 2, sign * 6);
		GravityWorkFace face;
		face.cell = layout.index(cell);
		face.axis = 0;
		face.sign = sign;
		face.areaFraction = 0.25;
		face.potentialFraction = 2.0 / 3;
		face.neighbor = cells.ranges()[1];
		face.massFlux = faces.ranges()[1];
		auto const canonical = gravity::fluxWork(block, {face}, field.handle(), 0, mass.handle(), 0, dt);
		auto const provisional = gravity::fluxWork(block, {face}, field.handle(), 0, mass.handle(), 0, dt, false);
		EXPECT_NEAR(units::value(canonical.work[face.cell]), -2.4, 1e-14);
		EXPECT_NEAR(units::value(provisional.work[face.cell]), -0.8, 1e-14);
		for (std::size_t i = 0; i < canonical.work.size(); ++i) if (i != face.cell) {
			EXPECT_EQ(canonical.work[i], units::EnergyDensity{});
			EXPECT_EQ(provisional.work[i], units::EnergyDensity{});
		}
	}
}

TEST_F(GravityFluxWork, PhysicalBoundaryUsesMeanAccelerationAndSignedEnergyLedger) {
	setPotential(2, 0, -2, 0);
	setPotential(4, 0, -4, 1);
	for (int sign : {-1, 1}) {
		mesh::Coordinates cell{}, coordinate{};
		cell[0] = sign > 0 ? 1 : 0;
		coordinate[0] = sign > 0 ? 2 : 0;
		setFlux(coordinate, 2, 0);
		GravityWorkFace face;
		face.cell = layout.index(cell);
		face.axis = 0;
		face.sign = sign;
		face.physical = true;
		face.massFlux = block.massFlux.slice(allFaceIndex(layout, 0, coordinate), 1);
		auto const result = gravity::fluxWork(block, {face}, field.handle(), 0, field.handle(), 1, mass.handle(), 0, dt);
		EXPECT_NEAR(units::value(result.work[face.cell]), -0.3, 1e-14);
		auto const transported = result.boundary.outward.potentialEnergy - result.boundary.inward.potentialEnergy;
		EXPECT_NEAR(units::value(transported), sign * 0.8 * units::value(layout.cellMeasure(block.cellWidth)) * (3 + sign * 0.375), 1e-14);
	}
}
