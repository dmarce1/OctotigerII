#include "octotigerII/subgrid/subgrid.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "octotigerII/problems.hpp"


namespace octotigerII {
Snapshot initialSnapshot(Config const& c, mesh::BlockLocation location) {
	using std::ldexp;

	Snapshot data_;
	c.validate();
	if (location.level != c.mesh.level) throw std::invalid_argument("Subgrid location/config mismatch");
	data_.location = location;
	data_.layout = mesh::MeshLayout(c.mesh.cells);
	auto const blockWidth = (c.mesh.upper - c.mesh.lower) * ldexp(Real(1), -c.mesh.level);
	data_.cellWidth = blockWidth / Real(c.mesh.cells);
	for (int axis = 0; axis < ndim; ++axis) {
		if (location.coordinates[axis] < 0 || location.coordinates[axis] >= (1 << c.mesh.level)) throw std::invalid_argument("Subgrid coordinate out of bounds");
		data_.lower[axis] = c.mesh.lower + Real(location.coordinates[axis]) * blockWidth;
	}
	data_.hydroEnabled = c.hydroEnabled();
	data_.radiationEnabled = c.radiationEnabled();
	data_.gravityEnabled = c.gravityEnabled();
	if constexpr (build::hydro) data_.hydro = hydro::Fields(data_.layout, data_.cellWidth, data_.lower);
	if constexpr (build::radiation) data_.radiation = radiation::Fields(data_.layout, data_.cellWidth, data_.lower);
	if constexpr (build::gravity) {
		data_.gravity = gravity::Fields(data_.layout, data_.cellWidth, data_.lower);
		if (!data_.hydroEnabled) data_.density = mesh::PatchData<units::Density>(data_.layout, data_.cellWidth, data_.lower);
	}
	initializeProblem(data_, c);
	return data_;
}

}	 // namespace octotigerII
