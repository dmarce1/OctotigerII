/** @file
 * @brief Application field schema over the mesh-independent typed store.
 * @ingroup storage
 */
#pragma once

#include "octotigerII/storage/columns.hpp"
#include "octotigerII/subgrid/subgrid.hpp"

namespace octotigerII {

// Application field registration is separate from the generic typed store.
/// Serializable handles for the enabled hydro, radiation, and gravity schemas.
/// Radiation stores E and ndim physical EnergyFlux quantities.
/// @ingroup storage
class FieldDirectory {
public:
	storage::ColumnHandle<hydro::ConservedState> hydro;
	storage::ColumnHandle<radiation::RadiationSystem::State> radiation;
	storage::ColumnHandle<gravity::State> gravity;
	storage::FieldHandle<units::Density> density;
	storage::ColumnHandle<hydro::ConservedFlux> hydroFlux;
	storage::ColumnHandle<radiation::RadiationSystem::Flux> radiationFlux;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & hydro & radiation & gravity & density & hydroFlux & radiationFlux;
	}
};

/// Registers the enabled application fields in a shared generic store.
/// Physics field names and schema belong here, outside StoragePartition.
/// @ingroup storage
class FieldRepository {
public:
	FieldRepository(Config const& config, storage::Layout const& layout, std::vector<storage::Locality> const& localities)
	  : store_(localities) {
		if (config.amr.enabled) {
			storage::Layout fluxLayout(std::vector<std::size_t>(layout.ranges().size(), boundaryFluxCount(config.mesh.cells)), localities.size());
			if constexpr (build::hydro) {
				hydroFlux_ = std::make_unique<storage::ColumnFields<hydro::ConservedFlux>>(fluxLayout, store_, "hydro.boundaryFlux");
				directory_.hydroFlux = hydroFlux_->handle();
			}
			if constexpr (build::radiation) {
				radiationFlux_ = std::make_unique<storage::ColumnFields<radiation::RadiationSystem::Flux>>(fluxLayout, store_, "radiation.boundaryFlux");
				directory_.radiationFlux = radiationFlux_->handle();
			}
		}
		if constexpr (build::hydro) {
			hydro_ = std::make_unique<storage::ColumnFields<hydro::ConservedState>>(layout, store_, "hydro");
			directory_.hydro = hydro_->handle();
		}
		if constexpr (build::radiation) {
			radiation_ = std::make_unique<storage::ColumnFields<radiation::RadiationSystem::State>>(layout, store_, "radiation");
			directory_.radiation = radiation_->handle();
		}
		if constexpr (build::gravity) {
			gravity_ = std::make_unique<storage::ColumnFields<gravity::State>>(layout, store_, "gravity");
			directory_.gravity = gravity_->handle();
			if constexpr (!build::hydro) {
				density_ = std::make_unique<storage::Field<units::Density>>(layout, store_, 2, "density");
				directory_.density = density_->handle();
			}
		}
	}

	/// Return the enabled application field handles.
	FieldDirectory const& directory() const {
		return directory_;
	}

private:
	storage::PartitionSet store_;
	FieldDirectory directory_;
	std::unique_ptr<storage::ColumnFields<hydro::ConservedState>> hydro_;
	std::unique_ptr<storage::ColumnFields<radiation::RadiationSystem::State>> radiation_;
	std::unique_ptr<storage::ColumnFields<gravity::State>> gravity_;
	std::unique_ptr<storage::Field<units::Density>> density_;
	std::unique_ptr<storage::ColumnFields<hydro::ConservedFlux>> hydroFlux_;
	std::unique_ptr<storage::ColumnFields<radiation::RadiationSystem::Flux>> radiationFlux_;
};

}	 // namespace octotigerII
