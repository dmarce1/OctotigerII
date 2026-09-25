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
	std::vector<storage::FieldHandle<units::Density>> species;
	std::vector<storage::FieldHandle<units::MassFlux>> speciesFlux;
	storage::FieldHandle<units::MassFlux> massFlux;
	storage::FieldHandle<units::VelocitySquared> oldPotential;
	storage::ColumnHandle<gravity::State> oldGravity;
	storage::FieldHandle<units::EnergyDensity> gravityKickWork;
	storage::ColumnHandle<radiation::RadiationSystem::State> radiation;
	storage::ColumnHandle<gravity::State> gravity;
	storage::FieldHandle<units::Density> density;
	storage::ColumnHandle<hydro::ConservedFlux> hydroFlux;
	storage::ColumnHandle<radiation::RadiationSystem::Flux> radiationFlux;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & species & speciesFlux & massFlux & oldPotential & oldGravity & gravityKickWork;
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
		unsigned const banks = config.amr.enabled && config.timestep.refinement ? (config.gravityEnabled() ? 4 : 3) : 2;
		if (config.massFractions.enabled) {
			for (auto const& s : config.massFractions.species) {
				species_.push_back(std::make_unique<storage::Field<units::Density>>(layout, store_, banks, "massFractions." + s.name));
				directory_.species.push_back(species_.back()->handle());
			}
		}
		if (config.hydroEnabled() && (config.massFractions.enabled || config.gravityEnabled())) {
			storage::Layout faceLayout(std::vector<std::size_t>(layout.ranges().size(), allFaceCount(config.mesh.cells)), localities.size());
			massFlux_ = std::make_unique<storage::Field<units::MassFlux>>(faceLayout, store_, 1, "hydro.massFlux");
			directory_.massFlux = massFlux_->handle();
		}
		if (config.hydroEnabled() && config.gravityEnabled()) {
			oldGravity_ = std::make_unique<storage::ColumnFields<gravity::State>>(layout, store_, "gravity.old");
			gravityKickWork_ = std::make_unique<storage::Field<units::EnergyDensity>>(layout, store_, 1, "gravity.kickWork");
			directory_.oldGravity = oldGravity_->handle();
			directory_.oldPotential = std::get<0>(directory_.oldGravity.fields);
			directory_.gravityKickWork = gravityKickWork_->handle();
		}
		if (config.amr.enabled) {
			storage::Layout fluxLayout(std::vector<std::size_t>(layout.ranges().size(), boundaryFluxCount(config.mesh.cells)), localities.size());
			for (auto const& s : config.massFractions.species) if (config.massFractions.enabled) {
				speciesFlux_.push_back(std::make_unique<storage::Field<units::MassFlux>>(fluxLayout, store_, 2, "massFractions.boundaryFlux." + s.name));
				directory_.speciesFlux.push_back(speciesFlux_.back()->handle());
			}
			if (build::hydro && config.hydroEnabled()) {
				hydroFlux_ = std::make_unique<storage::ColumnFields<hydro::ConservedFlux>>(fluxLayout, store_, "hydro.boundaryFlux");
				directory_.hydroFlux = hydroFlux_->handle();
			}
			if (build::radiation && config.radiationEnabled()) {
				radiationFlux_ = std::make_unique<storage::ColumnFields<radiation::RadiationSystem::Flux>>(fluxLayout, store_, "radiation.boundaryFlux");
				directory_.radiationFlux = radiationFlux_->handle();
			}
		}
		if (build::hydro && config.hydroEnabled()) {
			hydro_ = std::make_unique<storage::ColumnFields<hydro::ConservedState>>(layout, store_, "hydro", config.massFractions.enabled, banks);
			directory_.hydro = hydro_->handle();
			if (config.massFractions.enabled)
				for (std::size_t s = 0; s < config.massFractions.species.size(); ++s)
					if (!config.massFractions.species[s].tracer()) std::get<0>(directory_.hydro.fields).sumSources.push_back(directory_.species[s]);
		}
		if (build::radiation && config.radiationEnabled()) {
			radiation_ = std::make_unique<storage::ColumnFields<radiation::RadiationSystem::State>>(layout, store_, "radiation", false, banks);
			directory_.radiation = radiation_->handle();
		}
		if (build::gravity && config.gravityEnabled()) {
			gravity_ = std::make_unique<storage::ColumnFields<gravity::State>>(layout, store_, "gravity", false, banks);
			directory_.gravity = gravity_->handle();
			if (!config.hydroEnabled()) {
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
	std::vector<std::unique_ptr<storage::Field<units::Density>>> species_;
	std::vector<std::unique_ptr<storage::Field<units::MassFlux>>> speciesFlux_;
	std::unique_ptr<storage::Field<units::MassFlux>> massFlux_;
	std::unique_ptr<storage::ColumnFields<gravity::State>> oldGravity_;
	std::unique_ptr<storage::Field<units::EnergyDensity>> gravityKickWork_;
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
