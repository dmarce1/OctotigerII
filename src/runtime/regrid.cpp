/** @file
 * @brief Mesh replacement and conservative regrid energy recovery.
 */
#include "internal.hpp"
#include "octotigerII/amr/regridSelection.hpp"

namespace octotigerII {

void Runtime::Impl::install(std::vector<mesh::BlockLocation> const& leaves, amr::Hierarchy const& source) {
	bool const conserveEnergy = config.hydroEnabled() && config.gravityEnabled() && config.gravity.conserveRegridEnergy;
	auto nextTopology = std::make_unique<CartesianTopology>(config, localities.size(), leaves);
	auto nextFields = std::make_unique<FieldRepository>(config, nextTopology->storageLayout(), localities);
	auto const& directory = nextFields->directory();
	for (auto const& block : nextTopology->blocks()) {
		auto const snapshot = source.transfer(block.location);
		initializeSpecies(directory, block, snapshot);
		if (build::hydro && config.hydroEnabled()) initializeFields(directory.hydro, block.interior, snapshot.hydro);
		if (conserveEnergy) {
			// Reuse the inactive kick register for conservatively transferred E+rho*phi/2.
			// Hydro remains gas-only for mesh selection, EOS calls and shadow estimates.
			auto const energy = source.transferGasGravityEnergy(block.location);
			auto output = directory.gravityKickWork.output(block.interior, 0);
			std::copy(energy.begin(), energy.end(), output.data());
			directory.gravityKickWork.commit(block.interior, 0, output);
		}
		if (build::radiation && config.radiationEnabled()) initializeFields(directory.radiation, block.interior, snapshot.radiation);
		if (build::gravity && config.gravityEnabled()) {
			initializeFields(directory.gravity, block.interior, snapshot.gravity);
			if (!config.hydroEnabled()) {
				auto output = directory.density.output(block.interior, 0);
				std::copy(snapshot.density.values().begin(), snapshot.density.values().end(), output.data());
				directory.density.commit(block.interior, 0, output);
			}
		}
	}
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<hpx::id_type>> pending;
	for (std::size_t i = 0; i < localities.size(); ++i)
		pending.push_back(hpx::new_<LocalExecutor>(localities[i], config, nextTopology->blocks(), directory, i));
	auto nextExecutors = collect(pending);
	executors.swap(nextExecutors);
#else
	auto nextExecutor = std::make_unique<LocalExecutor>(config, nextTopology->blocks(), directory, 0);
	executor.swap(nextExecutor);
#endif
#if OCTOTIGERII_GRAVITY
	gravitySolver.reset();
#endif
	topology.swap(nextTopology);
	fields.swap(nextFields);
	bank = 0;
	++generation;
	gravityReady = false;
	regridEnergyPending = conserveEnergy;
}

void Runtime::Impl::recoverRegridEnergy(unsigned targetBank) {
	if (!regridEnergyPending) return;
#if OCTOTIGERII_HYDRO && OCTOTIGERII_GRAVITY
	auto const& directory = fields->directory();
	hydro::HydroSystem const gas(config);
	for (auto const& block : topology->blocks()) {
		auto const combined = directory.gravityKickWork.read(block.interior, 0).get();
		auto const gravity = directory.gravity.read(block.interior, targetBank).get();
		auto const input = directory.hydro.read(block.interior, targetBank).get();
		auto output = directory.hydro.output(block.interior, targetBank);
		for (std::size_t i = 0; i < block.interior.count; ++i) {
			auto state = input.at(i);
			state.totalEnergy() = combined.data()[i] - 0.5 * state.density() * gravity.at(i).potential();
			gas.synchronize(state);
			if (!gas.admissible(state)) throw std::runtime_error("Regrid gravity-energy recovery produced an inadmissible gas state");
			output.put(i, state);
		}
		directory.hydro.commit(block.interior, targetBank, output);
	}
#else
	(void) targetBank;
	throw std::logic_error("Gravity energy recovery requires hydro and gravity in this build");
#endif
}

bool Runtime::regrid(units::Time nextStep, bool force) {
	std::lock_guard guard(impl_->apiMutex);
	auto const& config = impl_->config;
	if (!config.amr.enabled) return false;
	if (impl_->gravityEnergyActive) throw std::logic_error("Cannot regrid inside a gravity energy step");
	if (impl_->regridEnergyPending) throw std::logic_error("Solve gravity before regridding again");
	if (!(nextStep >= units::Time{}) || !units::finite(nextStep)) throw std::invalid_argument("Invalid AMR lookahead timestep");
	bool due = force || !impl_->regridInitialized || impl_->time.step - impl_->lastRegridStep >= std::uint64_t(config.amr.regridEvery);
	for (int d = 0; d < ndim; ++d)
		if (impl_->travel[d] + impl_->signalSpeed[d] * nextStep > impl_->travelBudget[d] * (1 + 64 * epsilonR)) due = true;
	if (!due) return false;
	profiling::Elapsed profile("amr.regrid.wall_ns");
	bool changed = false;
	auto const horizon = Real(config.amr.regridEvery) * nextStep;
	auto const original = impl_->exportSnapshots();
	amr::Hierarchy const source(config, original);
	auto candidate = original;
	std::vector<mesh::BlockLocation> leaves;
	std::unique_ptr<amr::Hierarchy> nextShadow;
	std::array<units::Velocity, ndim> bufferSpeed{};
	// Every transfer reads the same immutable pre-regrid state. Build the full
	// Morton-ordered destination before replacing any published field handles.
	for (int pass = 0; pass <= config.amr.maxLevel + 1; ++pass) {
		auto shadow = std::make_unique<amr::Hierarchy>(config, candidate);
		auto const& estimate = pass == 0 && impl_->shadow ? *impl_->shadow : *shadow;
		auto const selected = amr::selectMesh(config, candidate, estimate, horizon, impl_->criteria, pass == 0 && impl_->regridInitialized, shadow.get());
		if (pass == 0) bufferSpeed = selected.signalSpeed;
		if (!selected.changed) {
			nextShadow = std::move(shadow);
			break;
		}
		leaves = selected.leaves;
		candidate.clear();
		for (auto const& leaf : leaves)
			candidate.push_back(source.transfer(leaf));
		changed = true;
		if (pass == config.amr.maxLevel + 1) throw std::runtime_error("AMR refinement failed to settle");
	}
	if (changed) {
		if (config.hydroEnabled() && config.gravityEnabled() && config.gravity.conserveRegridEnergy &&
			(!impl_->gravityReady || impl_->gravityTime != impl_->time.time))
			throw std::logic_error("Energy-conserving regrid requires synchronized gravity on the old mesh");
		impl_->install(leaves, source);
	}
	impl_->shadow = std::move(nextShadow);
	impl_->lastRegridStep = impl_->time.step;
	impl_->regridInitialized = true;
	impl_->travel.fill({});
	for (int d = 0; d < ndim; ++d)
		impl_->travelBudget[d] = config.amr.signalBuffer * bufferSpeed[d] * horizon;
	return changed;
}

} // namespace octotigerII
