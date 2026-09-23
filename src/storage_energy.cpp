#include "octotigerII/storage/field.hpp"
#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_PARTITIONED_VECTOR(OctoEnergy)
HPX_REGISTER_ACTION(octotigerII::storage::ReadAction<octotigerII::units::EnergyDensity>, octoEnergyRead)
HPX_REGISTER_ACTION(octotigerII::storage::WriteAction<octotigerII::units::EnergyDensity>, octoEnergyWrite)
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::ReadAction<octotigerII::units::EnergyDensity>, "octoEnergyRead", std::size_t(-1), std::size_t(-1))
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::WriteAction<octotigerII::units::EnergyDensity>, "octoEnergyWrite", std::size_t(-1), std::size_t(-1))
#endif
