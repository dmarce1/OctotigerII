#include "octotigerII/storage/field.hpp"
#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_PARTITIONED_VECTOR(OctoMomentum)
HPX_REGISTER_ACTION(octotigerII::storage::ReadAction<octotigerII::units::MomentumDensity>, octoMomentumRead)
HPX_REGISTER_ACTION(octotigerII::storage::WriteAction<octotigerII::units::MomentumDensity>, octoMomentumWrite)
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::ReadAction<octotigerII::units::MomentumDensity>, "octoMomentumRead", std::size_t(-1), std::size_t(-1))
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::WriteAction<octotigerII::units::MomentumDensity>, "octoMomentumWrite", std::size_t(-1), std::size_t(-1))
#endif
