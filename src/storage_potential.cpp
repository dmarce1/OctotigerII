#include "octotigerII/storage/field.hpp"
#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_PARTITIONED_VECTOR(OctoPotential)
HPX_REGISTER_ACTION(octotigerII::storage::ReadAction<octotigerII::units::VelocitySquared>, octoPotentialRead)
HPX_REGISTER_ACTION(octotigerII::storage::WriteAction<octotigerII::units::VelocitySquared>, octoPotentialWrite)
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::ReadAction<octotigerII::units::VelocitySquared>, "octoPotentialRead", std::size_t(-1), std::size_t(-1))
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::WriteAction<octotigerII::units::VelocitySquared>, "octoPotentialWrite", std::size_t(-1), std::size_t(-1))
#endif
