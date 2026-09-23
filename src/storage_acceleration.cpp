#include "octotigerII/storage/field.hpp"
#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_PARTITIONED_VECTOR(OctoAcceleration)
HPX_REGISTER_ACTION(octotigerII::storage::ReadAction<octotigerII::units::Acceleration>, octoAccelerationRead)
HPX_REGISTER_ACTION(octotigerII::storage::WriteAction<octotigerII::units::Acceleration>, octoAccelerationWrite)
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::ReadAction<octotigerII::units::Acceleration>, "octoAccelerationRead", std::size_t(-1), std::size_t(-1))
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::WriteAction<octotigerII::units::Acceleration>, "octoAccelerationWrite", std::size_t(-1), std::size_t(-1))
#endif
