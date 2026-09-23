#include "octotigerII/storage/field.hpp"
#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_PARTITIONED_VECTOR(OctoDensity)
HPX_REGISTER_ACTION(octotigerII::storage::ReadAction<octotigerII::units::Density>, octoDensityRead)
HPX_REGISTER_ACTION(octotigerII::storage::WriteAction<octotigerII::units::Density>, octoDensityWrite)
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::ReadAction<octotigerII::units::Density>, "octoDensityRead", std::size_t(-1), std::size_t(-1))
HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(
    octotigerII::storage::WriteAction<octotigerII::units::Density>, "octoDensityWrite", std::size_t(-1), std::size_t(-1))
#endif
