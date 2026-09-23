// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/storage/partition.hpp"

#ifdef OCTOTIGERII_WITH_HPX


namespace octotigerII::storage {

void retireField(hpx::id_type id, FieldId field) {
	hpx::get_ptr<StoragePartition>(id).get()->retire(field);
}

}	 // namespace octotigerII::storage


using StorageComponent = hpx::components::component<octotigerII::storage::StoragePartition>;
HPX_REGISTER_COMPONENT(StorageComponent, octotigerII_storage)
HPX_REGISTER_ACTION(octotigerII::storage::RetireAction, octoRetireField)
#define OCTO_STORAGE_TYPE(Type, Name)                                                                                                             \
	HPX_REGISTER_ACTION(octotigerII::storage::CreateAction<Type>, Name##Create)                                                                   \
	HPX_REGISTER_ACTION(octotigerII::storage::ReadAction<Type>, Name##Read)                                                                       \
	HPX_REGISTER_ACTION(octotigerII::storage::WriteAction<Type>, Name##Write)                                                                     \
	HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(octotigerII::storage::ReadAction<Type>, #Name "Read", std::size_t(-1), std::size_t(-1)) \
	HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DEFINITION(octotigerII::storage::WriteAction<Type>, #Name "Write", std::size_t(-1), std::size_t(-1))
#include "octotigerII/storage/types.def"
#undef OCTO_STORAGE_TYPE
#endif
