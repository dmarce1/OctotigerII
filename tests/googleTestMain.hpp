#pragma once
#include <gtest/gtest.h>
#include "runtimeMain.hpp"

// CMake discovery only lists tests; it must not start HPX or contact AGAS.
// GoogleTest consumes its flags before HPX sees the remaining runtime flags.
inline int googleTestMain(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	if (::testing::GTEST_FLAG(list_tests)) return RUN_ALL_TESTS();
	return runtimeMain(argc, argv, [](int, char**) { return RUN_ALL_TESTS(); });
}
