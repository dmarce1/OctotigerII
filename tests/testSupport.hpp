#pragma once
#include <gtest/gtest.h>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <vector>
#include "octotigerII/units/cgs.hpp"

namespace octotigerII::test {

// Unique across concurrent CTest processes, repeated runs and test shuffling.
class TemporaryDirectory {
public:

	TemporaryDirectory() {
		auto pattern = (std::filesystem::temp_directory_path() / "octotigerII-test-XXXXXX").string();
		std::vector<char> buffer(pattern.begin(), pattern.end());
		buffer.push_back('\0');
		auto* directory = ::mkdtemp(buffer.data());
		if (!directory) throw std::runtime_error("Cannot create test temporary directory");
		path = directory;
	}

	~TemporaryDirectory() {
		std::error_code error;
		std::filesystem::remove_all(path, error);
	}

	TemporaryDirectory(TemporaryDirectory const&) = delete;
	TemporaryDirectory& operator=(TemporaryDirectory const&) = delete;

	std::filesystem::path path;
};


template <typename State>
void expectStateNear(State const& actual, State const& expected, Real relative = 2e-13) {
	expected.forEach([&](auto field, auto q) {
		Real const a = units::value(actual.template get<field>()), b = units::value(q);
		EXPECT_TRUE(std::isfinite(a));
		EXPECT_NEAR(a, b, relative * std::max(Real(1), std::abs(b))) << "field " << int(field);
	});
}

} // namespace octotigerII::test
