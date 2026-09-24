#pragma once
#include <gtest/gtest.h>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <vector>
#include "octotigerII/units/cgs.hpp"
#include "octotigerII/config.hpp"

namespace octotigerII::test {
#ifndef OCTOII_TEST_PROBLEM
#define OCTOII_TEST_PROBLEM "sod"
#endif
inline constexpr char problem[] = OCTOII_TEST_PROBLEM;
inline constexpr bool hydro = std::string_view(problem) == "sod" || std::string_view(problem) == "collapse" || std::string_view(problem) == "polytrope" || std::string_view(problem) == "rayleigh-taylor";
inline constexpr bool radiation = std::string_view(problem) == "streaming";
inline constexpr bool gravity = std::string_view(problem) == "gravity-sphere" || std::string_view(problem) == "collapse" || std::string_view(problem) == "polytrope";
inline Config parseConfig(std::vector<std::string> arguments) {
 arguments.insert(arguments.begin(), std::string("--problem.name=") + problem);
 return octotigerII::parseConfig(arguments);
}


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
