// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <type_traits>
#include <utility>

#ifdef OCTOTIGERII_PROFILE_HPX
#include <hpx/config.hpp>
#include <hpx/functional.hpp>
#if defined(HPX_HAVE_APEX)
#include <apex_api.hpp>
#include <chrono>
#include <hpx/threading_base/thread_helpers.hpp>
#endif
#endif


namespace octotigerII::profiling {


/// Name the scheduled task, retaining HPX's yield/resume accounting.
/// Names must have static storage duration (normally string literals).
template <typename Function>
auto annotated(Function&& function, char const* name) {
#ifdef OCTOTIGERII_PROFILE_HPX
	return hpx::annotated_function(std::forward<Function>(function), name);
#else
	(void) name;
	return std::decay_t<Function>(std::forward<Function>(function));
#endif
}

/// Record a workload count or other scalar on the calling locality.
inline void sample(char const* name, double value) {
#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
	if (hpx::threads::get_self_ptr() != nullptr) apex::sample_value(name, value);
#else
	(void) name;
	(void) value;
#endif
}


/// A real APEX timer for synchronous computation. Do not span HPX waits,
/// yields, or task migration with this object; use Elapsed for those scopes.
/// Unlike scoped_annotation in HPX 1.11, this starts/stops a nested timer.
class Region {
public:

	explicit Region(char const* name) {
#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
		if (hpx::threads::get_self_ptr() != nullptr) timer_ = apex::start(name);
#else
		(void) name;
#endif
	}

	~Region() {
#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
		if (timer_ != nullptr) apex::stop(timer_);
#endif
	}

	Region(Region const&) = delete;
	Region& operator=(Region const&) = delete;

private:

#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
	apex::profiler* timer_ = nullptr;
#endif
};


/// Wall latency including waits and child work, reported as an APEX sampled
/// counter in nanoseconds. No active APEX timer is carried across suspension.
/// Use a name ending in .wall_ns; these samples are not CPU-time totals.
class Elapsed {
public:

	explicit Elapsed(char const* name) {
#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
		name_ = name;
		begin_ = Clock::now();
#else
		(void) name;
#endif
	}

	~Elapsed() {
#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
		sample(name_, std::chrono::duration<double, std::nano>(Clock::now() - begin_).count());
#endif
	}

	Elapsed(Elapsed const&) = delete;
	Elapsed& operator=(Elapsed const&) = delete;

private:

#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
	using Clock = std::chrono::steady_clock;
	char const* name_;
	Clock::time_point begin_;
#endif
};


}	 // namespace octotigerII::profiling
