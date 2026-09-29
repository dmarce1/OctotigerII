#include "octotigerII/profilingSnapshots.hpp"
#include "octotigerII/profiling.hpp"

#if defined(OCTOTIGERII_PROFILE_HPX) && defined(HPX_HAVE_APEX)
#include <apex_options.hpp>
#include <hpx/runtime_local/get_locality_id.hpp>
#include <hpx/runtime_distributed/find_all_localities.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/runtime_local/interval_timer.hpp>
#include <hpx/runtime_local/shutdown_function.hpp>
#include <hpx/synchronization/mutex.hpp>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <vector>

namespace octotigerII::profiling {
namespace {
// APEX's screen/CSV reduction can use MPI collectives. Independent locality
// timers must only dump local profiles. Restore normal final-output settings.
struct LocalOutput {
	bool screen = apex::apex_options::use_screen_output();
	int csv = apex::apex_options::use_csv_output();
	int profile = apex::apex_options::use_profile_output();
	LocalOutput() {
		apex::apex_options::use_screen_output(false);
		apex::apex_options::use_csv_output(0);
		apex::apex_options::use_profile_output(1);
	}
	~LocalOutput() {
		apex::apex_options::use_screen_output(screen);
		apex::apex_options::use_csv_output(csv);
		apex::apex_options::use_profile_output(profile);
	}
};

struct Snapshots {
	hpx::mutex mutex;
	bool stopped = false;
	unsigned long long sequence = 0;
	std::filesystem::path root, directory;
	std::string filename;
	std::deque<std::filesystem::path> completed;
	std::unique_ptr<hpx::util::interval_timer> timer;

	bool write() {
		std::lock_guard lock(mutex);
		if (stopped) return false;
		try {
			{
				LocalOutput options;
				apex::dump(false); // cumulative; never reset the measured totals
			}
			auto const source = root / filename;
			std::ifstream profile(source);
			std::string header;
			std::getline(profile, header);
			if (header.find("templated_functions_MULTI_TIME") == std::string::npos)
				throw std::runtime_error("APEX did not write a valid local profile");
			std::ostringstream name;
			name << "snapshot-" << std::setfill('0') << std::setw(8) << ++sequence << ".profile";
			auto const published = directory / name.str();
			auto const pending = directory / (name.str() + ".pending");
			std::filesystem::copy_file(source, pending, std::filesystem::copy_options::overwrite_existing);
			// Same-filesystem rename publishes only a complete copy. A killed dump
			// may damage the root profile, but leaves prior snapshots untouched.
			std::filesystem::rename(pending, published);
			completed.push_back(published);
			if (completed.size() > 2) {
				std::filesystem::remove(completed.front());
				completed.pop_front();
			}
		} catch (std::exception const& error) {
			std::cerr << "APEX periodic snapshots stopped: " << error.what() << '\n';
			stopped = true;
			return false; // Profiling output failure must not abort the simulation.
		}
		return true;
	}

	void stop() {
		timer->stop(true);
		// Wait cooperatively for any active dump and option restoration before
		// APEX finalization. Late timer callbacks see stopped and do no work.
		std::lock_guard lock(mutex);
		stopped = true;
	}
};
Snapshots& localState() { static Snapshots state; return state; }
}

void startSnapshots() {
	auto const* setting = std::getenv("OCTOTIGERII_PROFILE_INTERVAL_SECONDS");
	if (!setting) return;
	char* end = nullptr;
	double const seconds = std::strtod(setting, &end);
	if (end == setting || *end != '\0' || !std::isfinite(seconds) || seconds < 0 || seconds > 86400 || (seconds > 0 && seconds < .01)) {
		std::cerr << "APEX snapshots disabled: OCTOTIGERII_PROFILE_INTERVAL_SECONDS must be 0 or 0.01..86400\n";
		return;
	}
	if (seconds == 0 || apex::apex_options::disable()) return;
	if (apex::apex_options::use_final_output_only() || apex::apex_options::use_tau()) {
		std::cerr << "APEX snapshots disabled by APEX_FINAL_OUTPUT_ONLY or external TAU profiling\n";
		return;
	}
	// These modes can perform collectives or alter instrumentation when their
	// options change. Keep their settings intact; asynchronous dumps here are
	// deliberately limited to ordinary flat profiles.
	if (apex::apex_options::use_taskgraph_output() || apex::apex_options::use_tasktree_output() ||
		apex::apex_options::use_hatchet_output() || apex::apex_options::use_jupyter_support() ||
		apex::apex_options::task_scatterplot()) {
		std::cerr << "APEX snapshots disabled: periodic locality output requires flat profiling mode\n";
		return;
	}
	// Process lifetime; shutdown stops callbacks before APEX or HPX is torn down.
	auto& state = localState();
	if (state.timer) return;
	try {
		state.root = apex::apex_options::output_file_path();
		auto const id = std::to_string(hpx::get_locality_id());
		state.filename = "profile." + id + ".0.0";
		state.directory = state.root / "snapshots" / ("locality-" + id);
		std::filesystem::create_directories(state.directory);
		state.timer = std::make_unique<hpx::util::interval_timer>([] { return localState().write(); },
			static_cast<std::int64_t>(seconds * 1e6), "profiling.snapshot", true);
		hpx::register_pre_shutdown_function([] { localState().stop(); });
		state.timer->start(false);
		std::clog << "APEX cumulative snapshots every " << seconds << " s: " << state.directory << '\n';
	} catch (std::exception const& error) {
		std::cerr << "APEX snapshots disabled: " << error.what() << '\n';
	}
}
void stopLocalSnapshots() {
	auto& state = localState();
	if (state.timer) state.stop();
}
}
HPX_PLAIN_ACTION(octotigerII::profiling::stopLocalSnapshots, octo_stop_profiling_snapshots_action);
namespace octotigerII::profiling {
void stopSnapshots() {
	// Stop periodic work before distributed termination detection; waiting for
	// pre-shutdown alone lets a repeating timer keep the runtime alive forever.
	std::vector<hpx::future<void>> pending;
	for (auto const& locality : hpx::find_all_localities())
		pending.push_back(hpx::async<octo_stop_profiling_snapshots_action>(locality));
	for (auto& task : pending) task.get();
}
}
#else
namespace octotigerII::profiling {
void startSnapshots() {}
void stopSnapshots() {}
}
#endif
