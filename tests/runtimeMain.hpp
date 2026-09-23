#pragma once
#include <string>
#include <vector>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/hpx_init.hpp>
#endif

inline int runtimeMain(int argc, char** argv, int (*run)(int, char**)) {
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<char*> application{argv[0]}, runtime{argv[0]};
	for (int i = 1; i < argc; ++i)
		(std::string(argv[i]).starts_with("--hpx:") ? runtime : application).push_back(argv[i]);
	int const count = static_cast<int>(runtime.size());
	int const applicationCount = static_cast<int>(application.size());
	runtime.push_back(nullptr);
	application.push_back(nullptr);
	return hpx::init(
		[&](int, char**) {
			int const result = run(applicationCount, application.data());
			hpx::finalize();
			return result;
		},
		count, runtime.data());
#else
	return run(argc, argv);
#endif
}
