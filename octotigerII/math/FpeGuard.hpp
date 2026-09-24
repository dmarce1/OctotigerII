/** @file
 * @brief Scoped floating-point exception checking around numerical kernels.
 * @ingroup math
 */
#pragma once
#include <cfenv>

#if defined(__GLIBC__) || defined(__linux__)
extern "C" {
int fegetexcept();

int feenableexcept(int);

int fedisableexcept(int);
}
#define hasFpeEnableExcept 1
#else
#define hasFpeEnableExcept 0
#endif


/// Preserves the floating-point environment while checking sensitive calculations.
/// The guard is scoped to the executing thread; it does not validate remote work.
/// @ingroup math
class FpeGuard {
public:

	explicit FpeGuard(int mask = (FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID)) noexcept {
#if hasFpeEnableExcept
		previousEnabled = fegetexcept();
		if (previousEnabled < 0) {
			return;
		}
		enabledHere = (mask & FE_ALL_EXCEPT) & ~previousEnabled;
		if (enabledHere != 0) {
			// Pending masked exceptions belong to the caller. Clear them
			// before unmasking, then restore them after masking on exit.
			std::fegetexceptflag(&previousFlags, enabledHere);
			std::feclearexcept(enabledHere);
			feenableexcept(enabledHere);
		}
#else
		(void) mask;
#endif
	}

	FpeGuard(FpeGuard const&) = delete;

	FpeGuard(FpeGuard&&) = delete;

	~FpeGuard() noexcept {
#if hasFpeEnableExcept
		if (previousEnabled < 0) {
			return;
		}
		int const nowEnabled = fegetexcept();
		int const toDisable = nowEnabled & ~previousEnabled;
		int const toEnable = previousEnabled & ~nowEnabled;
		// Ordinary nested guards find identical masks and perform no writes.
		if (toDisable != 0) {
			fedisableexcept(toDisable);
		}
		if (toEnable != 0) {
			// Do not activate stale exceptions if code inside changed masks.
			std::feclearexcept(toEnable);
			feenableexcept(toEnable);
		}
		if (enabledHere != 0) {
			std::fesetexceptflag(&previousFlags, enabledHere);
		}
#endif
	}

	FpeGuard& operator=(FpeGuard const&) = delete;

	FpeGuard& operator=(FpeGuard&&) = delete;

private:

	int previousEnabled = -1;
	int enabledHere = 0;
	std::fexcept_t previousFlags{};
};
