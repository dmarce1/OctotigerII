# Project instructions

## HPX concurrency

Use HPX threading and synchronization facilities in HPX-enabled code. Standard
blocking synchronization can block HPX worker threads and interfere with its
scheduler.

- Use `hpx::mutex` and `hpx::shared_mutex` for locks.
- Use `hpx::threads::hardware_concurrency()` for hardware thread counts.
- Use HPX tasks, futures, condition variables, and waits when needed; do not
  introduce `std::thread`, `std::async`, or standard blocking synchronization
  primitives in HPX paths.
- Standard RAII lock wrappers such as `std::lock_guard`, `std::unique_lock`, and
  `std::shared_lock` may manage HPX mutexes.
- Standard mutexes and thread utilities are allowed only in the explicitly
  HPX-free serial build, behind the appropriate compile-time guard.
- Ensure every target using these guards receives the HPX compile definition
  and dependency, including when profiling is disabled.
