// CPU clocks for benchmarks; not part of the installed LCMS API.
#ifndef LCMS_BENCHMARK_CLOCK_HPP
// Keep the clock types usable from multiple benchmark translation units.
#define LCMS_BENCHMARK_CLOCK_HPP

#include <chrono>

namespace bench {

// Cumulative user + kernel CPU time of the calling thread. Subtract timestamps
// taken on the same thread. Sleeping and time spent descheduled are excluded.
struct ThreadCpuClock {
    using duration = std::chrono::nanoseconds;
    using rep = duration::rep;
    using period = duration::period;
    using time_point = std::chrono::time_point<ThreadCpuClock>;
    static constexpr bool is_steady = true;

    // Read the native CPU clock, throwing std::system_error on API failure.
    static time_point now();
};

// Cumulative user + kernel CPU time of all threads in the current process.
// Its elapsed time can exceed wall time when several threads run concurrently.
struct ProcessCpuClock {
    using duration = std::chrono::nanoseconds;
    using rep = duration::rep;
    using period = duration::period;
    using time_point = std::chrono::time_point<ProcessCpuClock>;
    static constexpr bool is_steady = true;

    // Read the native CPU clock, throwing std::system_error on API failure.
    static time_point now();
};

// Let the C++ runtime supply the platform's monotonic elapsed-time clock.
using WallClock = std::chrono::steady_clock;

// Convert an already-subtracted duration, avoiding floating-point timestamps.
template<class Rep, class Period>
double seconds(std::chrono::duration<Rep, Period> elapsed)
{
    return std::chrono::duration<double>(elapsed).count();
}

} // namespace bench
#endif
