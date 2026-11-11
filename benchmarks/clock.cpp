#include "clock.hpp"

#include <system_error>

#ifdef _WIN32
#ifndef NOMINMAX
// Windows min/max macros interfere with standard C++ library names.
#define NOMINMAX
#endif
#include <windows.h>
#include <cstdint>
#include <limits>
#include <stdexcept>
#else
#include <cerrno>
#include <time.h>
#endif

namespace {

#ifdef _WIN32
// Assemble FILETIME's two words without alignment or aliasing assumptions.
std::uint64_t ticks(const FILETIME& time)
{
    return (std::uint64_t(time.dwHighDateTime) << 32) | time.dwLowDateTime;
}

// Read execution time in 100-ns units; creation/exit timestamps are unused.
std::chrono::nanoseconds cpu_time(bool thread)
{
    FILETIME creation, exit, kernel, user;
    const BOOL ok = thread
        ? GetThreadTimes(GetCurrentThread(), &creation, &exit, &kernel, &user)
        : GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user);
    if (!ok)
        throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
                                thread ? "GetThreadTimes" : "GetProcessTimes");

    const auto kernel_ticks = ticks(kernel);
    const auto user_ticks = ticks(user);
    const auto limit = std::uint64_t(std::numeric_limits<std::chrono::nanoseconds::rep>::max()) / 100;
    if (kernel_ticks > limit || user_ticks > limit - kernel_ticks)
        throw std::overflow_error("CPU time exceeds nanosecond clock range");
    return std::chrono::nanoseconds(static_cast<std::chrono::nanoseconds::rep>(kernel_ticks + user_ticks) * 100);
}
#else
// Preserve integer seconds/nanoseconds until timestamps have been subtracted.
std::chrono::nanoseconds cpu_time(clockid_t clock)
{
    timespec time;
    if (clock_gettime(clock, &time) != 0)
        throw std::system_error(errno, std::generic_category(), "clock_gettime");
    return std::chrono::seconds(time.tv_sec) + std::chrono::nanoseconds(time.tv_nsec);
}
#endif

} // namespace

// The backend is private so benchmark call sites do not use native clock IDs.
bench::ThreadCpuClock::time_point bench::ThreadCpuClock::now()
{
#ifdef _WIN32
    return time_point(cpu_time(true));
#else
    return time_point(cpu_time(CLOCK_THREAD_CPUTIME_ID));
#endif
}

// Process accounting includes every process thread, including runtime helpers.
bench::ProcessCpuClock::time_point bench::ProcessCpuClock::now()
{
#ifdef _WIN32
    return time_point(cpu_time(false));
#else
    return time_point(cpu_time(CLOCK_PROCESS_CPUTIME_ID));
#endif
}
