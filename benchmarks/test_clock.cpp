// Check CPU accounting behavior rather than asserting a particular timer resolution.
#include "clock.hpp"
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <thread>
#include <type_traits>

static_assert(!std::is_same<bench::ThreadCpuClock::time_point,
                           bench::ProcessCpuClock::time_point>::value,
              "Thread and process timestamps must not be interchangeable");

// Keep the arithmetic workload observable while avoiding a shared hot variable.
static volatile std::uint64_t sink;

// Fail in release builds too; assertions may be disabled in benchmark builds.
static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

// Spend at least the requested amount of thread CPU time doing observable work.
static double burn_cpu(double seconds)
{
    const auto begin = bench::ThreadCpuClock::now();
    const auto deadline = bench::WallClock::now() + std::chrono::seconds(10);
    std::uint64_t value = 1;
    double elapsed;
    do {
        for (unsigned i = 0; i < 100000; ++i)
            value = value * 6364136223846793005ULL + 1;
        elapsed = bench::seconds(bench::ThreadCpuClock::now() - begin);
        require(bench::WallClock::now() < deadline, "CPU clock did not advance within 10 seconds");
    } while (elapsed < seconds);
    sink = value;
    return elapsed;
}

// A sleeping thread should accumulate little CPU; a worker must count toward
// process CPU time without being charged to the thread waiting for that worker.
int main()
{
    try {
        auto wall = bench::WallClock::now();
        auto thread = bench::ThreadCpuClock::now();
        auto process = bench::ProcessCpuClock::now();
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
        const double sleep_thread = bench::seconds(bench::ThreadCpuClock::now() - thread);
        const double sleep_process = bench::seconds(bench::ProcessCpuClock::now() - process);
        const double sleep_wall = bench::seconds(bench::WallClock::now() - wall);
        require(sleep_wall >= 0.20, "Wall clock failed to include sleep");
        require(sleep_thread >= 0 && sleep_thread < sleep_wall * 0.5, "Thread CPU clock counted sleep");
        require(sleep_process >= 0 && sleep_process < sleep_wall * 0.5, "Process CPU clock counted sleep");

        double worker_cpu = 0;
        thread = bench::ThreadCpuClock::now();
        process = bench::ProcessCpuClock::now();
        std::thread worker([&worker_cpu] { worker_cpu = burn_cpu(0.20); });
        worker.join();
        const double waiting_cpu = bench::seconds(bench::ThreadCpuClock::now() - thread);
        const double total_cpu = bench::seconds(bench::ProcessCpuClock::now() - process);
        require(worker_cpu >= 0.20, "Worker thread CPU clock did not advance");
        require(waiting_cpu >= 0 && waiting_cpu < worker_cpu * 0.5, "Worker CPU was charged to waiting thread");
        // Allow clock-accounting granularity at the two measurement boundaries.
        require(total_cpu >= worker_cpu * 0.7, "Process CPU clock missed worker execution");
        require(total_cpu < worker_cpu * 2.0 + 0.1, "Process CPU clock overcounted worker execution");
        std::printf("PASS sleep: wall=%.6f thread_cpu=%.6f process_cpu=%.6f; worker_cpu=%.6f waiting_cpu=%.6f total_cpu=%.6f\n",
                    sleep_wall, sleep_thread, sleep_process, worker_cpu, waiting_cpu, total_cpu);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "FAIL %s\n", error.what());
        return 1;
    }
}
