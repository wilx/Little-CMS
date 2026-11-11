Standalone LCMS benchmark clocks
================================

clock.hpp exposes three clocks with the std::chrono interface:

  bench::ThreadCpuClock  Calling thread's cumulative user + kernel execution time.
  bench::ProcessCpuClock All process threads' cumulative user + kernel time.
  bench::WallClock       Monotonic elapsed time, including waits/descheduling.

Usage:

  const auto begin = bench::ThreadCpuClock::now();
  work();
  const double elapsed = bench::seconds(bench::ThreadCpuClock::now() - begin);

Take both ThreadCpuClock readings on the same thread. Clock types have distinct
time_point types to prevent mixing process, thread and wall timestamps. CPU
timestamps use integer nanoseconds; subtraction precedes conversion to seconds.
The clocks measure time, not hardware instruction/cycle counters.

clock.cpp uses clock_gettime(CLOCK_THREAD_CPUTIME_ID / CLOCK_PROCESS_CPUTIME_ID)
on POSIX, and GetThreadTimes / GetProcessTimes on Windows. Windows CPU accounting
adds the user/kernel FILETIME values and converts their 100-nanosecond units.
WallClock delegates to the standard library's steady_clock. API failures throw
std::system_error; they never silently substitute wall time or zero CPU time.

Representational units do not guarantee measurement resolution. Windows samples
should batch enough work to exceed the effective CPU-accounting granularity.
CPU clocks exclude descheduling but do not remove cache, frequency or SMT effects.

Build separately from the LCMS library:

  cmake -S benchmarks -B /tmp/lcms-benchmarks -DCMAKE_BUILD_TYPE=Release
  cmake --build /tmp/lcms-benchmarks --parallel 4
  ctest --test-dir /tmp/lcms-benchmarks --output-on-failure

The same CMake project supports Windows with MSVC or MinGW and OpenMP.
On multi-configuration generators add --config Release to the build and
-C Release to ctest. test_clock checks sleep exclusion and verifies that worker
execution is counted by the process clock but not by the waiting thread's clock.

lab2rgb is the maintained version of the captured-fixture comparison driver.
Historical study sources, binaries, hashes and measured results remain untouched.
The fixture comes from rt-lab2rgbout-packet/fixtures (235840 pixels, 32 rows).
It uses LCMS 2.20's public API and loads each requested library independently.
The same NOOPTIMIZE/NOCACHE/BPC transform settings and bit-exact reference checks
are retained. This is an LCMS transform benchmark, not a full RawTherapee run.

Create libraries.tsv with one label, a TAB, and an absolute library path per line.
Paths can contain spaces; Windows builds accept DLLs through LoadLibraryA.
Use distinct DLL filenames for separate Windows variants and ensure dependencies
are available. Non-ASCII Windows paths must be representable in the active ANSI
code page. POSIX builds use dlopen. No LCMS library is linked into the executable.

  /tmp/lcms-benchmarks/lab2rgb libraries.tsv \
      rt-lab2rgbout-packet/fixtures 100 12345 1 > timings.csv

Arguments are library list, fixture directory, randomized rounds, random seed,
and optional fixture passes per sample (default 1). Each round measures every
variant at one and four threads; three warmup rounds are discarded. Use e.g.
100 passes for longer samples on a coarse CPU clock. The CSV records passes;
seconds columns cover the whole sample, and ns/pixel divides by pixels * passes.

Each worker owns the same rows across passes, and its CPU interval excludes the
final parallel-region barrier. The process interval includes parallel-region
entry/exit and any other process threads. Setup and correctness checks are not
timed. For controlled comparisons, configure affinity and the OpenMP environment
outside the executable, e.g. OMP_WAIT_POLICY=PASSIVE and OMP_PROC_BIND=TRUE.

Windows CPU time documentation:
https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-getthreadtimes
https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-getprocesstimes
