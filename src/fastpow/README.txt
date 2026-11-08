Arm scalar pow kernel used privately by Little CMS
================================================

Upstream: https://github.com/ARM-software/optimized-routines
Commit: 20859d3fd9e797dde592c489c8ef7d383f778620
Inputs: math/pow.c, math/pow_common.h, math/math_config.h,
        math/exp_data.c, math/pow_log_data.c.
License: MIT option of the upstream dual license; see LICENSE.
Original copyright and SPDX notices are retained in each derived file.

../cmspow.c retains the ordinary positive-normal logarithm/product/exponential
path, with the upstream coefficients and arithmetic. All other input and
result ranges call system pow. There is no new approximation or LUT transfer
curve. config.h is a minimal subset of the upstream configuration;
exp_data.h and log_data.h are the upstream tables with private static names.

The fast path requires positive, normal, finite x, 2^-65 <= |y| < 2^63,
and 2^-54 <= |y * log(x)| < 512 (tested on the computed high part).
Negative bases, zeros, subnormals, infinities, NaNs, very small or large
exponents, and results near overflow/underflow use system pow.

On glibc x86-64 with GCC >= 5 or Clang >= 7, a GNU IFUNC resolver chooses the
AVX2/FMA implementation once when the symbol is resolved. CPUs without both
features use system pow. Contraction follows the compiler and build settings;
the upstream explicit fma() calls remain fused independently of contraction.
Other platforms, fast-math builds, and builds defining CMS_DONT_USE_FAST_POW
use system pow. No public API, exported symbol, or structure layout is added.
The inline front end preserves compiler folding of constant powers and
constant-exponent squares, avoiding unnecessary kernel calls for these.

This is a scalar double-precision implementation, not a vector approximation.
Results can differ in the last bits from system pow; universal bit identity
with every compiler setting or platform/libm version is not promised.
