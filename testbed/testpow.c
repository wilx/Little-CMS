/* Private pow regression: the accelerated kernel and system fallback.
 * Distributed under the Little CMS MIT license (see ../LICENSE).
 */
#include "lcms2_internal.h"
#include <errno.h>
#include <fenv.h>
#include <float.h>
#include <stdint.h>

static double (*volatile SystemPow)(double, double) = pow;
static double (*volatile CandidatePow)(double, double) = _cmsPow;
static unsigned Failures, Cases;

// Copy a double representation to an integer for exact and ULP-distance comparisons.
static uint64_t Bits(double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return bits;
}

// Construct edge-case inputs, including subnormals and signaling NaNs, from raw bits.
static double FromBits(uint64_t bits)
{
    double x;
    memcpy(&x, &bits, sizeof(x));
    return x;
}

// Compare against system pow: exact result/status for required fallbacks, otherwise at most two ULPs.
static void Check(double x, double y, int exact)
{
    double reference, candidate;
    int flags, error;
    uint64_t a, b, distance;

    feclearexcept(FE_ALL_EXCEPT);
    errno = 0;
    reference = SystemPow(x, y);
    flags = fetestexcept(FE_ALL_EXCEPT);
    error = errno;
    feclearexcept(FE_ALL_EXCEPT);
    errno = 0;
    candidate = CandidatePow(x, y);
    ++Cases;

    /* Only fallback calls promise identical errno and exception flags on
       every supported libm. Fast-path results may differ by a last bit. */
    if (exact && (flags != fetestexcept(FE_ALL_EXCEPT) || error != errno))
        goto fail;
    a = Bits(reference);
    b = Bits(candidate);
    if (a == b) return;
    if (exact) goto fail;
    if (!isfinite(reference) || !isfinite(candidate) ||
        signbit(reference) != signbit(candidate)) goto fail;
    distance = a > b ? a - b : b - a;
    if (distance <= 2) return;
fail:
    if (Failures++ < 10)
        fprintf(stderr, "pow(%a, %a): system=%a candidate=%a exact=%d\n",
                x, y, reference, candidate, exact);
}

// Advance a deterministic xorshift64* generator for reproducible pow input coverage.
static uint64_t RandomBits(uint64_t* state)
{
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * UINT64_C(2685821657736338717);
}

// Run the pow edge cases and deterministic samples across rounding modes; return failure if any comparison fails.
int main(void)
{
    static const int modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    double edges[] = {0., -0., 1., -1., 2., -2., 0.5, -0.5,
        DBL_MIN, -DBL_MIN, DBL_MAX, -DBL_MAX,
        0x1p-65, 0x1p63, 0x1p-54, 512., -512., 1024., -1024.,
        INFINITY, -INFINITY, NAN, 0.41666772630749083, 2.4,
        0., 0., 0.};
    uint64_t state = 2204;
    unsigned mode, i, j;
    int oldMode = fegetround();

    /* This test manipulates IEEE binary64 representations, just as the fast
       kernel does; portable system-pow-only builds may use another format. */
    if (sizeof(double) != sizeof(uint64_t) || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024)
        return 0;

    edges[24] = FromBits(1);
    edges[25] = FromBits(UINT64_C(0x8000000000000001));
    edges[26] = FromBits(UINT64_C(0x7ff0000000000001));

    for (mode = 0; mode < sizeof(modes)/sizeof(modes[0]); ++mode) {
        if (fesetround(modes[mode]) != 0) continue;
        for (i = 0; i < sizeof(edges)/sizeof(edges[0]); ++i) {
            for (j = 0; j < sizeof(edges)/sizeof(edges[0]); ++j) {
                double x = edges[i], y = edges[j];
                /* Determine the initial fallback using bits, so even sNaNs
                   reach the evaluator before any floating-point operation. */
                uint64_t tx = Bits(x) >> 52, ty = (Bits(y) >> 52) & 0x7ff;
                int fallback = tx - 1 >= 0x7ff - 1 ||
                               ty - 0x3be >= 0x43e - 0x3be;
                Check(x, y, fallback);
            }
        }
        /* Positive normal inputs, ordinary and extreme result ranges, plus
           values adjacent to one and the sRGB exponent from the RT fixture. */
        for (i = 0; i < 10000; ++i) {
            double x = ldexp(1. + (RandomBits(&state) >> 11) * 0x1p-53,
                             (int)(RandomBits(&state) % 2046) - 1022);
            double y = ((int)(RandomBits(&state) % 20000001) - 10000000) / 1000000.;
            Check(x, y, 0);
            Check((double)(float)(0.0048 + (RandomBits(&state) >> 11) * 0x1p-53 * 0.2582),
                  0.41666772630749083, 0);
        }
        Check(nextafter(1., 0.), 0.41666772630749083, 0);
        Check(nextafter(1., 2.), 0.41666772630749083, 0);
    }
    fesetround(oldMode);
    printf("Private pow: %u cases, %u failures\n", Cases, Failures);
    return Failures != 0;
}
