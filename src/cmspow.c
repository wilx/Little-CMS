/*
 * Double-precision x^y function.
 *
 * Copyright (c) 2018-2026, Arm Limited.
 * SPDX-License-Identifier: MIT OR Apache-2.0 WITH LLVM-exception
 */

/* Adapted for Little CMS from Arm optimized-routines commit
   20859d3fd9e797dde592c489c8ef7d383f778620. MIT license option selected.
   See fastpow/LICENSE and fastpow/README.txt. */

#include "lcms2_internal.h"
#include <float.h>

/* GNU IFUNC is supported by the glibc ELF loader. Other platforms retain
   system pow. CPU detection runs once, before any calls to the kernel.
   Define CMS_DONT_USE_FAST_POW to disable this implementation at build time. */
#if !defined(CMS_DONT_USE_FAST_POW) && defined(__GLIBC__) && \
    defined(__x86_64__) && !defined(__FAST_MATH__) && \
    DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 && \
    ((defined(__GNUC__) && __GNUC__ >= 5) || \
     (defined(__clang__) && __clang_major__ >= 7))

#include "fastpow/config.h"
#include "fastpow/exp_data.h"
#include "fastpow/log_data.h"

/* Keep contraction local to this kernel, even when LCMS is compiled with
   -ffp-contract=off. The remaining library does not require AVX2 or FMA. */
#if defined(__clang__)
#define CMS_POW_TARGET __attribute__((target("avx2,fma")))
#define CMS_POW_CONTRACT _Pragma("clang fp contract(fast)")
#else
#define CMS_POW_TARGET __attribute__((target("avx2,fma"), optimize("fp-contract=fast")))
#define CMS_POW_CONTRACT
#endif

/*
 * Double-precision pow helper functions.
 *
 * Copyright (c) 2026, Arm Limited.
 * SPDX-License-Identifier: MIT OR Apache-2.0 WITH LLVM-exception
 */

static inline uint32_t
top12 (double x)
{
  return asuint64 (x) >> 52;
}

#define T cms_pow_log_data.tab
#define A cms_pow_log_data.poly
#define Ln2hi cms_pow_log_data.ln2hi
#define Ln2lo cms_pow_log_data.ln2lo
#define N (1 << POW_LOG_TABLE_BITS)
#define OFF 0x3fe6955500000000

/* Compute y+TAIL = log(x) where the rounded result is y and TAIL has about
   additional 15 bits precision.  IX is the bit representation of x, but
   normalized in the subnormal range using the sign bit for the exponent.  */
CMS_POW_TARGET static inline double
log_inline (uint64_t ix, double *tail)
{
  CMS_POW_CONTRACT
  double z, r, y, invc, logc, logctail, kd, hi, t1, t2, lo, lo1, lo2, p;
  uint64_t iz, tmp;
  int k, i;

  /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.
     The range is split into N subintervals.
     The ith subinterval contains z and c is near its center.  */
  tmp = ix - OFF;
  i = (tmp >> (52 - POW_LOG_TABLE_BITS)) % N;
  k = (int64_t) tmp >> 52; /* arithmetic shift.  */
  iz = ix - (tmp & 0xfffULL << 52);
  z = asdouble (iz);
  kd = (double) k;

  /* log(x) = k*Ln2 + log(c) + log1p(z/c-1).  */
  invc = T[i].invc;
  logc = T[i].logc;
  logctail = T[i].logctail;

  /* Note: 1/c is j/N or j/N/2 where j is an integer in [N,2N) and
     |z/c - 1| < 1/N, so r = z/c - 1 is exactly representible.  */
#if HAVE_FAST_FMA
  r = fma (z, invc, -1.0);
#else
  /* Split z such that rhi, rlo and rhi*rhi are exact and |rlo| <= |r|.  */
  double zhi = asdouble ((iz + (1ULL << 31)) & (-1ULL << 32));
  double zlo = z - zhi;
  double rhi = zhi * invc - 1.0;
  double rlo = zlo * invc;
  r = rhi + rlo;
#endif

  /* k*Ln2 + log(c) + r.  */
  t1 = kd * Ln2hi + logc;
  t2 = t1 + r;
  lo1 = kd * Ln2lo + logctail;
  lo2 = t1 - t2 + r;

  /* Evaluation is optimized assuming superscalar pipelined execution.  */
  double ar, ar2, ar3, lo3, lo4;
  ar = A[0] * r; /* A[0] = -0.5.  */
  ar2 = r * ar;
  ar3 = r * ar2;
  /* k*Ln2 + log(c) + r + A[0]*r*r.  */
#if HAVE_FAST_FMA
  hi = t2 + ar2;
  lo3 = fma (ar, r, -ar2);
  lo4 = t2 - hi + ar2;
#else
  double arhi = A[0] * rhi;
  double arhi2 = rhi * arhi;
  hi = t2 + arhi2;
  lo3 = rlo * (ar + arhi);
  lo4 = t2 - hi + arhi2;
#endif
  /* p = log1p(r) - r - A[0]*r*r.  */
#if POW_LOG_POLY_ORDER == 8
  p = (ar3
       * (A[1] + r * A[2]
	  + ar2 * (A[3] + r * A[4] + ar2 * (A[5] + r * A[6]))));
#endif
  lo = lo1 + lo2 + lo3 + lo4 + p;
  y = hi + lo;
  *tail = hi - y + lo;
  return y;
}

#undef N
#undef T
#define N (1 << EXP_TABLE_BITS)
#define InvLn2N cms_pow_exp_data.invln2N
#define NegLn2hiN cms_pow_exp_data.negln2hiN
#define NegLn2loN cms_pow_exp_data.negln2loN
#define Shift cms_pow_exp_data.shift
#define T cms_pow_exp_data.tab
#define C2 cms_pow_exp_data.poly[5 - EXP_POLY_ORDER]
#define C3 cms_pow_exp_data.poly[6 - EXP_POLY_ORDER]
#define C4 cms_pow_exp_data.poly[7 - EXP_POLY_ORDER]
#define C5 cms_pow_exp_data.poly[8 - EXP_POLY_ORDER]
#define C6 cms_pow_exp_data.poly[9 - EXP_POLY_ORDER]

/* The caller guarantees 2^-54 <= |x| < 512 and positive output. */
CMS_POW_TARGET static inline double
exp_normal(double x, double xtail)
{
  CMS_POW_CONTRACT
  uint64_t ki, idx, top, sbits;
  double kd, z, r, r2, scale, tail, tmp;
  /* exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].  */
  /* x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].  */
  z = InvLn2N * x;
#if TOINT_INTRINSICS
  kd = roundtoint (z);
  ki = converttoint (z);
#elif EXP_USE_TOINT_NARROW
  /* z - kd is in [-0.5-2^-16, 0.5] in all rounding modes.  */
  kd = eval_as_double (z + Shift);
  ki = asuint64 (kd) >> 16;
  kd = (double) (int32_t) ki;
#else
  /* z - kd is in [-1, 1] in non-nearest rounding modes.  */
  kd = eval_as_double (z + Shift);
  ki = asuint64 (kd);
  kd -= Shift;
#endif
  r = x + kd * NegLn2hiN + kd * NegLn2loN;
  /* The code assumes 2^-200 < |xtail| < 2^-8/N.  */
  r += xtail;
  /* 2^(k/N) ~= scale * (1 + tail).  */
  idx = 2 * (ki % N);
  top = ki << (52 - EXP_TABLE_BITS);
  tail = asdouble (T[idx]);
  /* This is only a valid scale when -1023*N < k < 1024*N.  */
  sbits = T[idx + 1] + top;
  /* exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).  */
  /* Evaluation is optimized assuming superscalar pipelined execution.  */
  r2 = r * r;
  /* Without fma the worst case error is 0.25/N ulp larger.  */
  /* Worst case error is less than 0.5+1.11/N+(abs poly error * 2^53) ulp.  */
#if EXP_POLY_ORDER == 4
  tmp = tail + r + r2 * C2 + r * r2 * (C3 + r * C4);
#elif EXP_POLY_ORDER == 5
  tmp = tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);
#elif EXP_POLY_ORDER == 6
  tmp = tail + r + r2 * (0.5 + r * C3) + r2 * r2 * (C4 + r * C5 + r2 * C6);
#endif
  scale = asdouble(sbits);
  return eval_as_double(scale + scale * tmp);
}

/* General-domain wrapper: unsupported inputs and result ranges use libm.
   Arithmetic and coefficient tables inside the normal kernel are unchanged. */
CMS_POW_TARGET static double
PowFMA(double x, double y)
{
  CMS_POW_CONTRACT
  uint64_t ix = asuint64(x), iy = asuint64(y);
  uint32_t topx = ix >> 52, topy = iy >> 52;
  if (topx - 0x001 >= 0x7ff - 0x001 ||
      (topy & 0x7ff) - 0x3be >= 0x43e - 0x3be)
    return pow(x,y);
  double lo;
  double hi = log_inline (ix, &lo);
  double ehi, elo;
#if HAVE_FAST_FMA
  ehi = y * hi;
  elo = y * lo + fma (y, hi, -ehi);
#else
  double yhi = asdouble (iy & -1ULL << 27);
  double ylo = y - yhi;
  double lhi = asdouble (asuint64 (hi) & -1ULL << 27);
  double llo = hi - lhi + lo;
  ehi = yhi * lhi;
  elo = ylo * lhi + y * llo; /* |elo| < |ehi| * 2^-25.  */
#endif
  uint32_t abstop = top12(ehi) & 0x7ff;
  if (abstop - top12(0x1p-54) >= top12(512.0) - top12(0x1p-54))
    return pow(x,y);
  return exp_normal(ehi,elo);
}

/* A local fallback avoids resolving libm's own IFUNC in our resolver. */
static double PowDefault(double x, double y)
{
    return pow(x, y);
}

typedef double (*PowEvaluator)(double, double);

static PowEvaluator ResolvePow(void)
{
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
        return PowFMA;
    return PowDefault;
}

double _cmsPowCore(double x, double y) __attribute__((ifunc("ResolvePow")));

#else

double _cmsPowCore(double x, double y)
{
    return pow(x, y);
}

#endif
