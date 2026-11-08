/*
 * Configuration for math routines.
 *
 * Copyright (c) 2017-2026, Arm Limited.
 * SPDX-License-Identifier: MIT OR Apache-2.0 WITH LLVM-exception
 */

/* Adapted for Little CMS from Arm optimized-routines commit
   20859d3fd9e797dde592c489c8ef7d383f778620. MIT license option selected.
   See LICENSE and README.txt. */

#ifndef CMS_POW_CONFIG_H
#define CMS_POW_CONFIG_H
#include <float.h>
#include <math.h>
#include <stdint.h>

#define HAVE_FAST_FMA 1
#define TOINT_INTRINSICS 0
#define ALIGN(X) __attribute__((aligned(X)))

static inline uint64_t
asuint64 (double f)
{
  union
  {
    double f;
    uint64_t i;
  } u = {f};
  return u.i;
}

static inline double
asdouble (uint64_t i)
{
  union
  {
    uint64_t i;
    double f;
  } u = {i};
  return u.f;
}

static inline double
eval_as_double (double x)
{
  return x;
}

#define EXP_TABLE_BITS 7
#define EXP_POLY_ORDER 5
/* Use polynomial that is optimized for a wider input range.  This may be
   needed for good precision in non-nearest rounding and !TOINT_INTRINSICS.  */
#define EXP_POLY_WIDE 0
/* Use close to nearest rounding toint when !TOINT_INTRINSICS.  This may be
   needed for good precision in non-nearest rouning and !EXP_POLY_WIDE.  */
#define EXP_USE_TOINT_NARROW 0
#define EXP2_POLY_ORDER 5
#define EXP2_POLY_WIDE 0
/* Wider exp10 polynomial necessary for good precision in non-nearest rounding
   and !TOINT_INTRINSICS.  */
#define EXP10_POLY_WIDE 0
struct exp_data
{
  double invln2N;
  double negln2hiN;
  double negln2loN;
  double poly[4]; /* Last four coefficients.  */
  double shift;

  double exp2_shift;
  double exp2_poly[EXP2_POLY_ORDER];

  double neglog10_2hiN;
  double neglog10_2loN;
  double exp10_poly[5];
  uint64_t tab[2*(1 << EXP_TABLE_BITS)];
  double invlog10_2N;
};

#define POW_LOG_TABLE_BITS 7
#define POW_LOG_POLY_ORDER 8
struct pow_log_data
{
  double ln2hi;
  double ln2lo;
  double poly[POW_LOG_POLY_ORDER - 1]; /* First coefficient is 1.  */
  /* Note: the pad field is unused, but allows slightly faster indexing.  */
  struct {double invc, pad, logc, logctail;} tab[1 << POW_LOG_TABLE_BITS];
};


#endif
