#ifndef VALIDITY_CHECK_MACROS_H
#define VALIDITY_CHECK_MACROS_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @file validity_check_macros.h
 *
 * Macros to check validity of floatds, doubles, etc.
 * Mainly so I can directly print out variables + variable names,
 * including filenames for tracebacks etc. Convenience functions,
 * basically.
 */

/*! Check that the value is a valid float.
 * Assume the argument given is positive and of type double. */
#define check_valid_float(v, check_grackle_limits)                             \
  ({                                                                           \
    if (v < 0.) {                                                              \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; negative: %.6e\n",      \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v > FLT_MAX) {                                                         \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; > FLT_MAX: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if ((v != 0.) && (v < FLT_MIN)) {                                          \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; < FLT_MIN: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    float vfloat = (float)v;                                                   \
    if (isinf(vfloat) || isnan(vfloat)) {                                      \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; nan/inf %.6e\n",        \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v != 0. && (fabs(v) > 1e30 || fabs(v) < 1e-30)) {                      \
      fflush(stdout);                                                          \
      fprintf(stdout, "WARNING: %s:%s:%d: " #v " has large exponent: %.6e\n",  \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      warnings++;                                                              \
    }                                                                          \
    if (check_grackle_limits) {                                                \
      if (v < 1.e-20) {                                                        \
        fflush(stdout);                                                        \
        fprintf(stderr,                                                        \
                "%s:%s:%d: " #v " below grackle TINY_NUMBER 1e-20: %.6e\n",    \
                __FILE__, __FUNCTION__, __LINE__, v);                          \
        abort();                                                               \
      }                                                                        \
    }                                                                          \
  })

#define FLOAT_TOLERANCE 1e-5
#define check_floats_equal(a, b)                                               \
  ({                                                                           \
    if (a == 0. && b == 0.) {                                                  \
    } else if (a == 0. && fabsf(b) > FLOAT_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else if (b == 0. && fabsf(a) > FLOAT_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else {                                                                   \
      if (1.f - fabsf(a / b) > FLOAT_TOLERANCE) {                              \
        error(#a " and " #b " are not equal: %.6e %.6e a/b=%.3e", a, b,        \
              a / b);                                                          \
      }                                                                        \
    }                                                                          \
  })

/*! Check that the value is a valid double.
 * Assume the argument given is positive and of type double. */
#define check_valid_double(v, check_grackle_limits)                            \
  ({                                                                           \
    if (v < 0.) {                                                              \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid double; negative: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (isinf(v) || isnan(v)) {                                                \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid double; nan/inf %.6e\n",       \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v > DBL_MAX) {                                                         \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid double; > DBL_MAX: %.6e\n",    \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if ((v != 0.) && (v < DBL_MIN)) {                                          \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; < DBL_MIN: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v != 0. && (fabs(v) > 1e290 || fabs(v) < 1e-290)) {                    \
      fflush(stdout);                                                          \
      fprintf(stdout, "WARNING: %s:%s:%d: " #v " has large exponent: %.6e\n",  \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      warnings++;                                                              \
    }                                                                          \
    if (check_grackle_limits) {                                                \
      if (v < 1.e-20) {                                                        \
        fflush(stdout);                                                        \
        fprintf(stderr,                                                        \
                "%s:%s:%d: " #v " below grackle TINY_NUMBER 1e-20: %.6e\n",    \
                __FILE__, __FUNCTION__, __LINE__, v);                          \
        abort();                                                               \
      }                                                                        \
    }                                                                          \
  })

#define DOUBLE_TOLERANCE 1e-15
#define check_doubles_equal(a, b)                                              \
  ({                                                                           \
    if (a == 0. && b == 0.) {                                                  \
    } else if (a == 0. && fabs(b) > DOUBLE_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else if (b == 0. && fabs(a) > DOUBLE_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else {                                                                   \
      if (1.f - fabs(a / b) > DOUBLE_TOLERANCE) {                              \
        error(#a " and " #b " are not equal: %.6e %.6e a/b=%.3e", a, b,        \
              a / b);                                                          \
      }                                                                        \
    }                                                                          \
  })

#endif
