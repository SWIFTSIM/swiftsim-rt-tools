#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "cell.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define HISTOGRAM_NBINS 128
#define HISTOGRAM_BIN_WIDTH 1.

#define check_is_NULL(arr)                                                     \
  ({                                                                           \
    if (*arr != NULL) {                                                        \
      printf(#arr "is not NULL\n");                                            \
      abort();                                                                 \
    }                                                                          \
  })

/**
 * Generate a histogram from given data.
 * Assumes data is 3D array.
 * Creates profiles starting at lower left corner of the array
 * (0, 0, 0).
 *
 * @param data the data to histogram
 * @param hist (return) the binned data
 * @param count (return) the event counts in each bin
 *
 * */
void get_histogram(float *data, float **hist, int **count) {

  float *h = malloc(HISTOGRAM_NBINS * sizeof(float));
  int *c = malloc(HISTOGRAM_NBINS * sizeof(int));

  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    h[i] = 0.f;
    c[i] = 0;
  }

  const float max_range = HISTOGRAM_BIN_WIDTH * HISTOGRAM_NBINS;

  for (int i = 0; i < NCELLS; i++) {
    float x = ((float)i + 0.5f) * CELL_WIDTH;
    float x2 = x * x;
    for (int j = 0; j < NCELLS; j++) {
      float y = ((float)j + 0.5f) * CELL_WIDTH;
      float y2 = y * y;
      for (int k = 0; k < NCELLS; k++) {
        float z = ((float)k + 0.5f) * CELL_WIDTH;
        float z2 = z * z;

        float d = sqrtf(x2 + y2 + z2);

        /* don't use anything further away than the box size */
        if (d > max_range)
          continue;

        int index = floor(d / HISTOGRAM_BIN_WIDTH);

        h[index] += data[get_array_index(i, j, k)];
        c[index] += 1;
      }
    }
  }

  check_is_NULL(hist);
  check_is_NULL(count);
  *hist = h;
  *count = c;
}

/**
 * Get the standard deviation of the previously profiled
 * data. The profiles should be spherical averaged values
 * starting from the lower left corner (0, 0, 0).
 *
 * @param data the data to histogram
 * @param mean the mean of the binned data
 * @param count the number of samples in any bin
 * @param std (return) standard deviation of the data
 *
 * */
void get_std(float *data, float *mean, int *count, float **std) {

  float *s = malloc(HISTOGRAM_NBINS * sizeof(float));
  int *c = malloc(HISTOGRAM_NBINS * sizeof(int));

  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    s[i] = 0.f;
    c[i] = 0;
  }

  const float max_range = HISTOGRAM_BIN_WIDTH * HISTOGRAM_NBINS;

  for (int i = 0; i < NCELLS; i++) {
    float x = ((float)i + 0.5f) * CELL_WIDTH;
    float x2 = x * x;
    for (int j = 0; j < NCELLS; j++) {
      float y = ((float)j + 0.5f) * CELL_WIDTH;
      float y2 = y * y;
      for (int k = 0; k < NCELLS; k++) {
        float z = ((float)k + 0.5f) * CELL_WIDTH;
        float z2 = z * z;

        float d = sqrtf(x2 + y2 + z2);

        /* don't use anything further away than the box size */
        if (d > max_range)
          continue;

        int index = floor(d / HISTOGRAM_BIN_WIDTH);
        float val = data[get_array_index(i, j, k)];

        float temp1 = (val - mean[index]);
        float temp2 = temp1 * temp1;

        s[index] += temp2;
        c[index] += 1;
      }
    }
  }

  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    /* Sanity check: Make sure we obtained the same count */
    if (c[i] != count[i]) {
      printf("Inconsistent counts: got=%d had=%d", c[i], count[i]);
      abort();
    }

    if (c[i] == 0)
      continue;
    if (c[i] == 1)
      continue;

    /* Make actual STD out of sum of squares */
    float temp = 1.f / ((float)c[i] - 1.f) * s[i];
    temp = sqrtf(temp);
    s[i] = temp;
  }

  check_is_NULL(std);
  *std = s;
  free(c);
}

/**
 * Get a profile
 *
 * @param profile (return) the pointer to the profile of the data
 * @param std (return) the pointer to the standard deviation of the data in the
 *profile
 * @param data the data to make a profile of
 *
 **/
void get_profile(float **profile, float **std, float *data) {

  float *hist = NULL;
  int *count = NULL;
  get_histogram(data, &hist, &count);

  float *prof = malloc(HISTOGRAM_NBINS * sizeof(float));
  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    prof[i] = 0.f;
  }

  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    float res = NAN;
    if (count[i] > 0)
      res = hist[i] / (float)count[i];
    prof[i] = res;
  }

  get_std(data, prof, count, std);

  check_is_NULL(profile);
  *profile = prof;

  free(hist);
  free(count);
}

#endif
