#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cell.h"

#define HISTOGRAM_NBINS 128
#define HISTOGRAM_BIN_WIDTH 1.


#define check_is_NULL(arr)            \
  ({                                  \
    if (*arr != NULL) {               \
      printf(#arr "is not NULL\n");   \
      abort();                        \
     }                                \
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

void get_histogram(float *data, float** hist, int** count){

  float *h = malloc(HISTOGRAM_NBINS * sizeof(float));
  int *c = malloc(HISTOGRAM_NBINS * sizeof(int));

  for (int i = 0; i < HISTOGRAM_NBINS; i++){
    h[i] = 0.f;
    c[i] = 0;
  }

  const float max_range = HISTOGRAM_BIN_WIDTH * HISTOGRAM_NBINS;

  for (int i = 0; i < NCELLS; i++){
    float x = ((float) i + 0.5f) * CELL_WIDTH;
    float x2 = x * x;
    for (int j = 0; j < NCELLS; j++) {
      float y = ((float) j + 0.5f) * CELL_WIDTH;
      float y2 = y * y;
      for (int k = 0; k < NCELLS; k++){
        float z = ((float) k + 0.5f) * CELL_WIDTH;
        float z2 = z * z;

        float d = sqrtf(x2 + y2 + z2);
        
        /* don't use anything further away than the box size */
        if (d > max_range) continue;

        int index = floor(d / HISTOGRAM_BIN_WIDTH);

        h[index] += data[get_array_index(i,j,k)];
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
 * Get a profile
 *
 * @param profile (return) the pointer to the profile of the data
 * @param data the data to make a profile of
 *
 **/
void get_profile(float **profile, float *data){

  float *hist = NULL;
  int *count = NULL;
  get_histogram(data, &hist, &count);

  float *prof = malloc(HISTOGRAM_NBINS * sizeof(float));
  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    prof[i] = 0.f;
  }

  for (int i = 0; i < HISTOGRAM_NBINS; i++){
    float res = NAN;
    if (count[i] > 0) res = hist[i] / (float) count[i];
    prof[i] = res;
  }

  check_is_NULL(profile);
  *profile = prof;
}

#endif
