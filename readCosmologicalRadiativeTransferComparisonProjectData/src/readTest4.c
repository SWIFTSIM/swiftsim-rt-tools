/* Read result data for test 4 */

#include "cell.h"
#include "histogram.h"
#include "ioutils.h"

/*! number of bins for histograms */
const int NBINS_HIST = 100;

/*! logarithmic min/max values for histogram */
const float logXHI_min = -2.5f;
const float logXHI_max = 0.f;
const float logT_min = 1.9f;
const float logT_max = 5.0f;

int main(int argc, char **argv) {

  if (argc != 2) {
    printf("Error: Wrong usage.\n");
    printf("Usage: ./readTestX testDataFile.bin\n");
    fflush(stdout);
    abort();
  }

  char *filename = argv[1];

  FILE *fp = io_open_file(filename);

  io_read_header(fp);

  float *xHI = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *T = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));

  io_read_scalar_field(fp, xHI);
  io_read_scalar_field(fp, T);

  /* Get slices */
  io_write_slice(filename, xHI, "xHI", 64);
  io_write_slice(filename, T, "T", 64);

  /* Get histograms */
  int *xHI_hist = NULL;
  get_log_histogram(xHI, &xHI_hist, NBINS_HIST, logXHI_min, logXHI_max);

  int *T_hist = NULL;
  get_log_histogram(T, &T_hist, NBINS_HIST, logT_min, logT_max);

  float *xHII = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  for (int i = 0; i < NCELLS * NCELLS * NCELLS; i++) {
    xHII[i] = (float)(1. - (double)xHI[i]);
  }
  int *xHII_hist = NULL;
  get_log_histogram(xHII, &xHII_hist, NBINS_HIST, logXHI_min, logXHI_max);

  /* Write histograms */
  io_write_histogram(filename, xHI_hist, NBINS_HIST, "logXHI", logXHI_min,
                     logXHI_max);
  io_write_histogram(filename, xHII_hist, NBINS_HIST, "logXHII", logXHI_min,
                     logXHI_max);
  io_write_histogram(filename, T_hist, NBINS_HIST, "logT", logT_min, logT_max);

  /* Cleanup */
  io_check_reached_EOF(fp);

  fclose(fp);
  free(xHI);
  free(xHII);
  free(T);
  free(xHI_hist);
  free(xHII_hist);
  free(T_hist);

  return 0;
}
