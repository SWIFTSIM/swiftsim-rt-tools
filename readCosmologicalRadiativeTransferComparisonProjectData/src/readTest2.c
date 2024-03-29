/* Read result data for test 2 */

#include "cell.h"
#include "histogram.h"
#include "ioutils.h"

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

  /* Intended to check data */
  /* for (int z = 0; z < NCELLS; z++) { */
  /*   io_write_slice(filename, xHI, "xHI", z); */
  /*   io_write_slice(filename, T, "T", z); */
  /* } */

  float *xHI_hist = NULL;
  float *xHI_std = NULL;
  get_profile(&xHI_hist, &xHI_std, xHI);
  io_write_profile(filename, xHI_hist, xHI_std, HISTOGRAM_NBINS, "xHI");

  float *T_hist = NULL;
  float *T_std = NULL;
  get_profile(&T_hist, &T_std, T);
  io_write_profile(filename, T_hist, T_std, HISTOGRAM_NBINS, "T");

  io_check_reached_EOF(fp);

  fclose(fp);
  free(xHI);
  free(xHI_hist);
  free(T);
  free(T_hist);

  return 0;
}
