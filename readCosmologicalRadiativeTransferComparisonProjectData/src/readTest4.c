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

  io_write_slice(filename, xHI, "xHI", 64);
  io_write_slice(filename, T, "T", 64);

  fclose(fp);
  free(xHI);
  free(T);

  return 0;
}
