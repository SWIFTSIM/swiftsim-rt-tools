/* Read result data for test 3 */

#include "cell.h"
#include "histogram.h"
#include "ioutils.h"

/**
 * Extract the profile along the axis of symmetry
 *
 * @param data data to extract profile from
 * @param profile (return) the returned profile
 * @param n size of data (assumed n x n x n) and profile (n)
 * @param
 */
void extract_line_profile(float *data, float *profile, int n) {

  for (int i = 0; i < n; i++)
    profile[i] = 0.f;

  for (int i = 0; i < n; i++) {

    int check = 0;
    for (int y = 63; y <= 64; y++) {
      for (int z = 63; z <= 64; z++) {
        int ind = get_array_index(i, y, z);
        profile[i] += data[ind];
        check += 1;
      }
    }

    profile[i] *= 0.25;
  }
}

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

  float *xHI_profile = malloc(NCELLS * sizeof(float));
  float *T_profile = malloc(NCELLS * sizeof(float));

  extract_line_profile(xHI, xHI_profile, NCELLS);
  extract_line_profile(T, T_profile, NCELLS);

  io_write_profile(filename, xHI_profile, /*std=*/NULL, NCELLS, "xHI");
  io_write_profile(filename, T_profile, /*std=*/NULL, NCELLS, "T");

  io_check_reached_EOF(fp);

  fclose(fp);
  free(xHI);
  free(T);
  free(xHI_profile);
  free(T_profile);

  return 0;
}
