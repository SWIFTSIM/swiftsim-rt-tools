/* Read result data for test 6 */

/* This script assumes that the original files have been renamed
 * as follows: Form the format <codename><outputnumber>.bin and
 * <codename>_add<outputnumber>.bin to <codename><some suffix>
 * and <codename><some suffix>-P2, respectively
 *
 * Example:
 * licorice.bin and licorice_add1.bin to Licorice_1Myr and Licorice_1Myr-P2
 *
 * This allows to take only one cmdline argument, and get the
 * second file by knowing that the suffix '-P2' will be added
 */

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

  float *xHI = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *P = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *T = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *rho = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *mach = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *xHII = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));

  float *xHI_profile = malloc(NCELLS * sizeof(float));
  float *T_profile = malloc(NCELLS * sizeof(float));
  float *P_profile = malloc(NCELLS * sizeof(float));
  float *rho_profile = malloc(NCELLS * sizeof(float));
  float *mach_profile = malloc(NCELLS * sizeof(float));
  float *xHII_profile = malloc(NCELLS * sizeof(float));



  /* Read first file */
  char filename[80] = "\0";
  strcpy(filename, argv[1]);

  FILE *fp = io_open_file(filename);

  io_read_header(fp); io_read_scalar_field(fp, xHI);
  io_read_scalar_field(fp, P);
  io_read_scalar_field(fp, T);
  io_check_reached_EOF(fp);
  fclose(fp);

  /* Read second file */
  char filename2[80] = "\0";
  strcpy(filename2, filename);
  strcat(filename2, "-P2");

  FILE *fp2 = io_open_file(filename2);

  io_read_header(fp2);
  io_read_scalar_field(fp2, rho);
  io_read_scalar_field(fp2, mach);
  io_read_scalar_field(fp2, xHII);
  io_check_reached_EOF(fp2);
  fclose(fp2);


  /* Write outputs now */
  io_write_slice(filename, xHI, "xHI", 64);
  io_write_slice(filename, xHII, "xHII", 64);
  io_write_slice(filename, rho, "n", 64);
  io_write_slice(filename, T, "T", 64);
  io_write_slice(filename, P, "P", 64);
  io_write_slice(filename, mach, "mach", 64);

  extract_line_profile(xHI, xHI_profile, NCELLS);
  extract_line_profile(xHII, xHII_profile, NCELLS);
  extract_line_profile(rho, rho_profile, NCELLS);
  extract_line_profile(T, T_profile, NCELLS);
  extract_line_profile(P, P_profile, NCELLS);
  extract_line_profile(mach, mach_profile, NCELLS);

  io_write_profile(filename, xHI_profile, /*std=*/NULL, NCELLS, "xHI");
  io_write_profile(filename, xHII_profile, /*std=*/NULL, NCELLS, "xHII");
  io_write_profile(filename, rho_profile, /*std=*/NULL, NCELLS, "n");
  io_write_profile(filename, T_profile, /*std=*/NULL, NCELLS, "T");
  io_write_profile(filename, P_profile, /*std=*/NULL, NCELLS, "P");
  io_write_profile(filename, mach_profile, /*std=*/NULL, NCELLS, "mach");


  /* Cleanup */

  free(xHI);
  free(xHI_profile);
  free(T);
  free(T_profile);
  free(P);
  free(P_profile);

  free(xHII);
  free(xHII_profile);
  free(rho);
  free(rho_profile);
  free(mach);
  free(mach_profile);

  return 0;
}
