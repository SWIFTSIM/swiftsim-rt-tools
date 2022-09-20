/* Read result data for test 5 */

/* This script assumes that the original files have been renamed
 * as follows: Form the format <codename><outputnumber>.bin and
 * <codename>_add<outputnumber>.bin to <codename><some suffix>
 * and <codename><some suffix>-P2, respectively
 *
 * Example:
 * susa1.bin and susa_add1.bin to RSPH_1Myr and RSPH-1Myr-P2
 *
 * This allows to take only one cmdline argument, and get the
 * second file by knowing that the suffix '-P2' will be added
 */

/* Some code outputs use a record length of 2 for whatever Lemmy
 * forsaken reason. So set that up here.
 * These outputs are for the codes:
 * - Zeus
 * - Flash
 * - Licorice
 *
 * For the other codes, you'll need to comment out this line.
 * Or use #define RECORDLEN 1
 */
/* #define RECORDLEN 2 */

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

  /* Read first file */
  char filename[80] = "\0";
  strcpy(filename, argv[1]);

  FILE *fp = io_open_file(filename);

  io_read_header(fp);
  fflush(stdout);

  float *xHI = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *P = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *T = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));

  io_read_scalar_field(fp, xHI);
  io_read_scalar_field(fp, P);
  io_read_scalar_field(fp, T);

  io_write_slice(filename, xHI, "xHI", 0);
  io_write_slice(filename, T, "T", 0);
  io_write_slice(filename, P, "P", 0);

  float *xHI_hist = NULL;
  float *xHI_std = NULL;
  float *T_hist = NULL;
  float *T_std = NULL;
  float *P_hist = NULL;
  float *P_std = NULL;

  get_profile(&xHI_hist, &xHI_std, xHI);
  get_profile(&T_hist, &T_std, T);
  get_profile(&P_hist, &P_std, P);

  io_write_profile(filename, xHI_hist, xHI_std,  NCELLS, "xHI");
  io_write_profile(filename, T_hist, T_std, NCELLS, "T");
  io_write_profile(filename, P_hist, P_std, NCELLS, "P");

  io_check_reached_EOF(fp);

  fclose(fp);
  free(xHI);
  free(xHI_hist);
  free(xHI_std);
  free(T);
  free(T_hist);
  free(T_std);
  free(P);
  free(P_hist);
  free(P_std);



  /* Read second file */
  char filename2[80] = "\0";
  strcpy(filename2, filename);
  strcat(filename2, "-P2");

  FILE *fp2 = io_open_file(filename2);

  io_read_header(fp2);

  float *rho = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *mach = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));
  float *xHII = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));

  io_read_scalar_field(fp2, rho);
  io_read_scalar_field(fp2, mach);
  io_read_scalar_field(fp2, xHII);

  /* Deliberately use filname of first file here */
  io_write_slice(filename, rho, "n", 0);
  io_write_slice(filename, mach, "mach", 0);
  io_write_slice(filename, xHII, "xHII", 0);

  float *rho_hist = NULL;
  float *rho_std = NULL;
  float *mach_hist = NULL;
  float *mach_std = NULL;
  float *xHII_hist = NULL;
  float *xHII_std = NULL;

  get_profile(&rho_hist, &rho_std, rho);
  get_profile(&mach_hist, &mach_std, mach);
  get_profile(&xHII_hist, &xHII_std, xHII);

  io_write_profile(filename, rho_hist, rho_std, NCELLS, "n");
  io_write_profile(filename, mach_hist, mach_std, NCELLS, "mach");
  io_write_profile(filename, xHII_hist, xHII_std,  NCELLS, "xHII");

  io_check_reached_EOF(fp2);

  fclose(fp2);

  free(xHII);
  free(xHII_hist);
  free(xHII_std);
  free(rho);
  free(rho_hist);
  free(rho_std);
  free(mach);
  free(mach_hist);
  free(mach_std);


  return 0;
}
