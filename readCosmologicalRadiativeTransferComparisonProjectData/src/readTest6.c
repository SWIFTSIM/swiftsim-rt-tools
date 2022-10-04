/* Read result data for test 6 */

/* TODO: check comments here */
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
#define RECORDLEN 2

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

  float *xHI_hist = NULL;
  float *xHI_std = NULL;
  float *T_hist = NULL;
  float *T_std = NULL;
  float *P_hist = NULL;
  float *P_std = NULL;

  get_profile(&xHI_hist, &xHI_std, xHI);
  get_profile(&T_hist, &T_std, T);
  get_profile(&P_hist, &P_std, P);
  io_check_reached_EOF(fp);
  fclose(fp);

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

  float *rho_hist = NULL;
  float *rho_std = NULL;
  float *mach_hist = NULL;
  float *mach_std = NULL;
  float *xHII_hist = NULL;
  float *xHII_std = NULL;

  get_profile(&rho_hist, &rho_std, rho);
  get_profile(&mach_hist, &mach_std, mach);
  get_profile(&xHII_hist, &xHII_std, xHII);
  io_check_reached_EOF(fp2);

  /* Write output now */
  /* Generate output filename */
  char outputfile[200] = "\0";
  strcpy(outputfile, filename);
  char addname[80];
  sprintf(addname, "_profiles.dat");
  strcat(outputfile, addname);

  FILE *outfp = fopen(outputfile, "w");
  if (outfp == NULL) {
    printf("couldn't open file '%s'\n", outputfile);
    abort();
  }
  fprintf(outfp,
          "# %10s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
          "xHI [1]", "xHI std", "xHII [1]", "xHII std", "n [cm^-3]", "n std",
          "T [K]", "T std", "P [g/cm/s^2]", "P std", "Mach [1]", "Mach std");
  for (int i = 0; i < HISTOGRAM_NBINS; i++) {
    fprintf(outfp,
            "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e "
            "%12.6e %12.6e %12.6e\n",
            xHI_hist[i], xHI_std[i], xHII_hist[i], xHII_std[i], rho_hist[i],
            rho_std[i], T_hist[i], T_std[i], P_hist[i], P_std[i], mach_hist[i],
            mach_std[i]);
  }

  fclose(outfp);

  printf("written file %s\n", outputfile);

  fclose(fp2);

  free(xHI);
  free(xHI_hist);
  free(xHI_std);
  free(T);
  free(T_hist);
  free(T_std);
  free(P);
  free(P_hist);
  free(P_std);

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
