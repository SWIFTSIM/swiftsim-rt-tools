/* Read result data for test 4 Initial Contition density field */

#include "cell.h"
#include "histogram.h"
#include "ioutils.h"

int main(int argc, char **argv) {

  if (argc != 2) {
    printf("Error: Wrong usage.\n");
    printf("Usage: ./readTest4IC density.bin\n");
    fflush(stdout);
    abort();
  }

  char *filename = argv[1];

  int header, footer;
  float redshift = -1.f;
  float dummy = -1.f;
  float buffer;
  float *density_field = malloc(NCELLS * NCELLS * NCELLS * sizeof(float));

  /* Read in data */
  FILE *fp = io_open_file(filename);

  header = 0;
  footer = 0;
  fread(&header, INTSIZE, 1, fp);
  fread(&redshift, FLOATSIZE, 1, fp);
  fread(&footer, INTSIZE, 1, fp);
  check_record(header, 4);
  check_record(footer, 4);

  header = 0;
  footer = 0;
  fread(&header, INTSIZE, 1, fp);
  fread(&dummy, FLOATSIZE, 1, fp);
  fread(&footer, INTSIZE, 1, fp);
  check_record(header, 4);
  check_record(footer, 4);

  /* printf("redshift %.3f dummy %.3f\n", redshift, dummy); */

  for (int i = 0; i < NCELLS * NCELLS * NCELLS; i++)
    density_field[i] = 0.f;

  for (int k = 0; k < 128; k++) {
    header = 0;
    footer = 0;
    fread(&header, INTSIZE, 1, fp);
    check_record(header, 4 * 128 * 128);
    for (int j = 0; j < 128; j++) {
      for (int i = 0; i < 128; i++) {
        fread(&buffer, INTSIZE, 1, fp);
        int ind = get_array_index(i, j, k);
        density_field[ind] = buffer;
      }
    }
    fread(&footer, INTSIZE, 1, fp);
    check_record(footer, 4 * 128 * 128);
  }

  /* For debugging checks */
  /* for (int z = 0; z < 128; z++) { */
  /*   io_write_slice(filename, density_field, "density", z); */
  /* } */

  /* Now write output */
  char outputfile[80] = "Iliev4DensityIC.dat";
  FILE *ofp = fopen(outputfile, "w");
  if (ofp == NULL) {
    printf("Error opening output file '%s'\n", outputfile);
    abort();
  }

  fprintf(ofp, "# cell index i,j,k, density in cm^-3\n");
  fprintf(ofp, "# redshift %.6g\n", redshift);
  for (int i = 0; i < 128; i++) {
    for (int j = 0; j < 128; j++) {
      for (int k = 0; k < 128; k++) {
        int ind = get_array_index(i, j, k);
        fprintf(ofp, "%4d %4d %4d %.6e\n", i, j, k, density_field[ind]);
      }
    }
  }

  io_check_reached_EOF(fp);

  /* cleanup */
  fclose(fp);
  free(density_field);

  return 0;
}
