#ifndef IOUTILS_H
#define IOUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // stat

#include "cell.h"

#define INTSIZE 4
#define FLOATSIZE 4

#define check_record(read, expect)                                             \
  ({                                                                           \
    if (read != expect) {                                                      \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: Error: wrong " #read "=%d expect %d\n",       \
              __FILE__, __FUNCTION__, __LINE__, read, expect);                 \
      abort();                                                                 \
    }                                                                          \
  })

/**
 * Open up a file pointer and do additional checks.
 *
 * @param filename name of file to open.
 * */
FILE *io_open_file(char *filename) {

  FILE *fp = fopen(filename, "rb");
  if (fp == NULL) {
    printf("error while opening file %s", filename);
    abort();
  }

  return fp;
}

/**
 * Read the header of the test binary file
 *
 * @param fp file pointer to check
 * */
void io_read_header(FILE *fp) {

  int header = 0, footer = 0;
  int gridDim[3] = {0, 0, 0};

  header = 0;
  footer = 0;
  fread(&header, INTSIZE, 1, fp);
  fread(gridDim, INTSIZE, 3, fp);
  fread(&footer, INTSIZE, 1, fp);
  /* Header and footer shall be enclosed by 12 */
  check_record(header, 12);
  check_record(footer, 12);

  if (gridDim[0] != NCELLS || gridDim[1] != NCELLS || gridDim[2] != NCELLS) {
    printf(
        "Got griddim %d x %d x %d, should be %d x %d x %d; something's wrong\n",
        gridDim[0], gridDim[1], gridDim[2], NCELLS, NCELLS, NCELLS);
    abort();
  }
}

/**
 * Read a 3D array of float scalars from binary file
 *
 * @param fp file pointer to check
 * @param buffer (return) where to write read in data into
 * */
void io_read_scalar_field(FILE *fp, float *buffer) {

  int header = 0, footer = 0;
  const int expect = 4 * NCELLS * NCELLS * NCELLS;

  fread(&header, INTSIZE, 1, fp);
  check_record(header, expect);
  fread(buffer, FLOATSIZE, NCELLS * NCELLS * NCELLS, fp);
  fread(&footer, INTSIZE, 1, fp);
  check_record(footer, expect);
}

/**
 * write a 2D slice of the data.
 * Assumes data is 3D array with size 128 * 128 * 128
 *
 * @param srcfilename filename that was read in to obtain data
 * @param data the data
 * @param descriptor additional descriptor to add to output filename
 * @param z z-index to write out
 */
void io_write_slice(char *srcfilename, float *data, char *descriptor, int z) {

  /* Generate output filename */
  char outputfile[200] = "\0";
  strcpy(outputfile, srcfilename);
  char addname[80];
  sprintf(addname, "_slice_%s_z=%d.dat", descriptor, z);
  strcat(outputfile, addname);

  FILE *fp = fopen(outputfile, "w");
  for (int i = 0; i < NCELLS; i++) {
    for (int j = 0; j < NCELLS; j++) {
      int ind = get_array_index(i, j, z);
      fprintf(fp, "%.5e", data[ind]);
      if (j < NCELLS - 1)
        fprintf(fp, ", ");
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  printf("written file %s\n", outputfile);
}

/**
 * write a profile to file.
 * Assumes data is 3D array with size 128 * 128 * 128
 *
 * @param srcfilename filename that was read in to obtain data
 * @param profile the data
 * @param std the standard deviation of the data. May be NULL
 * @param n number of elements in profile array
 * @param descriptor additional descriptor to add to output filename
 */
void io_write_profile(char *srcfilename, float *profile, float *std, int n,
                      char *descriptor) {

  /* Generate output filename */
  char outputfile[200] = "\0";
  strcpy(outputfile, srcfilename);
  char addname[80];
  sprintf(addname, "_profile_%s.dat", descriptor);
  strcat(outputfile, addname);

  FILE *fp = fopen(outputfile, "w");
  for (int i = 0; i < n; i++) {
    if (std != NULL) {
      fprintf(fp, "%.6e, %.6e\n", profile[i], std[i]);
    } else {
      fprintf(fp, "%.6e\n", profile[i]);
    }
  }

  fclose(fp);

  printf("written file %s\n", outputfile);
}

/**
 * Check that you reached the End of File.
 *
 * @param fp file pointer to check
 */
void io_check_reached_EOF(FILE *fp) {

  getc(fp);
  if (!feof(fp)) {
    printf("Error: Didn't reach EOF\n");
    abort();
  }
}

#endif
