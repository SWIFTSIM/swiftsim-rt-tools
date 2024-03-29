#ifndef IOUTILS_H
#define IOUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // stat

#include "cell.h"

#define INTSIZE 4
#define FLOATSIZE 4

/* Allow specific tests to modify how long you assume
 * the fortran record length is */
#ifndef RECORDLEN
#define RECORDLEN 1
#endif

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
    fprintf(stderr, "error while opening file %s\n", filename);
    fflush(stderr);
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

  int header[RECORDLEN];
  int footer[RECORDLEN];
  for (int i = 0; i < RECORDLEN; i++) {
    header[i] = 0;
    footer[i] = 0;
  }
  int gridDim[3] = {0, 0, 0};

  fread(header, INTSIZE, RECORDLEN, fp);
  fread(gridDim, INTSIZE, 3, fp);
  fread(footer, INTSIZE, RECORDLEN, fp);

  /* If != 12, we abort later. Print me some info to screen. */
  if (header[0] != 12 || footer[0] != 12) {
    printf("Debug output:");
    for (int i = 0; i < RECORDLEN; i++)
      printf("header[%d] = %d\n", i, header[i]);

    printf("griddim %d %d %d\n", gridDim[0], gridDim[1], gridDim[2]);

    for (int i = 0; i < RECORDLEN; i++)
      printf("footer[%d] = %d\n", i, footer[i]);
  }

  /* Header and footer shall be enclosed by 12 */
  check_record(header[0], 12);
  check_record(footer[0], 12);

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

  int header[RECORDLEN];
  int footer[RECORDLEN];
  for (int i = 0; i < RECORDLEN; i++) {
    header[i] = 0;
    footer[i] = 0;
  }
  const int expect = 4 * NCELLS * NCELLS * NCELLS;

  fread(header, INTSIZE, RECORDLEN, fp);
  check_record(header[0], expect);
  fread(buffer, FLOATSIZE, NCELLS * NCELLS * NCELLS, fp);
  fread(footer, INTSIZE, RECORDLEN, fp);
  check_record(footer[0], expect);
}

/**
 * write a 2D slice of the data.
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
  if (fp == NULL) {
    printf("couldn't open file '%s'\n", outputfile);
    abort();
  }
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
  if (fp == NULL) {
    printf("couldn't open file '%s'\n", outputfile);
    abort();
  }
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
 * write a histogram to file.
 *
 * @param srcfilename filename that was read in to obtain data
 * @param hist the histogram
 * @param n number of elements in profile array
 * @param descriptor additional descriptor to add to output filename
 * @param minval lower threshold for histogram values
 * @param maxval upper threshold for histogram values
 */
void io_write_histogram(char *srcfilename, int *hist, int n, char *descriptor,
                        float minval, float maxval) {

  /* Generate output filename */
  char outputfile[200] = "\0";
  strcpy(outputfile, srcfilename);
  char addname[80];
  sprintf(addname, "_histogram_%s.dat", descriptor);
  strcat(outputfile, addname);

  float dx = (maxval - minval) / (float)n;

  FILE *fp = fopen(outputfile, "w");
  if (fp == NULL) {
    printf("couldn't open file '%s'\n", outputfile);
    abort();
  }
  for (int i = 0; i < n; i++) {
    fprintf(fp, "%.6e, %6d\n", minval + (i + 0.5) * dx, hist[i]);
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
