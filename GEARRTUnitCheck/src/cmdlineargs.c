#include "cmdlineargs.h"
#include "error.h"
#include <string.h>

/**
 * @brief print the --help.
 */
void print_help_and_exit(void) {

  printf(
      "GEARRTUnitCheck\n"
      "A simple program to check units of initial conditions and simulation\n"
      "parameters for RT simulations by converting units and estimating\n"
      "typical values for gas/radiation quantities, and checking whether\n"
      "these values may lead to problems in terms of precision limits.\n\n"
      "Usage: \n\n"
      "  ./GEARRTUnitCheck --help\n"
      "  ./GEARRTUnitCheck -h\n\n"
      "      Print this message and exit.\n\n"
      "  ./GEARRTUnitCheck <IC_params_file.yml> "
      "<simulation_params_file.yml>\n\n"
      "      Run the check using <IC_params_file.yml> as the parameter \n"
      "      file of the initial conditions, and "
      "<simulation_params_file.yml>\n"
      "      as the parameter file intended to be used for the simulation\n"
      "      run.\n\n"
      "      See swiftsim-rt-tools/GEARRTUnitCheck/README.md for more \n"
      "      details on the parameter files and on how to generate the \n"
      "      IC_params_file.yml\n");

  fflush(stdout);
  exit(0);
}

/**
 * @brief Check that a file exists and can be read.
 *
 * @param filename filename to check.
 */
void verify_file_exists(char filename[MAX_FILENAME_SIZE]) {

  FILE *f = fopen(filename, "r");

  if (f == NULL)
    error("file '%s' not found.\n", filename);

  fclose(f);
}

/**
 * @brief read in cmdline arguments, and check that the cmdline
 * arguments are valid file names. Print help and exit otherwise.
 *
 * @param sim_run_params_filename (out) parameter file to be used in simulation
 * @param IC_params_filename (out) parameter file containing IC data
 */
void read_cmdlineargs(char sim_run_params_filename[MAX_FILENAME_SIZE],
                      char IC_params_filename[MAX_FILENAME_SIZE], int argc,
                      char *argv[]) {

  /* Initialize to NULL. */
  strcpy(sim_run_params_filename, "NULL\0");
  strcpy(IC_params_filename, "NULL\0");

  if (argc <= 1) {
    printf("ERROR: No cmdline args given. Printing help and exiting.\n\n");
    print_help_and_exit();
  }

  for (int arg = 1; arg < argc; arg++) {

    char *newarg = argv[arg];

    /* Help flag raised? */
    if ((strcmp(newarg, "-h") == 0) || strcmp(newarg, "--help") == 0) {
      print_help_and_exit();
    }

    if (strcmp(IC_params_filename, "NULL") == 0) {
      strcpy(IC_params_filename, newarg);
    } else if (strcmp(sim_run_params_filename, "NULL") == 0) {
      strcpy(sim_run_params_filename, newarg);
    } else {
      printf("Unrecognized argument '%s'\n", newarg);
    }
  }

  if (argc != 3) {
    /* Something is wrong. Print help and exit, if you didn't already
     * while checking for provided cmdline args. */
    printf("ERROR: Wrong number of cmdline args given. "
           "Printing help and exiting.\n\n");
    print_help_and_exit();
  }

  verify_file_exists(IC_params_filename);
  verify_file_exists(sim_run_params_filename);
}
