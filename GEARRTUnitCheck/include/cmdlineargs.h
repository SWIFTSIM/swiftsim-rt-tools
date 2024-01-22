#ifndef CMDLINEARGS_H
#define CMDLINEARGS_H

#define MAX_FILENAME_SIZE 200

void print_help_and_exit(void);
void verify_file_exists(char filename[MAX_FILENAME_SIZE]);
void read_cmdlineargs(char sim_run_params_filename[MAX_FILENAME_SIZE], char IC_params_filename[MAX_FILENAME_SIZE], int argc, char* argv[]);

#endif
