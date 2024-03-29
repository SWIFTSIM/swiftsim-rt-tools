#ifndef READ_PARAMS_H
#define READ_PARAMS_H

void params_read_file(struct swift_params *params, char *param_filename);
void params_read_simulation_params(struct swift_params *params,
                                   struct units *units,
                                   struct simulation_params *simulation_params);
void params_read_ic_params(struct swift_params *params, struct units *units,
                           struct simulation_params *simulation_params,
                           int verbose);

#endif
