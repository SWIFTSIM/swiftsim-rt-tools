#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC   = GEARRT_unit_checks
# source directory
SRCDIR = ../src
# manual definitions
DEFINES=

# define grackle root directory, if the variable isn't defined already.
# GRACKLE_SWIFT_ROOT=/path/to/your/grackle/install

# Grackle Dirs
# for 64 bit floats
GRACKLE_INCL =  -I$(GRACKLE_SWIFT_ROOT)/include
GRACKLE_LIBS =  -L$(GRACKLE_SWIFT_ROOT)/lib
CONFIG_BFLOAT= CONFIG_BFLOAT_8

# for 32 bit floats
# GRACKLE_INCL =  -I$(GRACKLE_FLOAT_ROOT)/include
# GRACKLE_LIBS =  -L$(GRACKLE_FLOAT_ROOT)/lib
# CONFIG_BFLOAT= CONFIG_BFLOAT_4

# GSL dirs
# GSL_ROOT=/path/to/your/gsl/install
GSL_INCL = -I$(GSL_ROOT)/include
GSL_LIBS = -L$(GSL_ROOT)/lib

# hdf5 dirs
# HDF5_ROOT=/path/to/your/hdf5/install
HDF5_INCL += -I$(HDF5_ROOT)/include
HDF5_LIBS=-L$(HDF5_ROOT)/lib

#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================

#Compiler to use
# gcc
CC = gcc
CCLIBS=-lgfortran
# intel oneAPI
# CC = icx
# CCLIBS=-lifcore -lifport
#
# additional options
# NOTE: some SWIFT functions need gcc extensions, so -pedantic will act up
OPTIONS = -std=c11 -DH5_USE_16_API -fPIC -D$(CONFIG_BFLOAT)
# optimization
# OPT = -O3
OPT = -Og
# dev options
# DEV =
DEV = -g -Wall -Wextra -fsignaling-nans -Warray-bounds 
# -Werror
#-fsanitize=address
# Warning: -ffinite-math-only assumes all math will be finite, 
# and won't catch infs and nans

# Additional libraries
LDFLAGS = $(GRACKLE_LIBS)
LDFLAGS += $(GSL_LIBS)
LDFLAGS += $(HDF5_LIBS)
# Additional header files
INCLUDES = -I. -I../include/
INCLUDES += -I../../grackle/include
INCLUDES += $(GRACKLE_INCL)
INCLUDES += $(GSL_INCL)
INCLUDES += -I$(HDF5_INCL)

# library flags
LIBS  = -lgrackle $(CCLIBS) -lhdf5 -lgsl -lgslcblas -lm



CFLAGS = $(DEFINES) $(OPT) $(OPTIONS) $(DEV) $(INCLUDES)




#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(SRCDIR):../include:../../grackle/include

# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
OBJECTS = main.o read_params.o units.o simulation_params.o

HEADERS = constants.h ionization_equilibrium.h mean_molecular_weight.h cross_sections.h \
		  my_grackle_utils.h blackbody.h cross_sections.h photon_interaction_rates.h \
		  parser.h error.h grackle_checks.h grackle_cooling_test.h grackle_heating_test.h \
		  read_params.h validity_check_macros.h units.h simulation_params.h rt_ngroups.h \
		  conversions.h checks.h cmdlineargs.h

OBJECTS = main.o read_params.o checks.o cmdlineargs.o conversions.o simulation_params.o units.o parser.o 


#==================================
# RECEPIES
#==================================

all: default

default: $(EXEC)

$(EXEC): $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $(OBJECTS) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $(LIBS)

main.o: main.c $(HEADERS)
	$(CC) -c $< $(CFLAGS) -o $@

read_params.o: read_params.c read_params.h error.h parser.h simulation_params.h
	$(CC) -c $< $(CFLAGS) -o $@

units.o: units.c units.h error.h validity_check_macros.h
	$(CC) -c $< $(CFLAGS) -o $@

simulation_params.o: simulation_params.c simulation_params.h conversions.h error.h validity_check_macros.h constants.h 
	$(CC) -c $< $(CFLAGS) -o $@

conversions.o: conversions.c conversions.h error.h validity_check_macros.h constants.h
	$(CC) -c $< $(CFLAGS) -o $@

checks.o: checks.c checks.h error.h validity_check_macros.h simulation_params.h units.h
	$(CC) -c $< $(CFLAGS) -o $@

parser.o: parser.c parser.h error.h
	$(CC) -c $< $(CFLAGS) -o $@

cmdlineargs.o: cmdlineargs.c cmdlineargs.h error.h
	$(CC) -c $< $(CFLAGS) -o $@


#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXEC)
