#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC   = iliev_test
# source directory
SRCDIR = .
# include directory
IDIR = ../include
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

#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================

#Compiler to use
CC = gcc 
# additional options
OPTIONS = -pedantic -std=c11 -DH5_USE_16_API -fPIC -D$(CONFIG_BFLOAT)
# optimization
# OPT = -O3
OPT = -Og
# dev options
# DEV =
DEV = -g -Wall -Wextra -fsignaling-nans -Warray-bounds -Werror -fsanitize=address
# Warning: -ffinite-math-only assumes all math will be finite, 
# and won't catch infs and nans

# linking flags
LDFLAGS = $(GRACKLE_LIBS)
LDFLAGS += $(GSL_LIBS)

# Additional header files
INCLUDES = -I$(IDIR) -I.
INCLUDES += $(GRACKLE_INCL)
INCLUDES += $(GSL_INCL)

# library flags
LIBS  = -lgrackle -lgfortran -lgsl -lblas -lm



CFLAGS = $(DEFINES) $(OPT) $(OPTIONS) $(DEV) $(INCLUDES)




#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(IDIR):$(SRCDIR)


# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
OBJECTS = main.o

HEADERS = constants.h ionization_equilibrium.h mean_molecular_weight.h cross_sections.h \
		  my_grackle_utils.h blackbody.h cross_sections.h photon_interaction_rates.h



#==================================
# RECEPIES
#==================================


default: $(EXEC)

$(EXEC): $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $(CFLAGS) $(INCLUDES) $(LDFLAGS) $< $(LIBS) 

main.o: main.c $(HEADERS)
	$(CC) -c $< $(CFLAGS) $(INCLUDES) -o $@



#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXEC)
