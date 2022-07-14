#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC   = GEARRT_unit_checks
# source directory
SRCDIR = .
# include directory
IDIR = ../grackle/include
# manual definitions
DEFINES=

# Grackle Dirs
# for 64 bit floats
GRACKLE_INCL =  -I$(GRACKLE_ROOT)/include
GRACKLE_LIBS =  -L$(GRACKLE_ROOT)/lib
CONFIG_BFLOAT= CONFIG_BFLOAT_8

# for 32 bit floats
# GRACKLE_INCL =  -I$(GRACKLE_FLOAT_ROOT)/include
# GRACKLE_LIBS =  -L$(GRACKLE_FLOAT_ROOT)/lib
# CONFIG_BFLOAT= CONFIG_BFLOAT_4

GSL_INCL = -I$(GSL_ROOT)/include
GSL_LIBS = -I$(GSL_ROOT)/lib

#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================

#Compiler to use
CC = gcc 
# additional options
# NOTE: some SWIFT functions need gcc extensions, so -pedantic will act up
OPTIONS = -std=c11 -DH5_USE_16_API -fPIC -D$(CONFIG_BFLOAT)
# optimization
# OPT = -O3
OPT = -Og
# dev options
# DEV =
DEV = -g -Wall -Wextra -fsignaling-nans -Warray-bounds -Werror -fsanitize=address
# Warning: -ffinite-math-only assumes all math will be finite, 
# and won't catch infs and nans

# Additional libraries
LIBS = $(GRACKLE_LIBS)
LIBS += $(GSL_LIBS)
# Additional header files
INCLUDES = -I$(IDIR) -I.
INCLUDES += $(GRACKLE_INCL)
INCLUDES += $(GSL_INCL)
# linking flags
LDFLAGS  = -lm -lgrackle -lgfortran -lgsl

CFLAGS = $(DEFINES) $(OPT) $(OPTIONS) $(DEV) $(INCLUDES) $(LIBS) $(LDFLAGS)




#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(IDIR):$(SRCDIR)


# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
OBJECTS = main.o

HEADERS = constants.h ionization_equilibrium.h mean_molecular_weight.h cross_sections.h \
		  my_grackle_utils.h blackbody.h cross_sections.h photon_interaction_rates.h \
		  parser.h error.h grackle_checks.h grackle_cooling_test.h



#==================================
# RECEPIES
#==================================

all: default

default: $(EXEC)

$(EXEC): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

main.o: main.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ -c $<


#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXEC)
