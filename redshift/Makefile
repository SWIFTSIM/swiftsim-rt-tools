#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC   = analytic
# source directory
SRCDIR = .
# include directory
IDIR = ./include
# manual definitions
DEFINES=

# Depending on the compiler, you need to link a different fortran library
# to successfully link grackle.
# gcc: -lgfortran
# intel: -lifcore

# for gcc 
# GRACKLE_FORTRAN_LIB=-lgfortran
# for intel
GRACKLE_FORTRAN_LIB=-lifcore -lifport

# define grackle root directory, if the variable isn't defined already.
# GRACKLE_SWIFT_ROOT=/path/to/your/grackle/install
# Grackle Dirs
# for 64 bit floats
GRACKLE_INCL =  -I$(GRACKLE_SWIFT_ROOT)/include
GRACKLE_LD =  -L$(GRACKLE_SWIFT_ROOT)/lib
GRACKLE_LIBS = -lgrackle $(GRACKLE_FORTRAN_LIB)
CONFIG_BFLOAT= CONFIG_BFLOAT_8

# for 32 bit floats
# GRACKLE_INCL =  -I$(GRACKLE_FLOAT_ROOT)/include
# GRACKLE_LD =  -L$(GRACKLE_FLOAT_ROOT)/lib
# GRACKLE_LIBS = -lgrackle -lgfortran
# CONFIG_BFLOAT= CONFIG_BFLOAT_4


# GSL dirs
# GSL_ROOT=/path/to/your/gsl/install
GSL_INCL = -I$(GSL_ROOT)/include
GSL_LD = -L$(GSL_ROOT)/lib
GSL_LIBS = -lblas -lgsl -lgslcblas


#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================

#Compiler to use
CC = gcc 
# additional options
# -pedantic # avoid complaints about error macro
OPTIONS = -std=c11 -DH5_USE_16_API -fPIC -D$(CONFIG_BFLOAT)
# optimization
# OPT = -O3
OPT = -Og
# dev options
# DEV =
DEV = -g -Wall -Wextra -fsignaling-nans -Warray-bounds -Werror
# Warning: -ffinite-math-only assumes all math will be finite, 
# and won't catch infs and nans

# linking flags
LDFLAGS = $(GRACKLE_LD)
LDFLAGS += $(GSL_LD)

# Additional header files
INCLUDES = -I$(IDIR) -I.
INCLUDES += $(GRACKLE_INCL)
INCLUDES += $(GSL_INCL)

# Additional libraries
LIBS  = $(GRACKLE_LIBS) $(GSL_LIBS) -lm

CFLAGS = $(DEFINES) $(OPT) $(OPTIONS) $(DEV) $(INCLUDES)




#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(IDIR):$(SRCDIR)


# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
OBJECTS = main_analytic.o main_finlator.o makeIC.o

HEADERS = constants.h cosmology.h redshift.h



#==================================
# RECEPIES
#==================================


default: $(EXEC) makeIC finlator

$(EXEC): $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $(CFLAGS) $(INCLUDES) $(LDFLAGS) $< $(LIBS) 

makeIC: IC.o $(HEADERS)
	$(CC) -o $@ $(CFLAGS) $(INCLUDES) $(LDFLAGS) $< $(LIBS) 

finlator: main_finlator.o $(HEADERS)
	$(CC) -o $@ $(CFLAGS) $(INCLUDES) $(LDFLAGS) $< $(LIBS) 

main_analytic.o: main_analytic.c $(HEADERS)
	$(CC) -c $< $(CFLAGS) $(INCLUDES) -o $@

main_finlator.o: main_finlator.c $(HEADERS)
	$(CC) -c $< $(CFLAGS) $(INCLUDES) -o $@

IC.o: makeIC.c $(HEADERS)
	$(CC) -c $< $(CFLAGS) $(INCLUDES) -o $@

#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXEC) makeIC finlator
	rm -f IC.o makeIC
