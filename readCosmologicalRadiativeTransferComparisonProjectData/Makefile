#==================================
# COMPILATION TIME PARAMETERS
#==================================

# source directory
SRCDIR = src
# manual definitions
DEFINES=

#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================

#Compiler to use
CC = gcc 
# additional options
# NOTE: some SWIFT functions need gcc extensions, so -pedantic will act up
OPTIONS = -std=c11
# optimization
# OPT = -O3
OPT = -Og
# dev options
# DEV =
DEV = -g -Wall -Wextra -fsignaling-nans -Warray-bounds -Werror  -fsanitize=address
# Warning: -ffinite-math-only assumes all math will be finite, 
# and won't catch infs and nans

# Additional libraries
LIBS = 
# Additional header files
INCLUDES = -I./include
# linking flags
LDFLAGS  = -lm 

CFLAGS = $(DEFINES) $(OPT) $(OPTIONS) $(DEV) $(INCLUDES) $(LIBS)




#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(SRCDIR):./include

# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
EXECS = $(notdir $(patsubst %.c, %, $(wildcard $(SRCDIR)/*.c)))
# OBJECTS = main.o

HEADERS = cell.h histogram.h ioutils.h



#==================================
# RECEPIES
#==================================

all: default

default: $(EXECS)

# $(EXEC): $(OBJECTS) $(HEADERS)
#     $(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readTest1: readTest1.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readTest2: readTest2.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readTest3: readTest3.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readTest4: readTest4.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readTest4IC: readTest4IC.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readTest5: readTest5.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)


#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXECS)
