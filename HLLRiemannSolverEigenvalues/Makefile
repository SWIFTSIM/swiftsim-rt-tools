#==================================
# COMPILATION TIME PARAMETERS
#==================================

#directory where source files are
SRCDIR=.


#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
F90 = gfortran

# flags for gfortran compiler
DEBUGFLAGS = -fbacktrace -g -O 
CHECKS = -fbounds-check -Wuninitialized -Wall -fcheck=all -ffpe-trap=zero,underflow,overflow,invalid -Werror
EXTRAS = -ffree-line-length-none -finit-real=nan -cpp
INCFLAGS = 
# INCFLAGS = -I/${OPENBLAS_ROOT}/include
LDFLAGS = 
# LDFLAGS = -L/${OPENBLAS_ROOT}/lib
LIBFLAGS = -llapack -lblas
F90FLAGS= $(DEBUGFLAGS) $(CHECKS) $(EXTRAS)



default: get_eigenvalues.o

get_eigenvalues.o: get_eigenvalues.f90
	$(F90) $(F90FLAGS) get_eigenvalues.f90 $(INCFLAGS) $(LDFLAGS) $(LIBFLAGS) -o get_eigenvalues.o


clean:
	rm -f get_eigenvalues.o
