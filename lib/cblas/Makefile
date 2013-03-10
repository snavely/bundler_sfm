# This Makefile compiles the CBLAS routines
#
CC=gcc
CFLAGS = -O3 -DADD_
FFLAGS = -O3

all: libcblas.a

errhand = cblas_globals.o cblas_xerbla.o xerbla.o

#
# All object files for double real precision
#
dlev3 = cblas_dgemm.o

#
# Make an archive file
#

.SUFFIXES: .o .c .f

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f

libcblas.a: $(dlev3) $(errhand)
	$(AR) rs $@ $(dlev3) $(errhand)
	cp $@ ../

#
# Erase all object and archive files
#
clean:
	rm -f *.o *~

