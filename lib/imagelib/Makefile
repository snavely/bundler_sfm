# Makefile for imagelib

GCC			= g++

CC=gcc
OPTFLAGS=-O3
OTHERFLAGS=-Wall

IMAGELIB_OBJS= affine.o bmp.o canny.o color.o fileio.o filter.o fit.o	\
	fmatrix.o homography.o horn.o image.o lerp.o morphology.o	\
	pgm.o poly.o qsort.o resample.o tps.o transform.o		\
	triangulate.o util.o

INCLUDE_PATH=-I../matrix

CFLAGS=$(OTHERFLAGS) $(OPTFLAGS) $(INCLUDE_PATH)

IMAGELIB=libimage.a

all: $(IMAGELIB)

$(IMAGELIB): $(IMAGELIB_OBJS)
	ar r $(IMAGELIB) $(IMAGELIB_OBJS)
	cp $(IMAGELIB) ..

clean:
	rm -f *.o *~ $(IMAGELIB)
