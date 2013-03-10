# Makefile for sfm driver

CC=gcc
OPTFLAGS=-O3
OTHERFLAGS=-Wall
INCLUDE_PATH=-I../matrix -I../imagelib -I../sba-1.5

CFLAGS=$(OTHERFLAGS) $(OPTFLAGS) $(INCLUDE_PATH)

LIBSFM=libsfmdrv.a

all: $(LIBSFM)

$(LIBSFM): sfm.o
	ar r $@ sfm.o
	cp $@ ..

clean:
	rm -f $(LIBSFM) *.o *~
