OBJS=dpmpar.o \
	enorm.o \
	fdjac2.o \
	lmdif1.o \
	lmdif.o \
	lmpar.o \
	qrfac.o \
	qrsolv.o

AR=ar
FC=gcc
FFLAGS=-O2 -ffast-math

all: libminpack.a

libminpack.a: $(OBJS)
	$(AR) rs $@ $(OBJS)
	cp $@ ../

clean:
	rm -f *.o *~ libminpack.a
