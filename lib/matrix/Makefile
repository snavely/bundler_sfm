# Makefile for matrix library

CC=gcc
OPTFLAGS=-O3
OTHERFLAGS=-Wall -D__NO_MKL__

INCLUDE_PATH = -I../imagelib -I../../include
CFLAGS = $(OPTFLAGS) $(OTHERFLAGS) $(INCLUDE_PATH)

TARGET = libmatrix.a
OBJS = matrix.o vector.o svd.o

all: $(TARGET)

$(TARGET): $(OBJS)
	ar r $@ $(OBJS)
	cp $@ ..

clean:
	rm -f *.o $(TARGET) *~
