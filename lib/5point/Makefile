# 5point library makefile

LIB = lib5point.a

TARGET = $(LIB)
OBJS = 5point.o poly1.o poly3.o

CC=gcc
OPTFLAGS=-O3
OTHERFLAGS=-Wall

MATRIX_PATH=../matrix
IMAGELIB_PATH=../imagelib

CPPFLAGS = $(OTHERFLAGS) $(OPTFLAGS) -I$(MATRIX_PATH) -I$(IMAGELIB_PATH)

all: $(TARGET)

$(TARGET): $(OBJS)
	ar r $@ $(OBJS)
	cp $@ ..

clean: 
	rm -f $(TARGET) *.o *~
