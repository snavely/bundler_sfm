/* 
 *  Copyright (c) 2008  Noah Snavely (snavely (at) cs.washington.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* fileio.c */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#ifdef WIN32
#include "types.h"
#endif

FILE *open_file(char *fname, char *flags) {
    FILE *f = fopen(fname, flags);
    if (f == NULL) {
	printf("Error opening file %s\n", fname);
	exit(1);
    }
    
    return f;
}

void read_word(u_int32_t *w, FILE *f) {
    fread(w, 4, 1, f);
}

void read_short(u_int16_t *s, FILE *f) {
    fread(s, 2, 1, f);
}

void read_byte(u_int8_t *b, FILE *f) {
    fread(b, 1, 1, f);
}

void read_float(float *fl, FILE *f) {   
    fread(fl, sizeof(float), 1, f);
}

void read_double(u_int8_t *d, FILE *f) {
    fread(d, sizeof(double), 1, f);
}

void write_word(u_int32_t *w, FILE *f) {
    fwrite(w, 4, 1, f);
}

void write_short(u_int16_t *s, FILE *f) {
    fwrite(s, 2, 1, f);
}

void write_byte(u_int8_t *b, FILE *f) {
    fwrite(b, 1, 1, f);
}

void write_float(float *fl, FILE *f) {
    fwrite(fl, sizeof(float), 1, f);
}
    

void write_double(u_int8_t *d, FILE *f) {
    fwrite(d, sizeof(double), 1, f);
}
