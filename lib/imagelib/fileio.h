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

#ifndef __fileio_h__
#define __fileio_h__

#include <sys/types.h>

/* Open a file, exiting on error */
FILE *open_file(char *fname, char *flags);

/* Read/write various sized items */
void read_word(u_int32_t *w, FILE *f);
void read_short(u_int16_t *s, FILE *f);
void read_byte(u_int8_t *b, FILE *f);
void read_float(float *fl, FILE *f);
void read_double(double *d, FILE *f);
void write_word(u_int32_t *w, FILE *f);
void write_short(u_int16_t *s, FILE *f);
void write_byte(u_int8_t *b, FILE *f);
void write_float(float *fl, FILE *f);
void write_double(double *d, FILE *f);

#endif /* __fileio_h__ */
