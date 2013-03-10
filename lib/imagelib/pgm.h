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

/* pgm.h */
/* Read PGM files */

#ifndef __pgm_h__
#define __pgm_h__

#ifdef __cplusplus
extern "C" {
#endif

/* Read a PGM file */
img_t *img_read_pgm_file(const char *filename);

/* Read the dimensions of an image */
int img_read_pgm_dimensions(char *filename, int *w, int *h);

#ifdef __cplusplus
}
#endif

#endif /* __pgm_h__ */
