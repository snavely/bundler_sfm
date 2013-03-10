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

/* dmap-io.h */
/* Routines for reading / writing / rendering distance maps */

#ifndef __dmap_io_h__
#define __dmap_io_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "dmap.h"
#include "image.h"
#include "vector.h"

/* Create a clean dmap */
img_dmap_t *img_dmap_new(int w, int h);

/* Produces an image of a distance map */
img_t *img_dmap_render(img_dmap_t *dmap);

/* Produces an image of the disparity between two images */
img_t *img_dmap_render_disparity(img_dmap_t *dmap);

/* Produces an image of the flow map */
img_t *img_dmap_render_flow(img_dmap_t *dmap);

/* Write a distance map to a file */
void img_dmap_write(FILE *f, img_dmap_t *dmap);

/* Read a distance map from a file */
img_dmap_t *img_dmap_read(FILE *f);
img_dmap_t *img_dmap_read_file(char *fname);

/* Free a distance map */
void img_dmap_free(img_dmap_t *dmap);

#ifdef __cplusplus
}
#endif

#endif /* __dmap_io_h__ */
