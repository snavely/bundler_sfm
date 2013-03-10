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

#ifndef __pyramid_io_h__
#define __pyramid_io_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "pyramid.h"
#include "distpyr.h"

/* ------- Image pyramids ------- */

/* Render all images in the pyramid into a single image */
img_t *img_pyr_render(img_pyr_t *pyr);

/* Free a gaussian pyramid */
void img_pyr_free(img_pyr_t *pyr);


/* ------- Distance pyramids ------- */

/* Write a distance pyramid to a file */
void img_write_distance_pyramid(FILE *f, img_dist_pyr_t *pyr);

/* Write a distance pyramid to a file */
void img_write_distance_pyramid_file(img_dist_pyr_t *pyr, char *fname);

/* Read a distance pyramid from a file */
img_dist_pyr_t *img_read_distance_pyramid(FILE *f);

/* Read a distance pyramid from a file */
img_dist_pyr_t *img_read_distance_pyramid_file(char *fname);

/* Convert a distance map to a distance pyramid */
img_dist_pyr_t *dmap2dpyr(img_dmap_t *dmap);

/* Free a distance pyramid */
void img_free_distance_pyramid(img_dist_pyr_t *pyr);

/* Render the distance pyramid as an image */
img_t *img_distance_pyramid_render(img_dist_pyr_t *pyr, int mode);

#ifdef __cplusplus
}
#endif

#endif /* __pyramid_io_h__ */
