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

/* vanish.c */
/* Compute vanishing points from an image */

#ifndef __vanish_h__
#define __vanish_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "image.h"

/* Find vanishing points using a Hough transform */
void compute_vps_hough(img_t *img);

/* Find vanishing points using a Hough transform */
void compute_vps_hough_edges(img_t *img, int num_edges, edge_link_t *edges);

#ifdef __cplusplus
}
#endif

#endif /* __vanish_h__ */
