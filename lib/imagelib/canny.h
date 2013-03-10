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

/* canny.h */
/* Code for canny edge detection */

#ifndef __canny_h__
#define __canny_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "image.h"

/* Edge link structure */
typedef struct edge_link_t {
    struct edge_link_t *next;
    double x, y;
    int idx;
} edge_link_t;

/* Link foreground pixels in binary image `img' */
void img_link_edges(img_t *img, int *num_links_out, edge_link_t **links);

/* Detect edges using Canny's algorithm */
img_t *img_canny_edge_detect(img_t *img, double sigma, double t1, double t2);

/* Divide edge links at places where the curvature is high */
void img_break_edges(int num_links_in,  edge_link_t **links_in,
		     int *num_links_out, edge_link_t **links_out, 
		     double dev_thresh);

/* Fit a line to the given edge points */
void edge_fit_line(edge_link_t *edge, double scale, 
		   double xbias, double ybias, double *params);

/* Remove small edge links */
void edge_remove_small(int num_links_in, edge_link_t **links_in,
		       int *num_links_out, edge_link_t **links_out,
		       int min_length);

/* Remove large edge links */
void edge_remove_large(int num_links_in, edge_link_t **links_in,
		       int *num_links_out, edge_link_t **links_out,
		       int max_length);

/* Free the given edge link */
void edge_free(edge_link_t *head);

#ifdef __cplusplus
}
#endif

#endif /* __canny_h__ */
