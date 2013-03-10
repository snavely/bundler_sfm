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

/* canny.c */
/* Code for canny edge detection */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "canny.h"
#include "filter.h"
#include "fit.h"
#include "image.h"
#include "lerp.h"
#include "matrix.h"
#include "morphology.h"
#include "util.h"

int img_find_next_point(img_t *img, img_t *marked, int x, int y, 
			int *x_next, int *y_next)
{
    int w = img->w, h = img->h;
    int dx, dy;
    
    for (dy = -1; dy <= 1; dy++) {
	for (dx = -1; dx <= 1; dx++) {
	    int on, is_marked;

	    if (dx == 0 && dy == 0)
		continue;

	    if (x+dx < 0 || x+dx >= w || y+dy < 0 || y+dy >= h) {
		continue;
	    } else {
		on = (img_get_pixel(img, x+dx, y+dy).r == 0xff);
		is_marked = (img_get_pixel(marked, x+dx, y+dy).r == 0xff);

		if (on && !is_marked) {
		    *x_next = x + dx;
		    *y_next = y + dy;
		    return 1;
		}
	    }
	}
    }

    return 0;
}

edge_link_t *img_track_edge(img_t *img, img_t *marked, int x, int y) {
    edge_link_t *head = (edge_link_t *) malloc(sizeof(edge_link_t));
    edge_link_t *tail = NULL;

    // printf("[img_track_edge] Tracking edge [%d,%d]\n", x, y);

    head->x = x;
    head->y = y;
    head->next = NULL;
    img_set_pixel(marked, x, y, 0xff, 0xff, 0xff);

    tail = head;

    /* Track the edge forward */
    // printf("[img_track_edge] Forward\n");

    while (1) {
	edge_link_t *next;
	int x_next, y_next;
	int found = img_find_next_point(img, marked, x, y, &x_next, &y_next);

	if (!found)
	    break;
	
	// printf("  x_next, y_next = %d, %d\n", x_next, y_next);

	next = (edge_link_t *) malloc(sizeof(edge_link_t));
	next->x = x_next;
	next->y = y_next;
	next->next = NULL;
	img_set_pixel(marked, x_next, y_next, 0xff, 0xff, 0xff);

	tail->next = next;
	tail = next;

	x = x_next;
	y = y_next;
    }

    /* Track the edge backward */
    // printf("[img_track_edge] Backward\n");

    x = iround(head->x);
    y = iround(head->y);

    while (1) {
	edge_link_t *prev;
	int x_prev, y_prev;
	int found = img_find_next_point(img, marked, x, y, &x_prev, &y_prev);

	if (!found)
	    break;

	// printf("  x_prev, y_prev = %d, %d\n", x_prev, y_prev);
	
	prev = (edge_link_t *) malloc(sizeof(edge_link_t));
	prev->x = x_prev;
	prev->y = y_prev;
	img_set_pixel(marked, x_prev, y_prev, 0xff, 0xff, 0xff);

	prev->next = head;
	head = prev;

	x = x_prev;
	y = y_prev;
    }

    return head;
}

/* Link foreground pixels in binary image `img' */
void img_link_edges(img_t *img, int *num_links_out, edge_link_t **links) 
{
    int w = img->w;
    int h = img->h;
    int x, y;
    int done;

    img_t *marked = img_new(w, h);
    int num_links = 0;
    
    // printf("[img_link_edges] Linking edges\n");

    do {
	done = 1;

	for (y = 0; y < h; y++) {
	    for (x = 0; x < w; x++) {
		int on = (img_get_pixel(img, x, y).r == 0xff);
		int is_marked = (img_get_pixel(marked, x, y).r == 0xff);

		if (on && !is_marked) {
		    /* Track this edge */
		    edge_link_t *edge = img_track_edge(img, marked, x, y);

		    num_links++;

		    if (num_links == 1) {
			*links = (edge_link_t *) malloc(sizeof(edge_link_t));
		    } else {
			*links = 
			    (edge_link_t *) realloc(*links,
						    num_links * 
						    sizeof(edge_link_t));
		    }

		    (*links)[num_links-1].next = edge;

		    done = 0;
		}
	    }
	}
    } while (!done);

    *num_links_out = num_links;

    img_free(marked);
}

/* Fit a line to the given edge points */
void edge_fit_line(edge_link_t *edge, double scale, 
		   double xbias, double ybias, double *params) 
{
    int count, num_pts;
    v2_t *pts;
    edge_link_t *link;
    
    num_pts = 0;
    for (link = edge->next; link != NULL; link = link->next) {
	num_pts++;
    }

    pts = (v2_t *) malloc(sizeof(v2_t) * num_pts);
    
    count = 0;
    for (link = edge->next; link != NULL; link = link->next) {
	pts[count] = 
	    v2_new((link->x - xbias) * scale, (link->y - ybias) * scale);

	count++;
    }

    fit_2D_line(num_pts, pts, params);

    free(pts);
}

/* Remove small edge links */
void edge_remove_small(int num_links_in, edge_link_t **links_in,
		       int *num_links_out, edge_link_t **links_out,
		       int min_length)
{
    int i;
    *links_out = (edge_link_t *) malloc(sizeof(edge_link_t) * num_links_in);

    *num_links_out = 0;

    for (i = 0; i < num_links_in; i++) {
	/* Count length of chain */
	edge_link_t *link = (*links_in)[i].next;
	int length = 0;

	while (link != NULL) {
	    length++;
	    link = link->next;
	}

	if (length > min_length) {
	    (*links_out)[*num_links_out] = (*links_in)[i];
	    (*num_links_out)++;
	}
    }
}

/* Remove large edge links */
void edge_remove_large(int num_links_in, edge_link_t **links_in,
		       int *num_links_out, edge_link_t **links_out,
		       int max_length)
{
    int i;
    *links_out = (edge_link_t *) malloc(sizeof(edge_link_t) * num_links_in);

    *num_links_out = 0;

    for (i = 0; i < num_links_in; i++) {
	/* Count length of chain */
	edge_link_t *link = (*links_in)[i].next;
	int length = 0;

	while (link != NULL) {
	    length++;
	    link = link->next;
	}

	if (length < max_length) {
	    (*links_out)[*num_links_out] = (*links_in)[i];
	    (*num_links_out)++;
	}
    }
}

void edge_free(edge_link_t *head) 
{
    edge_link_t *link = head->next;
    
    while (link != NULL) {
	edge_link_t *next = link->next;
	free(link);
	link = next;
    }
}

/* Divide edge links at places where the curvature is high */
void img_break_edges(int num_links_in,  edge_link_t **links_in,
		     int *num_links_out, edge_link_t **links_out, 
		     double dev_thresh) 
{
    int i;
    double dev_thresh_sq = dev_thresh * dev_thresh;
    int changed = 0;

    /* Copy the original links */
    *links_out = (edge_link_t *) malloc(sizeof(edge_link_t) * num_links_in);
    for (i = 0; i < num_links_in; i++) {
	(*links_out)[i] = (*links_in)[i];
    }
    *num_links_out = num_links_in;

    for (i = 0; i < num_links_in; i++) {
	/* Fit a line to the points in the edge */
	int num_pts = 0;
	v2_t pts[2];
	edge_link_t *link = (*links_in)[i].next;
	int /* j, */ count;
	double params[2];
	double max_dev;
	int max_idx;

	while (link != NULL) {
	    num_pts++;
	    link = link->next;
	}

	if (num_pts < 2)
	    continue;

	// pts = (v2_t *) malloc(sizeof(v2_t) * num_pts);

#if 0
	count = 0;
	link = &((*links_in)[i]);
	while (link != NULL) {
	    pts[count] = v2_new(link->x, link->y);
	    link = link->next;
	    count++;
	}

	fit_2D_line(num_pts, pts, params);
#endif

	link = (*links_in)[i].next;
	pts[0] = v2_new(link->x, link->y);

	while (link->next != NULL) {
	    link = link->next;
	}
	
	pts[1] = v2_new(link->x, link->y);

	fit_2D_line(2, pts, params);
	// fit_2D_line_orthogonal_regression(2, pts, params);

	/* Find the max deviation */
	max_dev = 0.0;
	max_idx = -1;

	count = 0;
	link = (*links_in)[i].next;
	while (link != NULL) {
	    double dev = link->x * params[0] + link->y * params[1] + 1.0;
	    dev = dev * dev / (params[0] * params[0] + params[1] * params[1]);

	    if (dev > max_dev) {
		max_dev = dev;
		max_idx = count;
	    }

	    link = link->next;
	    count++;
	}

	if (max_dev > dev_thresh_sq) {
	    edge_link_t *head_new = NULL;

	    /* Split the edge */
	    count = 0;

	    link = (*links_in)[i].next;
	    while (count < max_idx) {
		link = link->next;
		count++;
	    }

	    head_new = link->next;
	    link->next = NULL;
	    
	    if (head_new != NULL) {
		(*num_links_out)++;
		*links_out = 
		    (edge_link_t *) realloc(*links_out,
					    sizeof(edge_link_t) * 
					    (*num_links_out));
		(*links_out)[(*num_links_out)-1].next = head_new;

		changed = 1;
	    }
	}
    }

    if (changed) {
	edge_link_t *links_out_new;

	/* Repeat */
	num_links_in = *num_links_out;
	links_in = links_out;

	img_break_edges(num_links_in, links_in,
			num_links_out, &links_out_new, dev_thresh);

	free(*links_out);
	*links_out = links_out_new;
    }
}

/* Detect edges using Canny's algorithm */
img_t *img_canny_edge_detect(img_t *img, double sigma, double t1, double t2)
{
    int w = img->w, h = img->h;
    int x, y;

    img_t *edges, *edges_tmp, *edges_tmp2, *thin;

    float *gmag, *gtheta;

    img_t *gray = img_convert_grayscale(img);
    fimg_t *smoothed = img_smooth_float(gray, sigma, 0);
    img_free(gray);

#if 0
    img_t *smoothed_int = fimg2img(smoothed);
    img_write_bmp_file(smoothed_int, "smooth.bmp");
    img_free(smoothed_int);
#endif

    /* Compute gradient magnitudes */
    gmag = (float *) malloc(sizeof(float) * w * h);
    gtheta = (float *) malloc(sizeof(float) * w * h);

    for (y = 1; y < h - 1; y++) {
	for (x = 1; x < w - 1; x++) {
	    double dx = 0.5 * (fimg_get_pixel(smoothed, x+1, y) -
			       fimg_get_pixel(smoothed, x-1, y)) / 256;
	    double dy = 0.5 * (fimg_get_pixel(smoothed, x, y+1) -
			       fimg_get_pixel(smoothed, x, y-1)) / 256;

	    double mag = sqrt(dx * dx + dy * dy);
	    double theta = atan2(dy, dx);

	    gmag[y * w + x] = (float) mag;
	    gtheta[y * w + x] = (float) theta;
	}
    }

    fimg_free(smoothed);

    /* Find the maxima */
    edges_tmp = img_new(w, h);

#define RADIUS 0.8
    for (y = 10; y < h - 10; y++) {
	for (x = 10; x < w - 10; x++) {
	    double mag = gmag[y * w + x];
	    double theta = gtheta[y * w + x];

	    double x_left = x + RADIUS * cos(theta);
	    double y_left = y + RADIUS * sin(theta);

	    double x_rt = x - RADIUS * cos(theta);
	    double y_rt = y - RADIUS * sin(theta);

	    double g_left = func_lerpf(w, h, gmag, x_left, y_left);
	    double g_rt = func_lerpf(w, h, gmag, x_rt, y_rt);

	    if (mag > g_left && mag > g_rt) {
		// printf("[%d,%d] mag, [%0.3f,%0.3f] l, [%0.3f,%0.3f] "
		//        "r = %0.3f, %0.3f, %0.3f\n", 
		//        x, y, x_left, y_left, x_rt, y_rt, mag, g_left, g_rt);
		img_set_pixel(edges_tmp, x, y, 0xff, 0xff, 0xff);
	    } else {
		img_set_pixel(edges_tmp, x, y, 0, 0, 0);
	    }
	}
    }

    /* Now do the thresholding */
    edges_tmp2 = img_new(w, h);

    /* First pass */
    for (y = 10; y < h - 10; y++) {
	for (x = 10; x < w - 10; x++) {
	    if (img_get_pixel(edges_tmp, x, y).r == 0x0) 
		continue;
	    
	    if (gmag[y * w + x] > t1)
		img_set_pixel(edges_tmp2, x, y, 0xff, 0xff, 0xff);
	}
    }



    edges = img_new(w, h);

    /* Second pass */
    for (y = 10; y < h - 10; y++) {
	for (x = 10; x < w - 10; x++) {
	    if (img_get_pixel(edges_tmp2, x, y).r == 0xff) {
		int dx, dy;

		img_set_pixel(edges, x, y, 0xff, 0xff, 0xff);

		for (dy = -1; dy <= 1; dy++) {
		    for (dx = -1; dx <= 1; dx++) {
			int is_edge = 
			    img_get_pixel(edges_tmp, x+dx, y+dx).r == 1;

			if (is_edge && gmag[(y+dy)*w + (x+dx)] > t2) {
			    img_set_pixel(edges, x, y, 0xff, 0xff, 0xff);
			}
		    }
		}
	    }
	}
    }

    img_free(edges_tmp);
    img_free(edges_tmp2);

    free(gmag);
    free(gtheta);

    thin = img_thin(edges);    

    img_free(edges);

    return thin;
    // return edges;
}
