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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "canny.h"
#include "defines.h"
#include "filter.h"
#include "image.h"
#include "matrix.h"
#include "util.h"
#include "vector.h"

/* Find vanishing points using a Hough transform */
void compute_vps_hough(img_t *img) 
{
    int w = img->w, h = img->h;
    double w2 = 0.5 * img->w, h2 = 0.5 * img->h;

    double size = MAX(w, h);

    int x, y;
    int t, p;

    /* Smooth the image and compute gradients */
    img_t *gray = img_convert_grayscale(img);
    img_t *smoothed = img_smooth(gray, 0.5, 0);

    double *gmag = NULL;
    double *gtheta = NULL;
    
    double *hough = NULL;
    double *lines = NULL;
    double *line_weights = NULL;

    double sigma_sq;

    int num_lines = 0;

    int i, count;

    img_t *img_out;
    double hmin, hmax;
    int tmax, pmax;
    double theta_max, phi_max, xmax, ymax;

    img_write_bmp_file(smoothed, "smooth.bmp");

    gmag = (double *) malloc(sizeof(double) * w * h);
    gtheta = (double *) malloc(sizeof(double) * w * h);

#define GRADIENT_THRESHOLD 0.04
    for (y = 1; y < h - 1; y++) {
	for (x = 1; x < w - 1; x++) {
	    double dx = 0.5 * (img_get_pixel(smoothed, x+1, y).r -
			       img_get_pixel(smoothed, x-1, y).r) / 256;
	    double dy = 0.5 * (img_get_pixel(smoothed, x, y+1).r -
			       img_get_pixel(smoothed, x, y-1).r) / 256;

	    double mag = sqrt(dx * dx + dy * dy);
	    double theta = atan2(dy, dx);

	    if (mag > GRADIENT_THRESHOLD)
		printf("%0.3f, %0.3f\n", dx, dy);
	    
	    gmag[y * w + x] = mag;
	    gtheta[y * w + x] = theta;
	}
    }

    /* Threshold the gradients */
    for (y = 1; y < h - 1; y++) {
	for (x = 1; x < w - 1; x++) {
	    // printf("g[%d,%d] = %0.3f\n", x, y, gmag[y * w + x]);
	    
	    if (gmag[y * w + x] >= GRADIENT_THRESHOLD) 
		num_lines++;
	}
    }
    
    printf("[compute_vps_hough] Found %d edgelets\n", num_lines);

    /* Compute the line parameters at each image point */
    lines = (double *) malloc(sizeof(double) * 3 * num_lines);
    line_weights = (double *) malloc(sizeof(double) * num_lines);

    count = 0;
    for (y = 1; y < h - 1; y++) {
	for (x = 1; x < w - 1; x++) {
	    int idx = y * w + x;

	    if (gmag[idx] < GRADIENT_THRESHOLD) {
		continue;
	    } else {
		double x1 = (x - w2) / size, y1 = (y - h2) / size;
		double x2 = x1 - sin(gtheta[idx]), y2 = y1 + cos(gtheta[idx]);
	    
		double p[3] = { x1, y1, 1.0 };
		double q[3] = { x2, y2, 1.0 };

		double norm;

		// img_set_pixel(smoothed, x, y, 0xff, 0, 0);
		if (count % 30 == 0) {
		    img_draw_line(smoothed, x, y, 
				  x - iround(10 * sin(gtheta[idx])),
				  y + iround(10 * cos(gtheta[idx])), 0xff, 0, 0);
		}
		
		printf("theta[%d] = %0.3f\n", count, gtheta[idx]);

		matrix_cross(p, q, lines + 3 * count);
	    
		norm = matrix_norm(3, 1, lines + 3 * count);
		matrix_scale(3, 1, lines + 3 * count, 
			     1.0 / norm, lines + 3 * count);

		line_weights[count] = gmag[idx];

		printf("line[%d] = %0.3f, %0.3f, %0.3f\n", count,
		       lines[3 * count + 0], 
		       lines[3 * count + 1], 
		       lines[3 * count + 2]);

		count++;
	    }
	}
    }

    /* Now do a hough transform to find vanishing points */
#define HOUGH_SAMPLES_THETA 512
#define HOUGH_SAMPLES_PHI 512
#define HOUGH_SIGMA 0.015 // 0.01 // 0.005

    sigma_sq = HOUGH_SIGMA * HOUGH_SIGMA;

    hough = (double *) malloc(sizeof(double) * 
			      HOUGH_SAMPLES_THETA * HOUGH_SAMPLES_PHI);

    printf("[compute_vps_hough] Starting transform...\n");

    for (t = 0; t < HOUGH_SAMPLES_THETA; t++) {
	printf(".");
	fflush(stdout);

	for (p = 0; p < HOUGH_SAMPLES_PHI; p++) {
	    int hidx = t * HOUGH_SAMPLES_PHI + p;

	    double theta = M_PI * t / HOUGH_SAMPLES_THETA - 0.5 * M_PI;
	    double phi = M_PI * p / HOUGH_SAMPLES_PHI; // - 0.5 * M_PI;

	    /* Find homogeneous image coords of this cell */
	    double x_c = cos(theta) * sin(phi);
	    double y_c = sin(theta) * sin(phi);
	    double z_c = cos(phi);

	    double c[3] = { x_c, y_c, z_c };

	    hough[hidx] = 0.0;

	    for (i = 0; i < num_lines; i++) {
		double dot;
		double weight;

		matrix_product(1, 3, 3, 1, c, lines + 3 * i, &dot);
		weight = exp(-dot * dot / sigma_sq);

		/* Accumulate */
		hough[hidx] += weight * line_weights[i];
	    }
	}
    }
    printf("\n");

    /* Write an output image */
    hmin = DBL_MAX;
    hmax = 0.0;
    tmax = pmax = 0;

    for (t = 0; t < HOUGH_SAMPLES_THETA; t++) {
	for (p = 0; p < HOUGH_SAMPLES_PHI; p++) {
	    int hidx = t * HOUGH_SAMPLES_PHI + p;
	    if (hough[hidx] < hmin)
		hmin = hough[hidx];
	    if (hough[hidx] > hmax) {
		tmax = t;
		pmax = p;
		hmax = hough[hidx];
	    }
	}
    }

    theta_max = M_PI * tmax / HOUGH_SAMPLES_THETA - 0.5 * M_PI;
    phi_max = M_PI * pmax / HOUGH_SAMPLES_PHI; // - 0.5 * M_PI;
    xmax = size * cos(theta_max) * sin(phi_max) / cos(phi_max) + w2;
    ymax = size * sin(theta_max) * sin(phi_max) / cos(phi_max) + h2;

    img_draw_pt(smoothed, iround(xmax), iround(ymax), 4, 0x0, 0x0, 0xff);

    printf("[compute_vps_hough] hmin, hmax = %0.3f, %0.3f\n", hmin, hmax);
    printf("[compute_vps_hough] xmax, ymax = %0.3f, %0.3f\n", xmax, ymax);

    img_out = img_new(HOUGH_SAMPLES_PHI, HOUGH_SAMPLES_THETA);

    for (t = 0; t < HOUGH_SAMPLES_THETA; t++) {
	for (p = 0; p < HOUGH_SAMPLES_PHI; p++) {
	    int hidx = t * HOUGH_SAMPLES_PHI + p;
	    double hnorm = (hough[hidx] - hmin) / (hmax - hmin);
	    int c = iround(255.0 * hnorm);

	    img_set_pixel(img_out, t, p, c, c, c);
	}
    }
    
    img_write_bmp_file(img_out, "out.bmp");
    img_write_bmp_file(smoothed, "marked.bmp");

    img_free(gray);
    img_free(smoothed);

    free(gmag);
    free(gtheta);
    free(hough);
    free(lines);
    free(line_weights);
}

/* Find vanishing points using a Hough transform */
void compute_vps_hough_edges(img_t *img, int num_edges, edge_link_t *edges) 
{
    int w = img->w, h = img->h;
    double w2 = 0.5 * img->w, h2 = 0.5 * img->h;
    double size = MAX(w, h);
    double size2 = 0.5 * size;

    int i, p, t, round;
    double *line_params = NULL;
    int *line_used = NULL;
    double *hough = NULL;
    double sigma_sq;

    double hmin, hmax;
    int tmax, pmax;
    double theta_max, phi_max, xmax, ymax;

    img_t *img_out;

    line_params = (double *) malloc(sizeof(double) * 3 * num_edges);
    line_used = (int *) malloc(sizeof(int) * num_edges);

    for (i = 0; i < num_edges; i++) {
	double norm;

	edge_fit_line(&(edges[i]), 1.0 / size2, 0.5 * img->w, 0.5 * img->h, 
		      line_params + 3 * i);
	line_params[3 * i + 2] = 1.0;

	norm = matrix_norm(3, 1, line_params + 3 * i);
	matrix_scale(3, 1, 
		     line_params + 3 * i, 1.0 / norm, 
		     line_params + 3 * i);

	line_used[i] = 0;

	edges[i].idx = -1;
    }

    /* Now do a hough transform to find vanishing points */
    sigma_sq = HOUGH_SIGMA * HOUGH_SIGMA;

    hough = (double *) malloc(sizeof(double) * 
			      HOUGH_SAMPLES_THETA * HOUGH_SAMPLES_PHI);

    printf("[compute_vps_hough_edges] Starting transform...\n");

    for (round = 0; round < 3; round++) {
	double x_cmax, y_cmax, z_cmax;
	double c_max[3];
	char outfile[256];

	for (t = 0; t < HOUGH_SAMPLES_THETA; t++) {
	    printf(".");
	    fflush(stdout);

	    for (p = 0; p < HOUGH_SAMPLES_PHI; p++) {
		int hidx = t * HOUGH_SAMPLES_PHI + p;

		double theta = 2.0 * M_PI * t / HOUGH_SAMPLES_THETA;
		double phi = 0.5 * M_PI * p / (HOUGH_SAMPLES_PHI - 1);

		/* Find homogeneous image coords of this cell */
		double x_c = cos(theta) * sin(phi);
		double y_c = sin(theta) * sin(phi);
		double z_c = cos(phi);

		double c[3] = { x_c, y_c, z_c };

		hough[hidx] = 0.0;
		
		for (i = 0; i < num_edges; i++) {
		    double dot;
		    double weight;

		    if (line_used[i])
			continue;

#if 0
		    if (round == 0 && line_used[i])
			continue;
		
		    if (round == 1 && !line_used[i])
			continue;
#endif

		    matrix_product(1, 3, 3, 1, c, line_params + 3 * i, &dot);
		    weight = exp(-dot * dot / sigma_sq);

		    /* Accumulate */
		    hough[hidx] += weight /* * line_weights[i] */;
		}
	    }
	}
	printf("\n");

	/* Write an output image */
	hmin = DBL_MAX;
	hmax = 0.0;
	tmax = pmax = 0;

	for (t = 0; t < HOUGH_SAMPLES_THETA; t++) {
	    for (p = 0; p < HOUGH_SAMPLES_PHI; p++) {
		int hidx = t * HOUGH_SAMPLES_PHI + p;
		if (hough[hidx] < hmin)
		    hmin = hough[hidx];
		if (hough[hidx] > hmax) {
		    tmax = t;
		    pmax = p;
		    hmax = hough[hidx];
		}
	    }
	}

	printf("[compute_vps_hough_edge] tmax, pmax = %d, %d\n", tmax, pmax);
	
	theta_max = 2.0 * M_PI * tmax / HOUGH_SAMPLES_THETA; //  - 0.5 * M_PI;
	phi_max = 0.5 * M_PI * pmax / (HOUGH_SAMPLES_PHI - 1); // - 0.5 * M_PI;
	xmax = size2 * cos(theta_max) * sin(phi_max) / cos(phi_max) + w2;
	ymax = size2 * sin(theta_max) * sin(phi_max) / cos(phi_max) + h2;

	if (round == 0)
	    img_draw_pt(img, iround(xmax), iround(ymax), 
			8, 0x0, 0x0, 0xff);
	else if (round == 1)
	    img_draw_pt(img, iround(xmax), iround(ymax), 
			8, 0x0, 0xff, 0x0);
	else if (round == 2)
	    img_draw_pt(img, iround(xmax), iround(ymax), 
			8, 0xff, 0x0, 0x0);

	printf("[compute_vps_hough] hmin, hmax = %0.3f, %0.3f\n", hmin, hmax);
	printf("[compute_vps_hough] xmax, ymax = %0.3f, %0.3f\n", xmax, ymax);

	img_out = img_new(HOUGH_SAMPLES_PHI, HOUGH_SAMPLES_THETA);

	for (t = 0; t < HOUGH_SAMPLES_THETA; t++) {
	    for (p = 0; p < HOUGH_SAMPLES_PHI; p++) {
		int hidx = t * HOUGH_SAMPLES_PHI + p;
		double hnorm = (hough[hidx] - hmin) / (hmax - hmin);
		int c = iround(255.0 * hnorm);

		img_set_pixel(img_out, t, p, c, c, c);
	    }
	}

	sprintf(outfile, "out%03d.bmp", round);
	img_write_bmp_file(img_out, outfile);
	img_free(img_out);

	/* Find the lines associated with this vanishing point */
	x_cmax = cos(theta_max) * sin(phi_max);
	y_cmax = sin(theta_max) * sin(phi_max);
	z_cmax = cos(phi_max);

	c_max[0] = x_cmax;
	c_max[1] = y_cmax;
	c_max[2] = z_cmax;

	for (i = 0; i < num_edges; i++) {
	    double dot;
		
	    matrix_product(1, 3, 3, 1, c_max, line_params + 3 * i, &dot);

	    if (fabs(dot) < 4.0 * HOUGH_SIGMA) {
		line_used[i] = 1;
		edges[i].idx = round;
	    }
	}
    }
    
    img_write_bmp_file(img, "marked.bmp");
}
