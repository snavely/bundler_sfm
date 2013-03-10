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

/* resample.c */
/* Contains routines for resampling an image given some transformation */

#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <float.h>

#include "bmp.h"
#include "image.h"
#include "lerp.h"
#include "transform.h"
#include "util.h"
#include "vector.h"

#define BILINEAR_INTERPOLATION
// #define LINEAR_INTERPOLATION
// #define BICUBIC_INTERPOLATION

static color_t background_color = { 0x0, 0x0, 0x0, 0x0 };

/* Set the background color for unmapped pixels in resampled images */
void set_resample_background(int r, int g, int b) {
    background_color.r = r;
    background_color.g = g;
    background_color.b = b;
}


/* Use linear interpolation to compute the intensity of the point (x,y) 
 * in the given image */
double pixel_lerp_intensity(img_t *img, double x, double y) {
    int xf = (int) floor(x), yf = (int) floor(y);
    double xp = x - xf, yp = y - yf;

#ifdef LINEAR_INTERPOLATION
    double f[4];
    double I;

    color_t *pixels[4];
    fcolor_t col;

    img_get_pixel_square(img, xf, yf, pixels);

    if (pixels[0] != NULL) {
	f[0] = color_intensity(*(pixels[0]));
    }

    if (pixels[1] != NULL) {
	f[1] = color_intensity(*(pixels[1]));
    }

    if (pixels[2] != NULL) {
	f[2] = color_intensity(*(pixels[2]));
    }
    
    if (pixels[3] != NULL) {
	f[3] = color_intensity(*(pixels[3]));
    }

    /* Lerp */
    if (xp >= yp) {
	I = f[0] + (f[1] - f[0]) * xp + (f[3] - f[1]) * yp;
    } else {
	I = f[0] + (f[2] - f[0]) * yp + (f[3] - f[2]) * xp;
    }
#endif

#ifdef BILINEAR_INTERPOLATION
    double params[2] = { xp, yp };

	color_t p1 = img_get_pixel(img, xf, yf);
	color_t p2 = img_get_pixel(img, xf + 1, yf);
    color_t p3 = img_get_pixel(img, xf, yf + 1);
	color_t p4 = img_get_pixel(img, xf + 1, yf + 1);

	color_t pixels[4];
	double f[4];
	double I;

	pixels[0] = p1;
	pixels[1] = p2;
	pixels[2] = p3;
	pixels[3] = p4;

    f[0] = color_intensity(pixels[0]);
	f[1] = color_intensity(pixels[1]);
	f[2] = color_intensity(pixels[2]);
	f[3] = color_intensity(pixels[3]);

    I = LERP(xp, yp, f[0], f[1], f[2], f[3]);
#endif    

#ifdef BICUBIC_INTERPOLATION
    color_t pixels[16] = { img_get_pixel(img, xf - 1, yf - 1),
			   img_get_pixel(img, xf + 0, yf - 1),
			   img_get_pixel(img, xf + 1, yf - 1),
			   img_get_pixel(img, xf + 2, yf - 1),

			   img_get_pixel(img, xf - 1, yf + 0),
			   img_get_pixel(img, xf + 0, yf + 0),
			   img_get_pixel(img, xf + 1, yf + 0),
			   img_get_pixel(img, xf + 2, yf + 0),		   

			   img_get_pixel(img, xf - 1, yf + 1),
			   img_get_pixel(img, xf + 0, yf + 1),
			   img_get_pixel(img, xf + 1, yf + 1),
			   img_get_pixel(img, xf + 2, yf + 1),		   

			   img_get_pixel(img, xf - 1, yf + 2),
			   img_get_pixel(img, xf + 0, yf + 2),
			   img_get_pixel(img, xf + 1, yf + 2),
			   img_get_pixel(img, xf + 2, yf + 2) };
    
    double f[16];

    double I;

    int i;

    for (i = 0; i < 16; i++) {
	f[i] = color_intensity(pixels[i]);
    }

    I = bicubic_interpolate_2D(xp, yp, f);
#endif

    return I;
}

/* Use linear interpolation to compute the color of the point (x,y) 
 * in the given image */
fcolor_t pixel_lerp(img_t *img, double x, double y) {
    int xf = (int) floor(x), yf = (int) floor(y);
    double xp = x - xf, yp = y - yf;

#ifdef LINEAR_INTERPOLATION
    int fr[4], fg[4], fb[4];
    double rd, gd, bd;

    color_t *pixels[4];
    fcolor_t col;

    img_get_pixel_square(img, xf, yf, pixels);

    if (pixels[0] != NULL) {
	fr[0] = pixels[0]->r;
	fg[0] = pixels[0]->g;
	fb[0] = pixels[0]->b;
    }

    if (pixels[1] != NULL) {
	fr[1] = pixels[1]->r;
	fg[1] = pixels[1]->g;
	fb[1] = pixels[1]->b;
    }

    if (pixels[2] != NULL) {
	fr[2] = pixels[2]->r;
	fg[2] = pixels[2]->g;
	fb[2] = pixels[2]->b;
    }
    
    if (pixels[3] != NULL) {
	fr[3] = pixels[3]->r;
	fg[3] = pixels[3]->g;
	fb[3] = pixels[3]->b;
    }

    /* Lerp */
    if (xp >= yp) {
	rd = ((double) fr[0]) + ((double) (fr[1] - fr[0])) * xp + ((double) (fr[3] - fr[1])) * yp;
	gd = ((double) fg[0]) + ((double) (fg[1] - fg[0])) * xp + ((double) (fg[3] - fg[1])) * yp;
	bd = ((double) fb[0]) + ((double) (fb[1] - fb[0])) * xp + ((double) (fb[3] - fb[1])) * yp;
    } else {
	rd = ((double) fr[0]) + ((double) (fr[2] - fr[0])) * yp + ((double) (fr[3] - fr[2])) * xp;
	gd = ((double) fg[0]) + ((double) (fg[2] - fg[0])) * yp + ((double) (fg[3] - fg[2])) * xp;
	bd = ((double) fb[0]) + ((double) (fb[2] - fb[0])) * yp + ((double) (fb[3] - fb[2])) * xp;
    }

    col.r = (float) rd;
    col.g = (float) gd;
    col.b = (float) bd;
#endif

#ifdef BILINEAR_INTERPOLATION
    double params[2] = { xp, yp };
	color_t pixels[4];
    fcolor_t col;
	double rd, gd, bd;

    pixels[0] = img_get_pixel(img, xf, yf);
	pixels[1] = img_get_pixel(img, xf + 1, yf);
	pixels[2] = img_get_pixel(img, xf, yf + 1);
	pixels[3] = img_get_pixel(img, xf + 1, yf + 1);

    //double fr[4] = { pixels[0].r, pixels[1].r, pixels[2].r, pixels[3].r };
    //double fg[4] = { pixels[0].g, pixels[1].g, pixels[2].g, pixels[3].g };
    //double fb[4] = { pixels[0].b, pixels[1].b, pixels[2].b, pixels[3].b };

    rd = LERP(xp, yp, pixels[0].r, pixels[1].r, pixels[2].r, pixels[3].r); // nlerp2(params, fr);
    gd = LERP(xp, yp, pixels[0].g, pixels[1].g, pixels[2].g, pixels[3].g); // nlerp2(params, fr);
    bd = LERP(xp, yp, pixels[0].b, pixels[1].b, pixels[2].b, pixels[3].b); // nlerp2(params, fr);
    // u_int8_t r = rint(rd);

    //double gd = nlerp2(params, fg);
    // u_int8_t g = rint(gd);

    //double bd = nlerp2(params, fb);
    // u_int8_t b = rint(bd);

    col.r = (float) rd;
    col.g = (float) gd;
    col.b = (float) bd;
#endif    

#ifdef BICUBIC_INTERPOLATION
    color_t pixels[16] = { img_get_pixel(img, xf - 1, yf - 1),
			   img_get_pixel(img, xf + 0, yf - 1),
			   img_get_pixel(img, xf + 1, yf - 1),
			   img_get_pixel(img, xf + 2, yf - 1),

			   img_get_pixel(img, xf - 1, yf + 0),
			   img_get_pixel(img, xf + 0, yf + 0),
			   img_get_pixel(img, xf + 1, yf + 0),
			   img_get_pixel(img, xf + 2, yf + 0),		   

			   img_get_pixel(img, xf - 1, yf + 1),
			   img_get_pixel(img, xf + 0, yf + 1),
			   img_get_pixel(img, xf + 1, yf + 1),
			   img_get_pixel(img, xf + 2, yf + 1),		   

			   img_get_pixel(img, xf - 1, yf + 2),
			   img_get_pixel(img, xf + 0, yf + 2),
			   img_get_pixel(img, xf + 1, yf + 2),
			   img_get_pixel(img, xf + 2, yf + 2) };
    
    double fr[16];
    double fg[16];
    double fb[16];

    double rd, gd, bd;
    u_int8_t r, g, b;

    fcolor_t col;

    int i;

    for (i = 0; i < 16; i++) {
	fr[i] = pixels[i].r;
	fg[i] = pixels[i].g;
	fb[i] = pixels[i].b;
    }

    rd = bicubic_interpolate_2D(xp, yp, fr);
    // rd = bicuberp(xp, yp, fr);
    // r = rint(rd);

    gd = bicubic_interpolate_2D(xp, yp, fg);
    // gd = bicuberp(xp, yp, fg);
    // g = rint(gd);

    bd = bicubic_interpolate_2D(xp, yp, fb);
    // bd = bicuberp(xp, yp, fb);
    // b = rint(bd);

    col.r = (float) rd;
    col.g = (float) gd;
    col.b = (float) bd;
#endif

    return col;
}

/* Compute the pixel that would be at location (x, y) if the given
 * image were transformed with the inverse of the given transform */
fcolor_t pixel_transform(img_t *img, trans2D_t *Tinv, int x, int y) {
    double Tp[2];
    
    int w = img->w, h = img->h;

    transform_point(Tinv, x, y, Tp+0, Tp+1);

    /* Check if the result is in range */
    if (Tp[0] < Vx(img->origin) || Tp[1] < Vy(img->origin) || 
	Tp[0] > Vx(img->origin) + w - 1 || Tp[1] > Vy(img->origin) + h - 1) 
	{
	    fcolor_t black = { 0.0, 0.0, 0.0 };
	    return black;
	}
    
    return pixel_lerp(img, Tp[0] - Vx(img->origin), Tp[1] - Vy(img->origin));
}

/* Create a new image by applying transformation T to img and 
 * resampling */
img_t *img_resample(img_t *img, trans2D_t *T) {
    int w = img->w, h = img->h;
    int x, y;
    trans2D_t *Tinv = transform_invert(T);
    img_t *Timg = img_new(w, h);

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            double Tp[2];
	    fcolor_t c;

            /* Invert the point (x, y) */
            transform_point(Tinv, x, y, &Tp[0], &Tp[1]);
            
#if 1
            /* Check if the result is in range */
            if (Tp[0] < 0.0 || Tp[1] < 0.0 || Tp[0] > w - 1 || Tp[1] > h - 1) {
                img_set_pixel(Timg, x, y, 
			      background_color.r, background_color.g, background_color.b);

                /* img_nullify_pixel(Timg, x, y); */
            } else {
                /* Apply bilinear interpolation */
                c = pixel_lerp(img, Tp[0], Tp[1]);
                img_set_pixel(Timg, x, y, iround(c.r), iround(c.g), iround(c.b));
            }
#else
	    if (Tp[0] < 0.0) Tp[0] = 0.0;
	    else if (Tp[0] > w - 1) Tp[0] = w - 1;
	    if (Tp[1] < 0.0) Tp[1] = 0.0;
	    else if (Tp[1] > h - 1) Tp[1] = h - 1;

	    c = pixel_lerp(img, Tp[0], Tp[1]);
	    img_set_pixel(Timg, x, y, (int) rint(c.r), (int) rint(c.g), (int) rint(c.b));
#endif
        }
    }

    return Timg;
}

/* Create a new image by applying transformation T to img and 
 * resampling.  Resize the image so that the whole thing fits when
 * transformed. */
img_t *img_resample_bbox(img_t *img, trans2D_t *T) {
    int w = img->w, h = img->h;
    int x, y;
    trans2D_t *Tinv = transform_invert(T);
    int w_new, h_new, i;
    // double x_min = DBL_MAX, x_max = -DBL_MAX, y_min = DBL_MAX, y_max = -DBL_MAX;
    v2_t min = v2_new(DBL_MAX, DBL_MAX);
    v2_t max = v2_new(-DBL_MAX, -DBL_MAX);
    v2_t origin;
    img_t *Timg;

    /* Find the new dimensions of the window */
    v2_t crs[4]; /* Four corners of the original image */

    crs[0] = v2_new(0, 0);
    crs[1] = v2_new(0, h - 1);
    crs[2] = v2_new(w - 1, 0);
    crs[3] = v2_new(w - 1, h - 1);

    for (i = 0; i < 4; i++) {
	crs[i] = v2_add(crs[i], img->origin);
	crs[i] = transform_vector(T, crs[i]);
	min = v2_minimum(min, crs[i]);
	max = v2_maximum(max, crs[i]);
    }

    Vx(min) = floor(Vx(min));
    Vy(min) = floor(Vy(min));

    w_new = iround(floor(Vx(max) - Vx(min) + 1));
    h_new = iround(floor(Vy(max) - Vy(min) + 1));
    origin = min;

    Timg = img_new(w_new, h_new);
    Timg->origin = origin;

    for (y = 0; y < h_new; y++) {
        for (x = 0; x < w_new; x++) {
            double Tp[2];
	    fcolor_t c;

            /* Invert the point (x, y) - trans */
            transform_point(Tinv, x + Vx(origin), y + Vy(origin), &Tp[0], &Tp[1]);
            
#if 1
            /* Check if the result is in range */
            if (Tp[0] < Vx(img->origin) || Tp[1] < Vy(img->origin) || 
		Tp[0] > Vx(img->origin) + w - 1 || Tp[1] > Vy(img->origin) + h - 1) {

		/* pass */

	    } else {
		/* Check if the result is valid */
		int x_f = (int) (Tp[0] - Vx(img->origin));
		int x_c = x_f + 1;
		int y_f = (int) (Tp[1] - Vy(img->origin));
		int y_c = y_f + 1;

		if (img_pixel_is_valid(img, x_f, y_f) ||
		    img_pixel_is_valid(img, x_c, y_f) ||
		    img_pixel_is_valid(img, x_f, y_c) ||
		    img_pixel_is_valid(img, x_c, y_c)) {

		    /* Apply bilinear interpolation */
		    c = pixel_lerp(img, Tp[0] - Vx(img->origin), Tp[1] - Vy(img->origin));
		    img_set_pixel(Timg, x, y, iround(c.r), iround(c.g), iround(c.b));
		} else {
		    // img_nullify_pixel(Timg, x, y);
		}
		
            }
#else
	    if (Tp[0] < 0.0) Tp[0] = 0.0;
	    else if (Tp[0] > w - 1) Tp[0] = w - 1;
	    if (Tp[1] < 0.0) Tp[1] = 0.0;
	    else if (Tp[1] > h - 1) Tp[1] = h - 1;

	    c = pixel_lerp(img, Tp[0], Tp[1]);
	    img_set_pixel(Timg, x, y, c.r, c.g, c.b);
#endif
        }
    }

    /* Add the old origin to the image */
    // Timg->origin = v2_add(origin, img->origin);
    transform_free(Tinv);

    return Timg;
}

img_t *img_fix_radial_distortion(img_t *img, double k1, double k2, double f) {
    int x, y, w = img->w, h = img->h;
    img_t *Timg = img_new(w, h);

    /* Map each pixel to 
     *       x = x' / (1 + k1 * r^2 + k2 * r^4) 
     *       y = y' / (1 + k1 * r^2 + k2 * r^4) */

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    double xf = ((double) x - 0.5 * img->w) / f;
	    double yf = ((double) y - 0.5 * img->h) / f;
	    double rsq = xf * xf + yf * yf;
	    
	    double xt = xf * (1.0 + rsq * (k1 + rsq * k2));
	    double yt = yf * (1.0 + rsq * (k1 + rsq * k2));

	    fcolor_t c;

	    xt = f * xt + 0.5 * img->w;
	    yt = f * yt + 0.5 * img->h;

            if (xt < 0.0 || yt < 0.0 || xt > w - 1 || yt > h - 1) {
                img_set_pixel(Timg, x, y, 
			      background_color.r, background_color.g, background_color.b);
	    } else {
		c = pixel_lerp(img, xt, yt);
                img_set_pixel(Timg, x, y, iround(c.r), iround(c.g), iround(c.b));
	    }
	}
    }

    return Timg;
}


#define MAX_TRIES 256
void img_sample_random_pt(img_t *img, trans2D_t *Tinv, img_t *img_old, int nhood_radius, 
			  int *xout, int *yout, double grad_threshold) 
{
    int x, y;
    int tries = 0;

    int fox = (int) floor(Vx(img->origin));
    int foy = (int) floor(Vy(img->origin));

    while (1) {
	// double grad;

	x = fox + rand() % (img->w - 2 * nhood_radius) + nhood_radius;
	y = foy + rand() % (img->h - 2 * nhood_radius) + nhood_radius;

	/* Make sure the region is valid */
	if (!img_region_is_valid(img, 
				 x - fox - nhood_radius, x - fox + nhood_radius, 
				 y - foy - nhood_radius, y - foy + nhood_radius)) {
	    tries++;
	    continue;
	}

#if 0
	/* Throw out points that don't have a high enough gradient, 
	 * since we want to use high-frequency detail to align the
	 * images */

	grad = img_gradient(img, x - fox, y - foy);
	// printf("%0.3f\n", grad);

	if (grad < grad_threshold && tries < MAX_TRIES) {
	    tries++;
	    continue;
	}
#endif
	
	break;
    }

    *xout = x - fox;
    *yout = y - foy;
}

/* Returns true if the image becomes disconnected under the given
 * homography */
int img_disconnected_under_homography(img_t *img, double *H) {
    /* Check if any point in the img maps to infinity under H.  This
     * is true for p if Hp = [x y 0], so if w is the third row of H, 
     * w^T * p = 0.  We check if this line intersects the image */

    v2_t origin = img->origin;
    int w = img->w, h = img->h;

    v3_t line = v3_new(H[6], H[7], H[8]);

    /* The bottom line is where y = Vy(origin) */
    v3_t bottom_line = v3_cross(v3_new(0.0, Vy(origin), 1.0), 
				v3_new(1.0, Vy(origin), 1.0));
    v3_t top_line = v3_cross(v3_new(0.0, Vy(origin) + h - 1, 1.0),
			     v3_new(1.0, Vy(origin) + h - 1, 1.0));
    v3_t left_line = v3_cross(v3_new(Vx(origin), 0.0, 1.0),
			      v3_new(Vx(origin), 1.0, 1.0));
    v3_t right_line = v3_cross(v3_new(Vx(origin) + w - 1, 0.0, 1.0),
			       v3_new(Vx(origin) + w - 1, 1.0, 1.0));

    v3_t isect;
    
    /* Do four intersection tests */
    isect = v3_homogenize(v3_cross(line, bottom_line));
    if (Vx(isect) > 0 && Vx(isect) < w - 1)
	return 1;
    
    isect = v3_homogenize(v3_cross(line, top_line));
    if (Vx(isect) > 0 && Vx(isect) < w - 1)
	return 1;

    isect = v3_homogenize(v3_cross(line, left_line));
    if (Vy(isect) > 0 && Vy(isect) < h - 1)
	return 1;

    isect = v3_homogenize(v3_cross(line, right_line));
    if (Vy(isect) > 0 && Vy(isect) < h - 1)
	return 1;

    return 0;
}
