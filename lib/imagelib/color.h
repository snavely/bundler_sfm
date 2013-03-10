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

/* color.h */

#ifndef __color_h__
#define __color_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>

#ifdef WIN32
#include "types.h"
#endif

typedef struct {
    u_int8_t r, g, b, extra;
} color_t;

typedef struct {
    float r, g, b;
} fcolor_t;

/* Create a new color */
color_t color_new(u_int8_t r, u_int8_t g, u_int8_t b);
fcolor_t fcolor_new(float r, float g, float b);

/* Return the intensity of a color */
double color_intensity(color_t col);
double fcolor_intensity(fcolor_t col);

/* Return the squared distance between two colors */
double color_squared_distance(color_t a, color_t b);

/* Return the squared (weighted) distance between two colors */
double color_squared_weighted_distance(color_t a, color_t b);

/* Return the distance between two colors */
double color_distance(color_t a, color_t b);

/* Return the LUV coordinates of the given color */
void color_RGBtoLUV(unsigned char R, unsigned char G, unsigned char B, 
		   double *L, double *U, double *V);

/* Return the LUV coordinates of the given color */
void color_fRGBtoLUV(double R, double G, double B, 
		     double *L, double *U, double *V);

#ifdef __cplusplus
}
#endif

#endif /* __color_h__ */
