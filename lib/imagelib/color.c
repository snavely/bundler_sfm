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

/* color.c */

#include <math.h>

#include "color.h"
#include "defines.h"

/* Create a new color */
color_t color_new(u_int8_t r, u_int8_t g, u_int8_t b) 
{
    color_t c = { r, g, b, 0x0 };
    return c;
}

fcolor_t fcolor_new(float r, float g, float b) 
{
    fcolor_t c = { r, g, b };
    return c;
}

/* Return the intensity of a color */
double color_intensity(color_t col) {   
    return 0.30 * col.r + 0.59 * col.g + 0.11 * col.b;
    // return 1.0 * col.r;
}

double fcolor_intensity(fcolor_t col) {
    return 0.30 * col.r + 0.59 * col.g + 0.11 * col.b;
}

/* Return the squared distance between two colors */
double color_squared_distance(color_t a, color_t b) {
    return (SQ(b.r - a.r) + SQ(b.g - a.g) + SQ(b.b - a.b));
}

/* Return the squared (weighted) distance between two colors */
double color_squared_weighted_distance(color_t a, color_t b) {
    return (0.30 * SQ(b.r - a.r) + 0.59 * SQ(b.g - a.g) + 0.11 * SQ(b.b - a.b));
}

/* Return the distance between two colors */
double color_distance(color_t a, color_t b) {
    return sqrt(color_squared_distance(a, b));
}

    

const static double Xn= 0.95050;
const static double Yn= 1.00000;
const static double Zn= 1.08870;

const static double Un_prime= 0.19784977571475;
const static double Vn_prime= 0.46834507665248;
const static double Lt= 0.008856;

/* RGB to LUV conversion */
const static double XYZ[3][3] = { 
    { 0.4125,  0.3576,  0.1804 },
    { 0.2125,  0.7154,  0.0721 },
    { 0.0193,  0.1192,  0.9502 }
};

/* LUV to RGB conversion */
const static double RGB[3][3] = {
    { 3.2405, -1.5371, -0.4985 },
    { -0.9693,  1.8760,  0.0416 },
    { 0.0556, -0.2040,  1.0573 }
};

double sRGBConvert(double c) 
{
    double cd = c / 255.0;

    if (cd <= 0.04045)
	return cd / 12.92;
    else
	return pow((cd + 0.055) / 1.055, 2.4);    
}

/* Return the LUV coordinates of the given color */
void color_fRGBtoLUV(double R, double G, double B, 
		     double *L, double *U, double *V) 
{
    double x, y, z, L0, u_prime, v_prime, constant;

    /* Convert RGV to XYZ */
    double r = sRGBConvert(R);
    double g = sRGBConvert(G);
    double b = sRGBConvert(B);
    
    x = XYZ[0][0] * r + XYZ[0][1] * g + XYZ[0][2] * b;
    y = XYZ[1][0] * r + XYZ[1][1] * g + XYZ[1][2] * b;
    z = XYZ[2][0] * r + XYZ[2][1] * g + XYZ[2][2] * b;

    /* Convert XYZ to LUV... */

    /* Compute L* */
    L0= y / (Yn);
    
    if (L0 > Lt)
	*L = (float)(116.0 * (pow(L0, 1.0/3.0)) - 16.0);
    else
	*L = (float)(903.3 * L0);

    /* Compute u_prime and v_prime */
    constant= x + 15 * y + 3 * z;
    
    if(constant != 0) {
	u_prime= (4 * x) / constant;
	v_prime = (9 * y) / constant;
    } else { 
	u_prime= 4.0;    
	v_prime= 9.0/15.0;
    }
    
    /* Compute u* and v* */
    *U = (float) (13 * *L * (u_prime - Un_prime));
    *V = (float) (13 * *L * (v_prime - Vn_prime));
}

/* Return the LUV coordinates of the given color */
void color_RGBtoLUV(unsigned char R, unsigned char G, unsigned char B, 
		   double *L, double *U, double *V) 
{
    double x, y, z, L0, u_prime, v_prime, constant;

    /* Convert RGV to XYZ */
    double r = sRGBConvert((double) R);
    double g = sRGBConvert((double) G);
    double b = sRGBConvert((double) B);
    
    x = XYZ[0][0] * r + XYZ[0][1] * g + XYZ[0][2] * b;
    y = XYZ[1][0] * r + XYZ[1][1] * g + XYZ[1][2] * b;
    z = XYZ[2][0] * r + XYZ[2][1] * g + XYZ[2][2] * b;

    /* Convert XYZ to LUV... */

    /* Compute L* */
    L0= y / (Yn);
    
    if (L0 > Lt)
	*L = (float)(116.0 * (pow(L0, 1.0/3.0)) - 16.0);
    else
	*L = (float)(903.3 * L0);

    /* Compute u_prime and v_prime */
    constant= x + 15 * y + 3 * z;
    
    if(constant != 0) {
	u_prime= (4 * x) / constant;
	v_prime = (9 * y) / constant;
    } else { 
	u_prime= 4.0;    
	v_prime= 9.0/15.0;
    }
    
    /* Compute u* and v* */
    *U = (float) (13 * *L * (u_prime - Un_prime));
    *V = (float) (13 * *L * (v_prime - Vn_prime));
}
