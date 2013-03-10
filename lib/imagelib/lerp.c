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

/* lerp.c */
/* Interpolation routines */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <sys/types.h>

#include "lerp.h"

#ifdef WIN32
#include "types.h"
#endif

/* Apply the permutation p to the array x */
static void permute(int n, int *x, int *p) {
    int i;
    int *xnew = malloc(sizeof(int) * n);
    
    for (i = 0; i < n; i++)
	xnew[i] = x[p[i]];

    memcpy(x, xnew, sizeof(int) * n);
    free(xnew);
}

static double eval(int n, double *f, int *coords) {
    u_int32_t idx = 0;
    int i;

    for (i = 0; i < n; i++) {
	idx += coords[i] * (1LL << i);
    }

    return f[idx];
}

/* Interpolate 0 <= x <= 1 between the given function values at points
 * on the [0,1]^n grid.
 *
 * x: an array of length n
 * f: an array of length 2^n
 *
 * To see how this works, see 
 *    http://osl.iu.edu/~tveldhui/papers/MAScThesis/node33.html */
double lerp(int n, double *xin, double *f) {
    double *x = malloc(sizeof(double) * n);
    double *xsort = malloc(sizeof(double) * n);
    double sum = 0.0;
    int *p = malloc(sizeof(int) * n);
    int *pinv = malloc(sizeof(int) * n);
    int *f1, *f2;
    int i, j;
    
    memcpy(xsort, xin, sizeof(double) * n);
    memcpy(x, xin, sizeof(double) * n);

    /* Find the permutation p that sorts x */
    for (i = 0; i < n; i++) {
	double max = xsort[0];
	int maxidx = 0;
	for (j = 1; j < n; j++) {
	    if (xsort[j] > max) {
		max = xsort[j];
		maxidx = j;
	    }
	}

	p[i] = maxidx;
	pinv[maxidx] = i;
	x[i] = max;
	xsort[maxidx] = -DBL_MAX;
    }

    sum = f[0];

    f1 = malloc(sizeof(int) * n);
    f2 = malloc(sizeof(int) * n);

    /* Interpolate! */
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    f1[j] = (i < j) ? 0 : 1;
	    f2[j] = (i <= j) ? 0 : 1;
	}
	
	permute(n, f1, pinv);
	permute(n, f2, pinv);
	    
	sum += (eval(n, f, f1) - eval(n, f, f2)) * x[i];
    }

    free(f1);
    free(f2);
    free(xsort);
    free(p);
    free(pinv);
    free(x);

    return sum;
}

/* Bilinear interpolation on a 2D function */
double func_lerp(int w, int h, double *f, double x, double y) {
    int x_f = (int) floor(x);
    int y_f = (int) floor(y);

    double fsub[4] = 
	{ f[y_f * w + x_f],     f[y_f * w + x_f + 1],
	  f[(y_f+1) * w + x_f], f[(y_f+1) * w + x_f + 1] };

    double xin[2] = { x - x_f, y - y_f };

    return nlerp2(xin, fsub);
}

/* Bilinear interpolation on a 2D function */
double func_lerpf(int w, int h, float *f, double x, double y) {
    int x_f = (int) floor(x);
    int y_f = (int) floor(y);

    double fsub[4] = 
	{ f[y_f * w + x_f],     f[y_f * w + x_f + 1],
	  f[(y_f+1) * w + x_f], f[(y_f+1) * w + x_f + 1] };

    double xin[2] = { x - x_f, y - y_f };

    return nlerp2(xin, fsub);
}

/* Special case lerp for 2 dimensions */
double lerp2(double x, double y, int *f) {
    if (x >= y)
	return ((double) f[0]) + ((double) (f[1] - f[0])) * x + ((double) (f[3] - f[1])) * y;
    else
	return ((double) f[0]) + ((double) (f[2] - f[0])) * y + ((double) (f[3] - f[2])) * x;
}

double nlerp2(double *xin, double *f) {
    return (1.0 - xin[1]) * ((1.0 - xin[0]) * f[0] + xin[0] * f[1]) + 
	xin[1] * ((1.0 - xin[0]) * f[2] + xin[0] * f[3]);
}

double nlerp(int n, double *xin, double *f) {
    if (n == 1)
        return (1.0 - xin[0]) * f[0] + xin[0] * f[1];
    else
        return (1.0 - xin[0]) * nlerp(n - 1, xin+1, f) +
            xin[0] * nlerp(n - 1, xin+1, f + (1LL << (n - 1)));
}

static double cubic_weight(double x) {
    double x_2 = x + 2.0;
    double x_1 = x + 1.0;
    double x_0 = x + 0.0;
    double x_minus1 = x - 1.0;
    
    if (x_2 < 0.0)
	x_2 = 0.0;
    if (x_1 < 0.0)
	x_1 = 0.0;
    if (x_0 < 0.0)
	x_0 = 0.0;
    if (x_minus1 < 0.0)
	x_minus1 = 0.0;
    
    return (x_2 * x_2 * x_2 - 
	    4.0 * x_1 * x_1 * x_1 + 
	    6.0 * x_0 * x_0 * x_0 - 
	    4.0 * x_minus1 * x_minus1 * x_minus1) / 6.0;
}

double bicubic_interpolate_2D(double x, double y, double *f) {
    double sum = 0.0;
    int m, n;
    
    for (m = 0; m <= 3; m++) {
	for (n = 0; n <= 3; n++) {
	    double fc = f[m * 4 + n];
	    sum += fc * cubic_weight(m - 1 - y) * cubic_weight(n - 1 - x);
	}
    }

    return sum;
}

double cuberp(double x, double *f) {
    double M[4][4] = 
	{ { -1.0,  3.0, -2.0, 0.0 },
	  {  3.0, -6.0, -3.0, 1.0 },
	  { -3.0,  3.0,  6.0, 0.0 },
	  {  1.0,  0.0, -1.0, 0.0 } };

    double coeff[4] = { 0.0, 0.0, 0.0, 0.0 };
    int i, j;
    
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    coeff[i] += f[i] * M[j][i];
	}

	// coeff[i] /= 6.0;
    }

    return coeff[3] + x * (coeff[2] + x * (coeff[1] + x * coeff[0]));
}

double bicuberp(double x, double y, double *f) {
    double fcol[4];
    int i;
    
    for (i = 0; i < 4; i++)
	fcol[i] = cuberp(x, f + 4 * i);

    return cuberp(y, fcol);
}

// #define TEST

#ifdef TEST
int main(int argc, char **argv) {
    int n = 2;
    double f[4] = { 0.0, 1.0, 2.0, 3.0 };
    double x[2] = { atof(argv[1]), atof(argv[2]) };
    double result;
    
    result = lerp(n, x, f);
    printf("Result is: %f\n", result);

    result = lerp2(x[0], x[1], f);
    printf("Result is: %f\n", result);

    return 0;
}
#endif
