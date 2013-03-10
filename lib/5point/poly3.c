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

/* poly3.c */
/* 3-variable polynomial functions */

#include <stdio.h>
#include <math.h>

#include "poly3.h"

poly3_t poly3_new(double x, double y, double z, double unit) {
	poly3_t p;
	int i;

	for (i = 0; i < NUM_POLY3_COEFFS; i++)
		p.v[i] = 0.0;

	p.v[POLY3_X] = x;
	p.v[POLY3_Y] = y;
	p.v[POLY3_Z] = z;
	p.v[POLY3_UNIT] = unit;

	return p;
}


poly3_t poly3_add(poly3_t a, poly3_t b) {
	int i;
	poly3_t r;

	for (i = 0; i < NUM_POLY3_COEFFS; i++)
		r.v[i] = a.v[i] + b.v[i];

	return r;
}

poly3_t poly3_add3(poly3_t a, poly3_t b, poly3_t c) {
	int i;
	poly3_t r;

	for (i = 0; i < NUM_POLY3_COEFFS; i++)
		r.v[i] = a.v[i] + b.v[i] + c.v[i];

	return r;
}

poly3_t poly3_sub(poly3_t a, poly3_t b) {
	int i;
	poly3_t r;

	for (i = 0; i < NUM_POLY3_COEFFS; i++)
		r.v[i] = a.v[i] - b.v[i];

	return r;	
}

int poly3_degree(poly3_t a) {
	if (a.v[POLY3_X3] != 0.0 || a.v[POLY3_Y3] != 0.0 || a.v[POLY3_Z3] != 0.0 ||
            a.v[POLY3_X2Y] != 0.0 || a.v[POLY3_X2Z] != 0.0 || a.v[POLY3_XY2] != 0.0 || 
            a.v[POLY3_XZ2] != 0.0 || a.v[POLY3_Y2Z] != 0.0 || a.v[POLY3_YZ2] != 0.0 ||
            a.v[POLY3_XYZ] != 0.0)
            return 3;

	if (a.v[POLY3_X2] != 0.0 || a.v[POLY3_Y2] != 0.0 || a.v[POLY3_Z2] != 0.0 ||
            a.v[POLY3_XY] != 0.0 || a.v[POLY3_XZ] != 0.0 || a.v[POLY3_YZ] != 0.0)
            return 2;

	if (a.v[POLY3_X] != 0.0 || a.v[POLY3_Y] != 0.0 || a.v[POLY3_Z] != 0.0)
            return 1;

	return 0;
}

int poly3_mult_check(poly3_t a, poly3_t b) {
	int deg_a = poly3_degree(a);
	int deg_b = poly3_degree(b);

	if (deg_a + deg_b <= 3)
		return 1;
	else
		return 0;
}

poly3_t poly3_mult(poly3_t a, poly3_t b) {
	poly3_t r;

#if 0
	if (!poly3_mult_check(a, b)) {
		printf("[poly3_mult] Polynomials cannot be multiplied!\n");
		return poly3_new(0.0, 0.0, 0.0, 0.0);
	}
#endif

	r.v[POLY3_UNIT] = a.v[POLY3_UNIT] * b.v[POLY3_UNIT];
	r.v[POLY3_X] = a.v[POLY3_X] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X];
	r.v[POLY3_Y] = a.v[POLY3_Y] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Y];
	r.v[POLY3_Z] = a.v[POLY3_Z] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Z];

	r.v[POLY3_XY] = a.v[POLY3_XY] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_XY] + a.v[POLY3_X] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_X];
	r.v[POLY3_XZ] = a.v[POLY3_XZ] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_XZ] + a.v[POLY3_X] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_X];
	r.v[POLY3_YZ] = a.v[POLY3_YZ] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_YZ] + a.v[POLY3_Y] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_Y];

	r.v[POLY3_X2] = a.v[POLY3_X2] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X2] + a.v[POLY3_X] * b.v[POLY3_X];
	r.v[POLY3_Y2] = a.v[POLY3_Y2] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Y2] + a.v[POLY3_Y] * b.v[POLY3_Y];
	r.v[POLY3_Z2] = a.v[POLY3_Z2] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Z2] + a.v[POLY3_Z] * b.v[POLY3_Z];

	r.v[POLY3_X2Y] = a.v[POLY3_X2Y] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X2Y] + 
						a.v[POLY3_X2] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_X2] + 
						a.v[POLY3_XY] * b.v[POLY3_X] + a.v[POLY3_X] * b.v[POLY3_XY];

	r.v[POLY3_X2Z] = a.v[POLY3_X2Z] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X2Z] + 
						a.v[POLY3_X2] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_X2] + 
						a.v[POLY3_XZ] * b.v[POLY3_X] + a.v[POLY3_X] * b.v[POLY3_XZ];

	r.v[POLY3_XY2] = a.v[POLY3_XY2] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_XY2] + 
						a.v[POLY3_X] * b.v[POLY3_Y2] + a.v[POLY3_Y2] * b.v[POLY3_X] + 
						a.v[POLY3_XY] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_XY];

	r.v[POLY3_Y2Z] = a.v[POLY3_Y2Z] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Y2Z] + 
						a.v[POLY3_Y2] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_Y2] + 
						a.v[POLY3_YZ] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_YZ];

	r.v[POLY3_XZ2] = a.v[POLY3_XZ2] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_XZ2] + 
						a.v[POLY3_X] * b.v[POLY3_Z2] + a.v[POLY3_Z2] * b.v[POLY3_X] + 
						a.v[POLY3_XZ] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_XZ];

	r.v[POLY3_YZ2] = a.v[POLY3_YZ2] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_YZ2] + 
						a.v[POLY3_Y] * b.v[POLY3_Z2] + a.v[POLY3_Z2] * b.v[POLY3_Y] + 
						a.v[POLY3_YZ] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_YZ];

	r.v[POLY3_XYZ] = a.v[POLY3_XYZ] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_XYZ] + 
						a.v[POLY3_XY] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_XY] + 
						a.v[POLY3_XZ] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_XZ] + 
						a.v[POLY3_YZ] * b.v[POLY3_X] + a.v[POLY3_X] * b.v[POLY3_YZ];

	r.v[POLY3_X3] = a.v[POLY3_X3] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X3] + 
					a.v[POLY3_X2] * b.v[POLY3_X] + a.v[POLY3_X] * b.v[POLY3_X2];

	r.v[POLY3_Y3] = a.v[POLY3_Y3] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Y3] + 
					a.v[POLY3_Y2] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_Y2];

	r.v[POLY3_Z3] = a.v[POLY3_Z3] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Z3] + 
					a.v[POLY3_Z2] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_Z2];

	return r;
}

poly3_t poly3_mult11(poly3_t a, poly3_t b) {
	poly3_t r;
        int i;

        for (i = 0; i < NUM_POLY3_COEFFS; i++) 
            r.v[i] = 0.0;

	r.v[POLY3_UNIT] = a.v[POLY3_UNIT] * b.v[POLY3_UNIT];
	r.v[POLY3_X] = a.v[POLY3_X] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X];
	r.v[POLY3_Y] = a.v[POLY3_Y] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Y];
	r.v[POLY3_Z] = a.v[POLY3_Z] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Z];

	r.v[POLY3_XY] = a.v[POLY3_X] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_X];
	r.v[POLY3_XZ] = a.v[POLY3_X] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_X];
	r.v[POLY3_YZ] = a.v[POLY3_Y] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_Y];

	r.v[POLY3_X2] = a.v[POLY3_X] * b.v[POLY3_X];
	r.v[POLY3_Y2] = a.v[POLY3_Y] * b.v[POLY3_Y];
	r.v[POLY3_Z2] = a.v[POLY3_Z] * b.v[POLY3_Z];

	return r;
}


poly3_t poly3_mult21(poly3_t a, poly3_t b) {
	poly3_t r;

	r.v[POLY3_UNIT] = a.v[POLY3_UNIT] * b.v[POLY3_UNIT];
	r.v[POLY3_X] = a.v[POLY3_X] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_X];
	r.v[POLY3_Y] = a.v[POLY3_Y] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Y];
	r.v[POLY3_Z] = a.v[POLY3_Z] * b.v[POLY3_UNIT] + a.v[POLY3_UNIT] * b.v[POLY3_Z];

	r.v[POLY3_XY] = a.v[POLY3_XY] * b.v[POLY3_UNIT] + a.v[POLY3_X] * b.v[POLY3_Y] + a.v[POLY3_Y] * b.v[POLY3_X];
	r.v[POLY3_XZ] = a.v[POLY3_XZ] * b.v[POLY3_UNIT] + a.v[POLY3_X] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_X];
	r.v[POLY3_YZ] = a.v[POLY3_YZ] * b.v[POLY3_UNIT] + a.v[POLY3_Y] * b.v[POLY3_Z] + a.v[POLY3_Z] * b.v[POLY3_Y];

	r.v[POLY3_X2] = a.v[POLY3_X2] * b.v[POLY3_UNIT] + a.v[POLY3_X] * b.v[POLY3_X];
	r.v[POLY3_Y2] = a.v[POLY3_Y2] * b.v[POLY3_UNIT] + a.v[POLY3_Y] * b.v[POLY3_Y];
	r.v[POLY3_Z2] = a.v[POLY3_Z2] * b.v[POLY3_UNIT] + a.v[POLY3_Z] * b.v[POLY3_Z];

	r.v[POLY3_X2Y] = 
            a.v[POLY3_X2] * b.v[POLY3_Y] + 
            a.v[POLY3_XY] * b.v[POLY3_X];

	r.v[POLY3_X2Z] = 
            a.v[POLY3_X2] * b.v[POLY3_Z] + 
            a.v[POLY3_XZ] * b.v[POLY3_X];

	r.v[POLY3_XY2] = 
            a.v[POLY3_Y2] * b.v[POLY3_X] + 
            a.v[POLY3_XY] * b.v[POLY3_Y];

	r.v[POLY3_Y2Z] = 
            a.v[POLY3_Y2] * b.v[POLY3_Z] + 
            a.v[POLY3_YZ] * b.v[POLY3_Y];

	r.v[POLY3_XZ2] = 
            a.v[POLY3_Z2] * b.v[POLY3_X] + 
            a.v[POLY3_XZ] * b.v[POLY3_Z];

	r.v[POLY3_YZ2] = 
            a.v[POLY3_Z2] * b.v[POLY3_Y] + 
            a.v[POLY3_YZ] * b.v[POLY3_Z];

	r.v[POLY3_XYZ] = 
            a.v[POLY3_XY] * b.v[POLY3_Z] +
            a.v[POLY3_XZ] * b.v[POLY3_Y] +
            a.v[POLY3_YZ] * b.v[POLY3_X];

	r.v[POLY3_X3] = a.v[POLY3_X2] * b.v[POLY3_X];
	r.v[POLY3_Y3] = a.v[POLY3_Y2] * b.v[POLY3_Y];
	r.v[POLY3_Z3] = a.v[POLY3_Z2] * b.v[POLY3_Z];

	return r;
}

poly3_t poly3_scale(poly3_t a, double scale) {

	poly3_t r;
	int i;

	for (i = 0; i < NUM_POLY3_COEFFS; i++) 
		r.v[i] = scale * a.v[i];

	return r;
}

double poly3_eval(poly3_t a, double x, double y, double z) {
	double r = a.v[POLY3_UNIT];

	r += a.v[POLY3_X] * x;
	r += a.v[POLY3_Y] * y;
	r += a.v[POLY3_Z] * z;

	r += a.v[POLY3_XY] * x * y;
	r += a.v[POLY3_XZ] * x * z;
	r += a.v[POLY3_YZ] * y * z;

	r += a.v[POLY3_X2] * x * x;
	r += a.v[POLY3_Y2] * y * y;
	r += a.v[POLY3_Z2] * z * z;

	r += a.v[POLY3_X2Y] * x * x * y;
	r += a.v[POLY3_X2Z] * x * x * z;
	r += a.v[POLY3_XY2] * x * y * y;
	r += a.v[POLY3_Y2Z] * y * y * z;
	r += a.v[POLY3_XZ2] * x * z * z;
	r += a.v[POLY3_YZ2] * y * z * z;
	r += a.v[POLY3_XYZ] * x * y * z;

	r += a.v[POLY3_X3] * x * x * x;
	r += a.v[POLY3_Y3] * y * y * y;
	r += a.v[POLY3_Z3] * z * z * z;

	return r;
}

double poly3_get(poly3_t a, poly3_coeff_t idx) {
	return a.v[idx];
}
