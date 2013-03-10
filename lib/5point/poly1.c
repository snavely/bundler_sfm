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

/* poly1.c */
/* 10-degree polynomial in one variable */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "poly1.h"

poly1_t poly1_new3(double a, double b, double c, double d) {
    return poly1_new4(0.0, a, b, c, d);
}

poly1_t poly1_new4(double a, double b, double c, double d, double e) {
	int i;
	poly1_t r;

	for (i = 0; i <= MAX_DEGREE; i++) {
		r.v[i] = 0.0;
	}

        r.v[0] = e;
	r.v[1] = d;
	r.v[2] = c;
	r.v[3] = b;
	r.v[4] = a;

	return r;
}

poly1_t poly1_add(poly1_t a, poly1_t b) {
	int i;
	poly1_t r;

	for (i = 0; i <= MAX_DEGREE; i++) {
		r.v[i] = a.v[i] + b.v[i];
	}

	return r;
}

poly1_t poly1_add3(poly1_t a, poly1_t b, poly1_t c) {
	int i;
	poly1_t r;

	for (i = 0; i <= MAX_DEGREE; i++) {
		r.v[i] = a.v[i] + b.v[i] + c.v[i];
	}

	return r;
}

poly1_t poly1_sub(poly1_t a, poly1_t b) {
	int i;
	poly1_t r;

	for (i = 0; i <= MAX_DEGREE; i++) {
		r.v[i] = a.v[i] - b.v[i];
	}

	return r;
}

int poly1_degree(poly1_t a) {
	int i;

	for (i = MAX_DEGREE; i >= 0; i--) {
            if (fabs(a.v[i]) > 0.0)
                return i;
	}

	return 0;
}

int poly1_mult_check(poly1_t a, poly1_t b) {
	int deg_a = poly1_degree(a);
	int deg_b = poly1_degree(b);

	if (deg_a + deg_b <= MAX_DEGREE)
		return 1;
	else
		return 0;
}

poly1_t poly1_mult(poly1_t a, poly1_t b) {
	poly1_t r;
	int i, j;

#if 0
	if (!poly1_mult_check(a, b)) {
		printf("[poly1_mult] Polynomials cannot be multiplied!\n");
		return poly1_new3(0.0, 0.0, 0.0, 0.0);
	}
#endif

	r = poly1_new3(0.0, 0.0, 0.0, 0.0);

	for (i = 0; i <= MAX_DEGREE; i++) {
		for (j = 0; j <= MAX_DEGREE; j++) {
			int place = i + j;
			if (place > MAX_DEGREE) {
				assert(a.v[i] == 0.0 || b.v[j] == 0.0);
				continue;
			}

			r.v[place] += a.v[i] * b.v[j];
		}
	}

	return r;
}

poly1_t poly1_scale(poly1_t a, double scale) {
    int i;
    poly1_t r;

    for (i = 0; i <= MAX_DEGREE; i++) {
        r.v[i] = scale * a.v[i];
    }

    return r;
}

poly1_t poly1_normalize(poly1_t a) {
    int d = poly1_degree(a);
    if (a.v[d] != 0)
        return poly1_scale(a, 1.0 / a.v[d]);
    else
        return a;
}

double poly1_eval(poly1_t a, double x) {
    double p = 1.0;
    double r = 0.0;
    int i;

    for (i = 0; i <= MAX_DEGREE; i++) {
        r += p * a.v[i];
        p = p * x;
    }

    return r;
}
