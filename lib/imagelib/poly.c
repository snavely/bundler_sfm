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

/* poly.c */
/* Routines for dealing with polynomials */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "poly.h"
// #include "root/root.h"

/* Allocate and return a new polynomial with degree n */
poly_t *poly_new(int n) {
    if (n <= 0) {
	printf("[poly_new] Degree of a polynomial must be >= 0\n");
	return NULL;
    } else {
	poly_t *p = (poly_t *)malloc(sizeof(poly_t));
	int i;

	p->n = n;
	p->coeffs = (double *)malloc(sizeof(double) * (n+1));

	for (i = 0; i < n + 1; i++) {
	    p->coeffs[i] = 0.0;
	}

	return p;
    }
}

/* Free a polynomial */
void poly_free(poly_t *p) {
    free(p->coeffs);
    free(p);
}

/* Get a coefficient from a polynomial */
double poly_get_coeff(poly_t *p, int idx) {
    if (idx > p->n) 
	return 0.0;
    else 
	return p->coeffs[idx];
}

/* Set a polynomial coefficient */
void poly_set_coeff(poly_t *p, int idx, double c) {
    if (idx > p->n) {
	printf("[poly_set_coeff] Tried to set out-of-bounds coefficient\n");
	return;
    } else {
	p->coeffs[idx] = c;
    }
}

/* Add two polynomials */
poly_t *poly_sum(poly_t *p, poly_t *q) {
    int n = MAX(p->n, q->n);
    poly_t *r = poly_new(n);
    int i;

    for (i = 0; i < n + 1; i++) {
	double c = poly_get_coeff(p, i) + poly_get_coeff(q, i);
	poly_set_coeff(r, i, c);
    }

    return r;
}

/* Subtract two polynomials */
poly_t *poly_diff(poly_t *p, poly_t *q) {
    int n = MAX(p->n, q->n);
    poly_t *r = poly_new(n);
    int i;

    for (i = 0; i < n + 1; i++) {
	double c = poly_get_coeff(p, i) - poly_get_coeff(q, i);
	poly_set_coeff(r, i, c);
    }

    return r;
}

/* Multiply two polynomials and return the result */
poly_t *poly_product(poly_t *p, poly_t *q) {
    int n = p->n + q->n;
    poly_t *r = poly_new(n);
    int i, j;
    
    /* Multiply */
    for (i = 0; i < p->n + 1; i++)
	for (j = 0; j < q->n + 1; j++)
	    r->coeffs[i+j] += p->coeffs[i] * q->coeffs[j];

    return r;
}

/* Compute the derivative of a polynomial */
poly_t *poly_deriv(poly_t *p) {
    if (p->n == 1) {
	/* Derivative is zero */
	return poly_new(1);
    } else {
	poly_t *p_d = poly_new(p->n - 1);
	int i = 0;

	for (i = 0; i < p->n; i++) {
	    p_d->coeffs[i] = (i + 1) * p->coeffs[i + 1];
	}

	return p_d;
    }
}

/* Evaluate the polynomial at the given point */
double poly_eval(poly_t *p, double x) {
    int i;
    int n = p->n;
    double accum = p->coeffs[n];

    for (i = n - 1; i >= 0; i--) {
	accum = accum * x + p->coeffs[i];
    }

    return accum;
}

/* Use Newton's method to find the root of a polynomial.  The
 * polynomial is represented by n coefficients.  The coefficients are
 * in order from low power to high power */
double poly_find_root(poly_t *p, double x0, double thres) {
    /* Evaluate the polynomial at the current x */
    int trials = 0;
    double x_curr = x0;
    double p_x;
    double roots[10];

    poly_t *p_deriv = poly_deriv(p);

    p_x = poly_eval(p, x_curr);

#define MAX_TRIALS 1024
    /* Iterate until convergence */
    while (fabs(p_x) > thres && trials < MAX_TRIALS) {
	double dp_dx = poly_eval(p_deriv, x_curr);
	
	x_curr = x_curr - p_x / dp_dx;
	p_x = poly_eval(p, x_curr);

	trials++;
    }

#if 0
    if (trials >= MAX_TRIALS) {
	printf("[poly_find_root] Maximum number of trials exceeded.\n");
    }
#endif

    poly_free(p_deriv);

    // find_roots(p->n, p->coeffs, roots);

    return x_curr;
}
