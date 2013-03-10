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

/* poly.h */
/* Routines for dealing with polynomials */

#ifndef __poly_h__
#define __poly_h__

typedef struct {
    int n;           /* Degree of the polynomial */
    double *coeffs;  /* Array of coefficients */
} poly_t;

/* Allocate and return a new polynomial with degree n */
poly_t *poly_new(int n);

/* Free a polynomial */
void poly_free(poly_t *p);

/* Get a coefficient from a polynomial */
double poly_get_coeff(poly_t *p, int idx);

/* Set a polynomial coefficient */
void poly_set_coeff(poly_t *p, int idx, double c);

/* Add two polynomials */
poly_t *poly_sum(poly_t *p, poly_t *q);

/* Subtract two polynomials */
poly_t *poly_diff(poly_t *p, poly_t *q);

/* Multiply two polynomials and return the result */
poly_t *poly_product(poly_t *p, poly_t *q);

/* Compute the derivative of a polynomial */
poly_t *poly_deriv(poly_t *p);

/* Evaluate the polynomial at the given point */
double poly_eval(poly_t *p, double x);

/* Use Newton's method to find the root of a polynomial.  The
 * polynomial is represented by n coefficients.  The coefficients are
 * in order from low power to high power */
double poly_find_root(poly_t *p, double x0, double thres);

#endif /* __poly_h__ */
