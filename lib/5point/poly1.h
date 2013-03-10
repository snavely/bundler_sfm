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

/* poly1.h */
/* 10-degree polynomial in one variable */

#ifndef __POLY1_H__
#define __POLY1_H__

#define MAX_DEGREE 10

typedef struct {
	double v[MAX_DEGREE+1];
} poly1_t;

poly1_t poly1_new3(double a, double b, double c, double d);
poly1_t poly1_new4(double a, double b, double c, double d, double e);
poly1_t poly1_add(poly1_t a, poly1_t b);
poly1_t poly1_add3(poly1_t a, poly1_t b, poly1_t c);
poly1_t poly1_sub(poly1_t a, poly1_t b);
poly1_t poly1_mult(poly1_t a, poly1_t b);
poly1_t poly1_scale(poly1_t a, double scale);
poly1_t poly1_normalize(poly1_t a);

double poly1_eval(poly1_t a, double x);

#endif /* __POLY1_H__ */
