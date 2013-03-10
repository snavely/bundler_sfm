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

/* util.h */
/* Various utility functions */

#ifndef __util_h__
#define __util_h__

#ifdef __cplusplus
extern "C" {
#endif

/* Returns the floor of the log (base 2) of the given number */
int ilog2(int n);

/* Returns a random double between 0.0 and 1.0 */
double rand_unit();

/* Returns a random double between min and max */
double rand_double(double min, double max);

/* Clamps the value x to lie within min and max */
double clamp(double x, double min, double max);

/* Returns true if n is a power of two */
int is_power_of_two(int n);

/* Returns the smallest power of two no smaller than the input */
int least_larger_power_of_two(int n);

/* Return the closest integer to x, rounding up */
int iround(double x);

#ifdef __cplusplus
}
#endif

#endif /* __util_h__ */
