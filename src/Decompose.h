/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
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

/* Decompose.h */
/* Decompose a homography into R and T */

#ifndef __decompose_h__
#define __decompose_h__

#include "vector.h"

/* H is a map from image 2 to image 1 */
bool DecomposeHomography(double *H, double f1, double f2, 
                         double *Ra, double *ta, 
                         double *Rb, double *tb, v2_t p, v2_t q);

bool ComputeFundamentalMatrix(double f1, double f2, 
                              double *R2, double *t2, double *F);

#endif /* __decompose_h__ */
