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

#ifndef __triangle_h__
#define __triangle_h__

#include "vector.h"

/* Compute the nearest point on the segment (s1, s2) to point p */
v3_t segment_pt_nearest_pt(v3_t s1, v3_t s2, v3_t p);

/* Compute the nearest point on the triangle (t1, t2, t3) to point p */
v3_t triangle_pt_nearest_pt(v3_t t1, v3_t t2, v3_t t3, v3_t p);

#endif /* __triangle_h__ */
