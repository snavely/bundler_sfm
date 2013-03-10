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

/* triangle.c */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "triangle.h"
#include "vector.h"

#define SGN(x) ((x) < 0 ? -1 : 1)

/* Returns 1 if p lies on the segment defined by points s1 and s2 */
int pt_on_segment(v3_t p, v3_t s1, v3_t s2) {
    v3_t d1 = v3_unit(v3_sub(s1, p));
    v3_t d2 = v3_unit(v3_sub(s2, p));
    return (v3_dotp(d1, d2) == -1.0);
}

/* Returns 1 if p lies inside the triangle defined by t1, t2, t3.
 * Point p is assumed to be coplanar with the triangle */
int pt_in_triangle(v3_t p, v3_t t1, v3_t t2, v3_t t3) {
    v3_t c1, c2, c3;
    double d1, d2, d3;

    /* Check if the point lies on any of the segments */
    if (pt_on_segment(p, t1, t2))
        return 1;
    if (pt_on_segment(p, t1, t3))
        return 1;
    if (pt_on_segment(p, t2, t3))
        return 1;

    /* Check that the cross products all point in the same direction */    
    c1 = v3_cross(v3_sub(p, t1), v3_sub(t2, t1));
    c2 = v3_cross(v3_sub(p, t2), v3_sub(t3, t2));
    c3 = v3_cross(v3_sub(p, t3), v3_sub(t1, t3));
    
    d1 = v3_dotp(c1, c2);
    d2 = v3_dotp(c1, c3);
    d3 = v3_dotp(c2, c3);
    
    return (SGN(d1) == SGN(d2) && SGN(d1) == SGN(d3));
}

v3_t segment_pt_nearest_pt(v3_t s1, v3_t s2, v3_t p) {
    v3_t diff = v3_sub(s2, s1);
    double t = v3_dotp(diff, v3_sub(p, s1)) / v3_magsq(diff);
    
    if (t >= 0.0 && t <= 1.0)
	return v3_add(s1, v3_scale(t, v3_sub(s2, s1)));
    else if (t < 0.0)
	return s1;
    else
	return s2;
}

v3_t triangle_pt_nearest_pt(v3_t t1, v3_t t2, v3_t t3, v3_t p) {
    v3_t norm, v12, v13, q;
    double dotp;
    double t;
    double d1, d2, d3;

    /* Compute the normal to the plane of the triangle 
     * using the cross product */
    v12 = v3_unit(v3_sub(t2, t1));
    v13 = v3_unit(v3_sub(t3, t1));
    
    /* Make sure v12 and v13 are linearly independent */
    dotp = v3_dotp(v12, v13);
    if (fabs(dotp) == 1.0) {
        /* t1, t2, t3 are colinear */
        printf("Error for now; make line test later\n");
    }

    norm = v3_unit(v3_cross(v12, v13));

    /* Find the intersection of p + t * norm with the plane 
     * of the triangle */
    t = v3_dotp(v3_sub(p, t1), norm);

    /* q is the projection of p on the plane */
    q = v3_sub(p, v3_scale(t, norm));
    
    /* Check if q lies within the triangle */
    if (pt_in_triangle(q, t1, t2, t3))
        return q;

    /* We know that the nearest point must be on the boundary
     * of the triangle.  Find out which line segment we should 
     * check */
    d1 = v3_magsq(v3_sub(p, t1));
    d2 = v3_magsq(v3_sub(p, t2));
    d3 = v3_magsq(v3_sub(p, t3));
    
    if (d1 >= d2 && d1 >= d3)
        return segment_pt_nearest_pt(t2, t3, p);
    else if (d2 >= d1 && d2 >= d3)
        return segment_pt_nearest_pt(t1, t3, p);
    else if (d3 >= d1 && d3 >= d2)
        return segment_pt_nearest_pt(t1, t2, p);
    else
        printf("Contradiction reached\n");

    return v3_new(0.0, 0.0, 0.0);
}
