/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
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

/* BoundingBox.cpp */
/* Routines for bounding boxes */

#include <float.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "BoundingBox.h"

#include "defines.h"

void BoundingBox::Print()
{
    printf("[(%0.3f, %0.3f) x (%0.3f, %0.3f)]\n", 
	   m_xmin, m_ymin, m_xmax, m_ymax);
    
}

/* Returns the area of the bounding box */
double BoundingBox::Area()
{   
    return (m_xmax - m_xmin) * (m_ymax - m_ymin);
}

BoundingBox BoundingBox::Intersect(const BoundingBox &bbox) const
{
	return BoundingBox(MAX(m_xmin, bbox.m_xmin), MAX(m_ymin, bbox.m_ymin), 
					   MIN(m_xmax, bbox.m_xmax), MIN(m_ymax, bbox.m_ymax));
}

/* Returns true iff the box contains the given point */
bool BoundingBox::Contains(double x, double y) 
{
    return (x >= m_xmin && x <= m_xmax && y >= m_ymin && y <= m_ymax);
}

bool BoundingBox::Contains(const BoundingBox &bbox)
{
    return (bbox.m_xmin >= m_xmin && 
	    bbox.m_xmax <= m_xmax &&
	    bbox.m_ymin >= m_ymin && 
	    bbox.m_ymax <= m_ymax);    
}

void BoundingBox::Scale(double scale) 
{
    m_xmin *= scale;
    m_xmax *= scale;
    m_ymin *= scale;
    m_ymax *= scale;    
}

/* Create a bounding box for a set of points */
BoundingBox CreateBoundingBox(const std::vector<v2_t> &points) 
{
    int num_points = (int) points.size();
    BoundingBox bb;

    bb.m_xmin = DBL_MAX;
    bb.m_xmax = -DBL_MAX;
    bb.m_ymin = DBL_MAX;
    bb.m_ymax = -DBL_MAX;

    if (num_points == 0) {
	printf("[CreateBoundingBox] No points given!\n");
	// bb.m_xmin = bb.m_xmax = bb.m_ymin = bb.m_ymin = 0.0;
	return bb;
    }

    for (int i = 0; i < num_points; i++) {
	bb.m_xmin = MIN(bb.m_xmin, Vx(points[i]));
	bb.m_xmax = MAX(bb.m_xmax, Vx(points[i]));
	bb.m_ymin = MIN(bb.m_ymin, Vy(points[i]));
	bb.m_ymax = MAX(bb.m_ymax, Vy(points[i]));
    }

    return bb;
}

/* Return the union of two bounding boxes */
BoundingBox BoundingBoxUnion(const BoundingBox &bb1, const BoundingBox &bb2)
{
    BoundingBox bb;
    bb.m_xmin = MIN(bb1.m_xmin, bb2.m_xmin);
    bb.m_xmax = MAX(bb1.m_xmax, bb2.m_xmax);

    bb.m_ymin = MIN(bb1.m_ymin, bb2.m_ymin);
    bb.m_ymax = MAX(bb1.m_ymax, bb2.m_ymax);    

    return bb;
}


