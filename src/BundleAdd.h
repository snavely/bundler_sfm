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

/* BundleAdd.h */

#ifndef __bundle_add_h__
#define __bundle_add_h__

#include "vector.h"
#include "sfm.h"

/* Triangulate two points */
v3_t Triangulate(v2_t p, v2_t q, 
                 camera_params_t c1, camera_params_t c2, 
                 double &proj_error, bool &in_front, double &angle,
                 bool explicit_camera_centers);

#endif /* __bundle_add_h__ */
