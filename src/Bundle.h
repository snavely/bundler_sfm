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

/* Bundle.h */
/* Utility functions for bundle adjustment */

#ifndef __bundle_h__
#define __bundle_h__

#include "ImageData.h"
#include "Geometry.h"

#include "sfm.h"

/* Compute the angle between two rays */
double ComputeRayAngle(v2_t p, v2_t q, 
		       const camera_params_t &cam1, 
		       const camera_params_t &cam2);

/* Check cheirality for a camera and a point */
bool CheckCheirality(v3_t p, const camera_params_t &camera);

void GetIntrinsics(const camera_params_t &camera, double *K);

double GetCameraDistance(camera_params_t *c1, camera_params_t *c2);

/* Use a 180 rotation to fix up the intrinsic matrix */
void FixIntrinsics(double *P, double *K, double *R, double *t);

void InitializeCameraParams(const ImageData &data, camera_params_t &camera);

#endif /* __bundle_h__ */
