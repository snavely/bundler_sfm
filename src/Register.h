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

/* Register.h */
/* Compute relationships between images */

#ifndef __register_h__
#define __register_h__

#include <vector>



#include "keys.h"

enum MotionModel {
    MotionRigid,
    MotionHomography,
};

/* Estimate a transform between two sets of keypoints */
std::vector<int> EstimateTransform(const std::vector<Keypoint> &k1, 
				   const std::vector<Keypoint> &k2, 
				   const std::vector<KeypointMatch> &matches, 
				   MotionModel mm,
				   int nRANSAC, double RANSACthresh, 
				   double *Mout);



#endif /* __register_h__ */
