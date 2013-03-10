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

/* Register.cpp */
/* Compute relationships between images */

#include <assert.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "keys.h"
#include "Register.h"



#include "homography.h"
#include "horn.h"
#include "matrix.h"
#include "tps.h"
#include "vector.h"

static int CountInliers(const std::vector<Keypoint> &k1, 
			const std::vector<Keypoint> &k2, 
			std::vector<KeypointMatch> matches,
			double *M, double thresh, std::vector<int> &inliers);

static int LeastSquaresFit(const std::vector<Keypoint> &k1, 
			   const std::vector<Keypoint> &k2, 
			   std::vector<KeypointMatch> matches, MotionModel mm,
			   const std::vector<int> &inliers, double *M);

/* Estimate a transform between two sets of keypoints */
std::vector<int> EstimateTransform(const std::vector<Keypoint> &k1, 
				   const std::vector<Keypoint> &k2, 
				   const std::vector<KeypointMatch> &matches, 
				   MotionModel mm,
				   int nRANSAC, double RANSACthresh, 
				   double *Mout) 
{
    int min_matches = -1;
    switch (mm) {
	case MotionRigid:
	    min_matches = 3;
	    break;
	case MotionHomography:
	    min_matches = 4;
	    break;
    }

    int *match_idxs = new int[min_matches];

    int num_matches = (int) matches.size();
    int max_inliers = 0;
    double Mbest[9];
    
    if (num_matches < min_matches) {
	std::vector<int> empty;
	printf("Cannot estimate rigid transform\n");
	return empty;
    }

    v3_t *r_pts = new v3_t[min_matches];
    v3_t *l_pts = new v3_t[min_matches];
    double *weight = new double[min_matches];

    for (int round = 0; round < nRANSAC; round++) {
	for (int i = 0; i < min_matches; i++) {
	    bool found;
	    int idx;
	    
	    do {
		found = true;
		idx = rand() % num_matches;
		
		for (int j = 0; j < i; j++) {
		    if (match_idxs[j] == idx) {
			found = false;
			break;
		    }
		}
	    } while (!found);

	    match_idxs[i] = idx;
	}

	/* Solve for the motion */
		
	for (int i = 0; i < min_matches; i++) {
	    int idx1 = matches[match_idxs[i]].m_idx1;
	    int idx2 = matches[match_idxs[i]].m_idx2;
	    
	    Vx(l_pts[i]) = k1[idx1].m_x;
	    Vy(l_pts[i]) = k1[idx1].m_y;
	    Vz(l_pts[i]) = 1.0;
		    
	    Vx(r_pts[i]) = k2[idx2].m_x;
	    Vy(r_pts[i]) = k2[idx2].m_y;
	    Vz(r_pts[i]) = 1.0;

	    weight[i] = 1.0;
	}

	double Mcurr[9];

	switch (mm) {
	    case MotionRigid: {
		double R[9], T[9], Tout[9], scale;
		align_horn(min_matches, r_pts, l_pts, R, T, Tout, &scale, weight);
		memcpy(Mcurr, Tout, 9 * sizeof(double));
		break;
	    }
		
	    case MotionHomography: {
		align_homography(min_matches, r_pts, l_pts, Mcurr, 0);
		break;
	    }
	}
		

	std::vector<int> inliers;
	int num_inliers = CountInliers(k1, k2, matches, Mcurr, 
				       RANSACthresh, inliers);

	if (num_inliers > max_inliers) {
	    max_inliers = num_inliers;
	    memcpy(Mbest, Mcurr, 9 * sizeof(double));
	}
    }

    std::vector<int> inliers;
    CountInliers(k1, k2, matches, Mbest, RANSACthresh, inliers);
    memcpy(Mout, Mbest, 9 * sizeof(double));
    LeastSquaresFit(k1, k2, matches, mm, inliers, Mout);

    // memcpy(Mout, Mbest, 9 * sizeof(double));

    delete [] match_idxs;
    delete [] r_pts;
    delete [] l_pts;
    delete [] weight;

    return inliers;
}

static int CountInliers(const std::vector<Keypoint> &k1, 
			const std::vector<Keypoint> &k2, 
			std::vector<KeypointMatch> matches,
			double *M, double thresh, std::vector<int> &inliers)
{
    inliers.clear();
    int count = 0;

    for (unsigned int i = 0; i < matches.size(); i++) {
	/* Determine if the ith feature in f1, when transformed by M,
	 * is within RANSACthresh of its match in f2 (if one exists)
	 *
	 * if so, increment count and append i to inliers */

	double p[3];

	p[0] = k1[matches[i].m_idx1].m_x;
	p[1] = k1[matches[i].m_idx1].m_y;
	p[2] = 1.0;

	double q[3];
	matrix_product(3, 3, 3, 1, M, p, q);

	double qx = q[0] / q[2];
	double qy = q[1] / q[2];

	double dx = qx - k2[matches[i].m_idx2].m_x;
	double dy = qy - k2[matches[i].m_idx2].m_y;
	
	double dist = sqrt(dx * dx + dy * dy);
	
	if (dist <= thresh) {
	    count++;
	    inliers.push_back(i);
	}
    }

    return count;
}

static int LeastSquaresFit(const std::vector<Keypoint> &k1, 
			   const std::vector<Keypoint> &k2, 
			   std::vector<KeypointMatch> matches, MotionModel mm,
			   const std::vector<int> &inliers, double *M)
{
    v3_t *r_pts = new v3_t[inliers.size()];
    v3_t *l_pts = new v3_t[inliers.size()];
    double *weight = new double[inliers.size()];

    /* Compute residual */
    double error = 0.0;
    for (int i = 0; i < (int) inliers.size(); i++) {
	int idx1 = matches[inliers[i]].m_idx1;
	int idx2 = matches[inliers[i]].m_idx2;
	
	double r[3], l[3];
	l[0] = k1[idx1].m_x;
	l[1] = k1[idx1].m_y;
	l[2] = 1.0;
		    
	r[0] = k2[idx2].m_x;
	r[1] = k2[idx2].m_y;
	r[2] = 1.0;	

	double rp[3];
	matrix_product(3, 3, 3, 1, M, l, rp);
	
	rp[0] /= rp[2];
	rp[1] /= rp[2];
	
	double dx = rp[0] - r[0];
	double dy = rp[1] - r[1];
	
	error += dx * dx + dy * dy;
    }

    printf("[LeastSquaresFit] Residual error (before) is %0.3e\n", error);    


    for (int i=0; i < (int) inliers.size(); i++) {
	int idx1 = matches[inliers[i]].m_idx1;
	int idx2 = matches[inliers[i]].m_idx2;
	
	Vx(l_pts[i]) = k1[idx1].m_x;
	Vy(l_pts[i]) = k1[idx1].m_y;
	Vz(l_pts[i]) = 1.0;
		    
	Vx(r_pts[i]) = k2[idx2].m_x;
	Vy(r_pts[i]) = k2[idx2].m_y;
	Vz(r_pts[i]) = 1.0;

	weight[i] = 1.0;
    }
    
    switch (mm) {
	case MotionRigid: {
	    double R[9], T[9], Tout[9], scale;
	    align_horn((int) inliers.size(), r_pts, l_pts, R, T, Tout, &scale, weight);
	    memcpy(M, Tout, 9 * sizeof(double));
	    break;
	}
	
	case MotionHomography: {
	    align_homography((int) inliers.size(), r_pts, l_pts, M, 1);
	    break;
	}
    }

    /* Compute residual */
    error = 0.0;
    for (int i = 0; i < (int) inliers.size(); i++) {
	int idx1 = matches[inliers[i]].m_idx1;
	int idx2 = matches[inliers[i]].m_idx2;
	
	double r[3], l[3];
	l[0] = k1[idx1].m_x;
	l[1] = k1[idx1].m_y;
	l[2] = 1.0;
		    
	r[0] = k2[idx2].m_x;
	r[1] = k2[idx2].m_y;
	r[2] = 1.0;	

	double rp[3];
	matrix_product(3, 3, 3, 1, M, l, rp);
	
	rp[0] /= rp[2];
	rp[1] /= rp[2];
	
	double dx = rp[0] - r[0];
	double dy = rp[1] - r[1];
	
	error += dx * dx + dy * dy;
    }
	
    printf("[LeastSquaresFit] Residual error (after) is %0.3e\n", error);    

    delete [] r_pts;
    delete [] l_pts;
    delete [] weight;

    return 0;
}


