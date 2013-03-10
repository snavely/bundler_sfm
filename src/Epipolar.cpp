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

/* Epipolar.cpp */
/* Routines for computing epipolar geometry */

#include <assert.h>
#include <math.h>
#include <string.h>

#include "Epipolar.h"

#ifndef __DEMO__
#include "5point.h"
#endif

#include "defines.h"
#include "fmatrix.h"
#include "matrix.h"
#include "triangulate.h"
#include "vector.h"

/* Estimate an E-matrix from a given set of point matches */
std::vector<int> EstimateEMatrix(const std::vector<Keypoint> &k1, 
                                 const std::vector<Keypoint> &k2, 
                                 std::vector<KeypointMatch> matches, 
                                 int num_trials, double threshold, 
                                 double f1, double f2, 
                                 double *E, double *F)
{
    int num_keys1 = k1.size();
    int num_keys2 = k2.size();

    std::vector<Keypoint> k1_norm, k2_norm;
    k1_norm.resize(num_keys1);
    k2_norm.resize(num_keys2);
    
    for (int i = 0; i < num_keys1; i++) {
        Keypoint k;
        k.m_x = k1[i].m_x / f1;
        k.m_y = k1[i].m_y / f1;
        k1_norm[i] = k;
    }

    for (int i = 0; i < num_keys2; i++) {
        Keypoint k;
        k.m_x = k2[i].m_x / f2;
        k.m_y = k2[i].m_y / f2;
        k2_norm[i] = k;
    }

    double scale = 0.5 * (f1 + f2);

    std::vector<int> inliers = 
        EstimateFMatrix(k1_norm, k2_norm, matches, num_trials, 
                        threshold / (scale * scale), E, true);
    
    double K1_inv[9] = { 1.0 / f1, 0.0, 0.0, 
                         0.0, 1.0 / f1, 0.0,
                         0.0, 0.0, 1.0 };
    double K2_inv[9] = { 1.0 / f2, 0.0, 0.0, 
                         0.0, 1.0 / f2, 0.0,
                         0.0, 0.0, 1.0 };
    
    double tmp[9];
    matrix_product(3, 3, 3, 3, K1_inv, E, tmp);
    matrix_product(3, 3, 3, 3, K2_inv, tmp, F);
    
    return inliers;
}

#ifndef __DEMO__
/* Estimate relative pose from a given set of point matches */
int EstimatePose5Point(const std::vector<Keypoint> &k1, 
                       const std::vector<Keypoint> &k2, 
                       std::vector<KeypointMatch> matches, 
                       int num_trials, double threshold, 
                       double *K1, double *K2, 
                       double *R, double *t)
{
    int num_pts = (int) matches.size();

    v2_t *k1_pts = new v2_t[num_pts];
    v2_t *k2_pts = new v2_t[num_pts];

    for (int i = 0; i < num_pts; i++) {
	int idx1 = matches[i].m_idx1;
	int idx2 = matches[i].m_idx2;

	k1_pts[i] = v2_new(k1[idx1].m_x, k1[idx1].m_y);
	k2_pts[i] = v2_new(k2[idx2].m_x, k2[idx2].m_y);
    }

    int num_inliers = compute_pose_ransac(num_pts, k1_pts, k2_pts, 
                                          K1, K2, threshold, num_trials, R, t);

    delete [] k1_pts;
    delete [] k2_pts;

    return num_inliers;
}
#endif

/* Estimate an F-matrix from a given set of point matches */
std::vector<int> EstimateFMatrix(const std::vector<Keypoint> &k1, 
				 const std::vector<Keypoint> &k2, 
				 std::vector<KeypointMatch> matches, 
				 int num_trials, double threshold, 
				 double *F, bool essential)
{
    int num_pts = (int) matches.size();

    /* num_pts should be greater than a threshold */
    if (num_pts < 20) {
        std::vector<int> inliers;
        return inliers;
    }

    v3_t *k1_pts = new v3_t[num_pts];
    v3_t *k2_pts = new v3_t[num_pts];

    v3_t *k1_pts_in = new v3_t[num_pts];
    v3_t *k2_pts_in = new v3_t[num_pts];

    for (int i = 0; i < num_pts; i++) {
        int idx1 = matches[i].m_idx1;
        int idx2 = matches[i].m_idx2;

        assert(idx1 < (int) k1.size());
        assert(idx2 < (int) k2.size());

        k1_pts[i] = v3_new(k1[idx1].m_x, k1[idx1].m_y, 1.0);
        k2_pts[i] = v3_new(k2[idx2].m_x, k2[idx2].m_y, 1.0);
    }

    estimate_fmatrix_ransac_matches(num_pts, k2_pts, k1_pts, 
        num_trials, threshold, 0.95,
        (essential ? 1 : 0), F);

    /* Find the inliers */
    std::vector<int> inliers;

    for (int i = 0; i < num_pts; i++) {
        double dist = fmatrix_compute_residual(F, k2_pts[i], k1_pts[i]);
        if (dist < threshold) {
            inliers.push_back(i);
        }
    }

    /* Re-estimate using inliers */
    int num_inliers = (int) inliers.size();

    for (int i = 0; i < num_inliers; i++) {
        k1_pts_in[i] = k1_pts[inliers[i]]; // v3_new(k1[idx1]->m_x, k1[idx1]->m_y, 1.0);
        k2_pts_in[i] = k2_pts[inliers[i]]; // v3_new(k2[idx2]->m_x, k2[idx2]->m_y, 1.0);
    }

    // printf("[1] num_inliers = %d\n", num_inliers);

#if 0
    double F0[9];
    double e1[3], e2[3];
    estimate_fmatrix_linear(num_inliers, k2_pts_in, k1_pts_in, F0, e1, e2);

    inliers.clear();
    for (int i = 0; i < num_pts; i++) {
        double dist = fmatrix_compute_residual(F0, k2_pts[i], k1_pts[i]);
        if (dist < threshold) {
            inliers.push_back(i);
        }
    }
    num_inliers = inliers.size();
    // printf("[2] num_inliers = %d\n", num_inliers);

    // matrix_print(3, 3, F0);
#else
    double F0[9];
    memcpy(F0, F, sizeof(double) * 9);
#endif

    if (!essential) {
        /* Refine using NLLS */
        for (int i = 0; i < num_inliers; i++) {
            k1_pts_in[i] = k1_pts[inliers[i]];
            k2_pts_in[i] = k2_pts[inliers[i]];
        }

        refine_fmatrix_nonlinear_matches(num_inliers, k2_pts_in, k1_pts_in, 
            F0, F);
    } else {
        memcpy(F, F0, sizeof(double) * 9);
    }

#if 0
    if (essential) {
        /* Compute the SVD of F */
        double U[9], S[3], VT[9];
        dgesvd_driver(3, 3, F, U, S, VT);
        double E0[9] = { 1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 0.0 };

        double tmp[9];
        matrix_product(3, 3, 3, 3, U, E0, tmp);
        matrix_product(3, 3, 3, 3, tmp, VT, F);
    }
#endif

    inliers.clear();
    for (int i = 0; i < num_pts; i++) {
        double dist = fmatrix_compute_residual(F, k2_pts[i], k1_pts[i]);
        if (dist < threshold) {
            inliers.push_back(i);
        }
    }
    num_inliers = (int) inliers.size();

    delete [] k1_pts;
    delete [] k2_pts;
    delete [] k1_pts_in;
    delete [] k2_pts_in;

    return inliers;
}

std::vector<int> EstimateFMatrix(const std::vector<KeypointWithDesc> &k1, 
				 const std::vector<KeypointWithDesc> &k2, 
				 std::vector<KeypointMatch> matches, 
				 int num_trials, double threshold, 
				 double *F, bool essential)
{
    int num_keys1 = (int) k1.size();
    int num_keys2 = (int) k2.size();
    
    std::vector<Keypoint> k1_prime, k2_prime;

    k1_prime.resize(num_keys1);
    k2_prime.resize(num_keys2);

    for (int i = 0; i < num_keys1; i++) {
        Keypoint k(k1[i].m_x, k1[i].m_y);
        k1_prime[i] = k;
    }
    
    for (int i = 0; i < num_keys2; i++) {
        Keypoint k(k2[i].m_x, k2[i].m_y);
        k2_prime[i] = k;
    }

    return EstimateFMatrix(k1_prime, k2_prime, matches, num_trials, 
                           threshold, F, essential);
}
