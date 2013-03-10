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

/* RelativePose.cpp */
/* Functions for computing the relative pose of two images */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BundlerApp.h"

#include "Decompose.h"
#include "Epipolar.h"
#include "Register.h"
#include "Bundle.h"

#include "sfm.h"
#include "fmatrix.h"
#include "matrix.h"
#include "triangulate.h"

bool BundlerApp::EstimateRelativePose(int i1, int i2, 
                                      camera_params_t &camera1, 
                                      camera_params_t &camera2)
{
    // int num_images = GetNumImages();
    MatchIndex list_idx;

    if (i1 < i2)
        list_idx = GetMatchIndex(i1, i2);
    else
        list_idx = GetMatchIndex(i2, i1);

    // int num_matches = (int) m_match_lists[list_idx].size();
    std::vector<KeypointMatch> &matches = m_matches.GetMatchList(list_idx);
    int num_matches = (int) matches.size();

    double f1 = m_image_data[i1].m_init_focal;
    double f2 = m_image_data[i2].m_init_focal;

    double E[9], F[9];
    std::vector<int> inliers;

    if (!m_optimize_for_fisheye) {
        inliers = 
            EstimateEMatrix(m_image_data[i1].m_keys, m_image_data[i2].m_keys, 
                            matches,
                            4 * m_fmatrix_rounds, // 8 * m_fmatrix_rounds,
                            m_fmatrix_threshold * m_fmatrix_threshold, 
                            f1, f2, E, F);
    } else {
        /* FIXME */
        inliers = 
            EstimateEMatrix(m_image_data[i1].m_keys, m_image_data[i2].m_keys, 
                            matches,
                            4 * m_fmatrix_rounds, // 8 * m_fmatrix_rounds,
                            m_fmatrix_threshold * m_fmatrix_threshold, 
                            f1, f2, E, F);        
    }

    if ((int) inliers.size() == 0)
        return false;

    int num_inliers = (int) inliers.size();
    printf("  Found %d / %d inliers (%0.3f%%)\n", num_inliers, num_matches,
           100.0 * num_inliers / num_matches);


    /* Estimate a homography with the inliers */
    std::vector<KeypointMatch> match_inliers;
    for (int i = 0; i < num_inliers; i++) {
        match_inliers.push_back(matches[inliers[i]]);
    }

    int num_match_inliers = (int) match_inliers.size();

    double H[9];
    std::vector<int> Hinliers = 
        EstimateTransform(m_image_data[i1].m_keys, m_image_data[i2].m_keys,
                          match_inliers, MotionHomography,
                          128 /*m_homography_rounds*/, 
                          6.0 /*m_homography_threshold*/, H);

    printf("  Found %d / %d homography inliers (%0.3f%%)\n",
           (int) Hinliers.size(), num_inliers, 
           100.0 * Hinliers.size() / num_inliers);

    bool initialized = false;
    if ((int) Hinliers.size() > 0) {
        matrix_print(3, 3, H);
        printf("\n");
    
        if ((double) Hinliers.size() / num_inliers >= 0.75 /*0.85*/) {
            KeypointMatch &match0 = matches[Hinliers[0]];
            v2_t p10 = v2_new(m_image_data[i1].m_keys[match0.m_idx1].m_x,
                              m_image_data[i1].m_keys[match0.m_idx1].m_y);
            v2_t p20 = v2_new(m_image_data[i2].m_keys[match0.m_idx2].m_x,
                              m_image_data[i2].m_keys[match0.m_idx2].m_y);

            double R1[9], t1[3], R2[9], t2[3];
            bool success = 
                DecomposeHomography(H, f1, f2, R1, t1, R2, t2, p10, p20);
        
            if (success) {
                printf("[SifterApp::BundleTwoFrame] Using homography "
                       "for initialization\n");

                /* Decide which solution to use */
                double F1h[9], F2h[9];
                ComputeFundamentalMatrix(f1, f2, R1, t1, F1h);
                ComputeFundamentalMatrix(f1, f2, R2, t2, F2h);
                
                double F1hT[9], F2hT[9];
                matrix_transpose(3, 3, F1h, F1hT);
                matrix_transpose(3, 3, F2h, F2hT);

                int num_inliers1 = 0, num_inliers2 = 0;

                for (int i = 0; i < num_match_inliers; i++) {
                    const KeypointMatch &match = match_inliers[i];
                    const Keypoint &k1 = m_image_data[i1].m_keys[match.m_idx1];
                    const Keypoint &k2 = m_image_data[i2].m_keys[match.m_idx2];

                    v3_t rt = v3_new(k1.m_x, k1.m_y, 1.0);
                    v3_t lft = v3_new(k2.m_x, k2.m_y, 1.0);
                    double r1a = fmatrix_compute_residual(F1h, lft, rt);
                    double r1b = fmatrix_compute_residual(F1hT, rt, lft);

                    double r2a = fmatrix_compute_residual(F2h, lft, rt);
                    double r2b = fmatrix_compute_residual(F2hT, rt, lft);

                    if (r1a < m_fmatrix_threshold && r1b < m_fmatrix_threshold)
                        num_inliers1++;

                    if (r2a < m_fmatrix_threshold && r2b < m_fmatrix_threshold)
                        num_inliers2++;
                }

                initialized = true;

                double *R, *t;
                printf("  H1: %d inliers, H2: %d inliers\n", 
                       num_inliers1, num_inliers2);
                if (num_inliers1 > num_inliers2) {
                    R = R1;
                    t = t1;
                } else {
                    R = R2;
                    t = t2;
                }

                memcpy(camera2.R, R, sizeof(double) * 9);
                matrix_transpose_product(3, 3, 3, 1, R, t, camera2.t);
                matrix_scale(3, 1, camera2.t, -1.0, camera2.t);
            }
        }
    }

    if (!initialized) {
        KeypointMatch &match = matches[inliers[0]];
        v2_t p1 = v2_new(m_image_data[i1].m_keys[match.m_idx1].m_x / f1,
                         m_image_data[i1].m_keys[match.m_idx1].m_y / f1);
        v2_t p2 = v2_new(m_image_data[i2].m_keys[match.m_idx2].m_x / f2,
                         m_image_data[i2].m_keys[match.m_idx2].m_y / f2);
    
        double R[9], t[3];
        int success = find_extrinsics_essential(E, p1, p2, R, t);
    
        if (!success) {
            return false;
        }

        memcpy(camera2.R, R, sizeof(double) * 9);

        matrix_transpose_product(3, 3, 3, 1, R, t, camera2.t);
        matrix_scale(3, 1, camera2.t, -1.0, camera2.t);
    }

    return true;
}

bool BundlerApp::EstimateRelativePose2(int i1, int i2, 
                                       camera_params_t &camera1, 
                                       camera_params_t &camera2)
{
    // int num_images = GetNumImages();
    MatchIndex list_idx;

    if (i1 < i2)
        list_idx = GetMatchIndex(i1, i2); // i1 * num_images + i2;
    else
        list_idx = GetMatchIndex(i2, i1); // i2 * num_images + i1;

    std::vector<KeypointMatch> &matches = m_matches.GetMatchList(list_idx);
    // int num_matches = (int) m_match_lists[list_idx].size();
    int num_matches = (int) matches.size();

    // double f1 = m_image_data[i1].m_init_focal;
    // double f2 = m_image_data[i2].m_init_focal;
    double K1[9], K2[9];
    GetIntrinsics(camera1, K1);
    GetIntrinsics(camera2, K2);

    double R0[9], t0[3];
    int num_inliers = 0;

    if (!m_optimize_for_fisheye) {
        num_inliers = 
            EstimatePose5Point(m_image_data[i1].m_keys, 
                               m_image_data[i2].m_keys, 
                               matches,
                               512, /* m_fmatrix_rounds, 8 * m_fmatrix_rounds */
                               0.25 * m_fmatrix_threshold, // 0.003, // 0.004 /*0.001,*/ // /*0.5 **/ m_fmatrix_threshold, 
                               K1, K2, R0, t0);
    } else {
        std::vector<Keypoint> k1 = m_image_data[i1].UndistortKeysCopy();
        std::vector<Keypoint> k2 = m_image_data[i2].UndistortKeysCopy();

        num_inliers = 
            EstimatePose5Point(k1, k2, matches,
                               1024, /*512*/ /* m_fmatrix_rounds, 8 * m_fmatrix_rounds */
                               0.25 * m_fmatrix_threshold, // 0.004, /*0.001,*/ // /*0.5 **/ m_fmatrix_threshold, 
                               K1, K2, R0, t0);
    }
    
    if (num_inliers == 0)
        return false;

    printf("  Found %d / %d inliers (%0.3f%%)\n", num_inliers, num_matches,
           100.0 * num_inliers / num_matches);

    bool initialized = false;
    if (!initialized) {
        memcpy(camera2.R, R0, sizeof(double) * 9);

        matrix_transpose_product(3, 3, 3, 1, R0, t0, camera2.t);
        matrix_scale(3, 1, camera2.t, -1.0, camera2.t);
    }

    return true;
}
