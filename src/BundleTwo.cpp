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

/* BundleTwo.cpp */
/* Two-frame bundle adjustment */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>

#include "BundlerApp.h"
#include "Bundle.h"
#include "BundleAdd.h"
#include "BundleUtil.h"
#include "Decompose.h"
#include "TwoFrameModel.h"

#include "sfm.h"

#include "defines.h"
#include "filter.h"
#include "horn.h"
#include "matrix.h"
#include "qsort.h"
#include "triangulate.h"
#include "util.h"
#include "vector.h"

#ifndef WIN32
#include <ext/hash_map>
#else
#include <hash_map>
#endif

/* Stolen from TwoFrameModel.cpp */
static void WriteVector(FILE *f, int n, const double *v)
{
    for (int i = 0; i < n; i++) {
        fprintf(f, "%0.16e ", v[i]);
    }
    fprintf(f, "\n");
}

static void WriteCamera(FILE *f, const camera_params_t &camera)
{
    WriteVector(f, 9, camera.R);
    WriteVector(f, 3, camera.t);
    fprintf(f, "%0.16e\n", camera.f);
}

void TwoFrameModel::WriteWithProjections(FILE *f, 
                                         const std::vector<TrackData> &tracks,
                                         int i1, int i2,
                                         const ImageData &img1,
                                         const ImageData &img2) const
{
    std::vector<v2_t> proj1, proj2;
    
    proj1.resize(m_num_points);
    proj2.resize(m_num_points);

    fprintf(f, "%d\n", m_num_points);

    fprintf(f, "%0.9f\n", m_angle);

    for (int i = 0; i < m_num_points; i++) {
        int tr;
            
        if (m_tracks == NULL)
            tr = -1;
        else
            tr = m_tracks[i];

        fprintf(f, "%d %0.16e %0.16e %0.16e\n", tr, 
                Vx(m_points[i]), Vy(m_points[i]), Vz(m_points[i]));
        
        int k1 = -1, k2 = -1;
        if (tr != -1) {
            const ImageKeyVector &v = tracks[tr].m_views;
            
            int num_views = (int) v.size();
            for (int j = 0; j < num_views; j++) {
                if (v[j].first == i1)
                    k1 = v[j].second;
                else if (v[j].first == i2)
                    k2 = v[j].second;
            }
        } else {
            k1 = m_keys1[i];
            k2 = m_keys2[i];
        }        

        assert(k1 >= 0 && k2 >= 0);
            
        double x1 = img1.m_keys[k1].m_x;
        double y1 = img1.m_keys[k1].m_y;
        double x2 = img2.m_keys[k2].m_x;
        double y2 = img2.m_keys[k2].m_y;

        if (img1.m_fisheye)
            img1.UndistortPoint(x1, y1, x1, y1);
        if (img2.m_fisheye)
            img2.UndistortPoint(x2, y2, x2, y2);

        fprintf(f, "%0.9e %0.9e %0.9e %0.9e\n", x1, y1, x2, y2);

        proj1[i] = v2_new(x1, y1);
        proj2[i] = v2_new(x2, y2);
    }
    
    WriteCamera(f, m_camera0);
    WriteCamera(f, m_camera1);
    
    WriteVector(f, 9, m_C0);
    WriteVector(f, 9, m_C1);

#if 0
    /* Compute the translation using DLT */
    double pts1[30] = {
        0.0210,    0.3174,    1.0000,
        -0.3350,    0.0522,    1.0000,
        -0.1398,    0.3074,    1.0000,
        -0.8325,   -0.1650,    1.0000,
        -1.2247,    0.4144,    1.0000,
        -0.4274,    0.6393,    1.0000,
        -0.5732,   -0.0541,    1.0000,
        -0.2511,   -0.0100,    1.0000,
        -0.0131,    0.4415,    1.0000,
        0.1108,    0.1131,    1.0000
    };

    double pts2[30] = {
        0.1993,   -0.1316,    1.0000,
        15.3933,    1.5673,    1.0000,
        0.4956,   -0.0401,    1.0000,
        0.8608,    0.3763,    1.0000,
        0.5506,    0.1301,    1.0000,
        0.5927,   -0.3441,    1.0000,
        2.5845,    0.5888,    1.0000,
        -40.0390,  -13.0980,    1.0000,
        0.2053,   -0.0880,    1.0000,
        0.0701,    0.5576,    1.0000,
    };

    double R1[9] = {
        0.8440,    0.1616,    0.5113,
        0.2310,    0.7509,   -0.6187,
        -0.4840,    0.6403,    0.5965,
    };

    double t1[3] = {
        -0.3775, -0.2959, -1.4751
    };
    
    double R2[9] = {
        0.8189,    0.4004,    0.4113,
        -0.0417,    0.7561,   -0.6532,
        -0.5725,    0.5177,    0.6358,
    };

    double t2[3] = {
        -0.2340, 0.1184, 0.3148
    };


    double R12[9];
    matrix_transpose_product2(3, 3, 3, 3, R2, R1, R12);

    double A[30];

    for (int i = 0; i < 10; i++) {
        double *p1 = pts1 + 3 * i;
        double *p2 = pts2 + 3 * i;

        double p1b[3];
        matrix_product(3, 3, 3, 1, R12, p1, p1b);
        p1[0] = p1b[0] / p1b[2];
        p1[1] = p1b[1] / p1b[2];
        
        double *r = A + i * 3;
        
        r[0] = -p1[1] + p2[1];
        r[1] = p1[0] - p2[0];
        r[2] = -p1[0] * p2[1] + p1[1] * p2[0];
    }
    
    double x[3];
    matrix_minimum_unit_norm_solution(10, 3, A, x);

    double c1[3], c2[3];
    matrix_transpose_product(3, 3, 3, 1, R1, t1, c1);
    matrix_scale(3, 1, c1, -1.0, c1);
    matrix_transpose_product(3, 3, 3, 1, R2, t2, c2);
    matrix_scale(3, 1, c2, -1.0, c2);

    double c[3];
    matrix_diff(3, 1, 3, 1, c2, c1, c);

    double t[3];
    matrix_product(3, 3, 3, 1, (double *) R2, c, t);
    matrix_scale(3, 1, t, -1.0, t);

    double norm = matrix_norm(3, 1, t);
    matrix_scale(3, 1, t, 1.0 / norm, t);

    double dot;
    matrix_product(1, 3, 3, 1, x, t, &dot);
    printf("%0.3f %0.3f %0.3f -> %0.3f %0.3f %0.3f [ %0.3f ]\n",
           x[0], x[1], x[2], t[0], t[1], t[2], dot);
#endif

#if 0
    int num_eqns = m_num_points;
    int num_vars = 3;
    double *A = new double[num_eqns * num_vars];

    for (int i = 0; i < m_num_points; i++) {
        v2_t p1 = v2_scale(-1.0 / m_camera0.f, proj1[i]);
        v2_t p2 = v2_scale(-1.0 / m_camera1.f, proj2[i]);
        
        double p1a[3] = { Vx(p1), Vy(p1), 1.0 };
        double p1b[3];
        matrix_product(3, 3, 3, 1, m_camera1.R, p1a, p1b);
        // p1 = v2_new(p1b[0] / p1b[2], p1b[1] / p1b[2]);
        p1 = v2_new(p1b[0], p1b[1]);
        
        double *r = A + i * num_vars;
        
        r[0] = -Vy(p1) + Vy(p2) * p1b[2];
        r[1] = Vx(p1) - Vx(p2) * p1b[2];
        r[2] = -Vx(p1) * Vy(p2) + Vy(p1) * Vx(p2);
    }
    
    double x[3];
    matrix_minimum_unit_norm_solution(num_eqns, num_vars, A, x);

    double c[3];
    matrix_diff(3, 1, 3, 1, (double *) m_camera1.t, (double *) m_camera0.t, c);

    double t[3];
    matrix_product(3, 3, 3, 1, (double *) m_camera1.R, c, t);
    matrix_scale(3, 1, t, -1.0, t);

    double norm = matrix_norm(3, 1, t);
    matrix_scale(3, 1, t, 1.0 / norm, t);

    double dot;
    matrix_product(1, 3, 3, 1, x, t, &dot);
    printf("%0.3f %0.3f %0.3f -> %0.3f %0.3f %0.3f [ %0.3f ]\n",
           x[0], x[1], x[2], t[0], t[1], t[2], dot);

    delete [] A;
#endif
}

#if 0
void WriteModelsProjections(ModelMap &models, int num_images, 
                            const std::vector<TrackData> &track_data,
                            std::vector<ImageData> &image_data, char *out_file)
{
    bool preload_keys = true;
    
    if (num_images > 30000)
        preload_keys = false;

    if (preload_keys) {
        for (int i = 0; i < num_images; i++) {
            image_data[i].LoadKeys(false);
        }
    }

    FILE *f = fopen(out_file, "w");
    if (f == NULL) {
        printf("[WriteModelsProjections] "
               "Error opening file %s for reading\n", out_file);
    } else {
        fprintf(f, "%d\n", num_images);

        // FIXME LOOP
        for (int i = 0; i < num_images; i++) {
            ImageData &img_i = image_data[i];
            
            if (!preload_keys)
                image_data[i].LoadKeys(false);

            for (int j = i+1; j < num_images; j++) {
                if (models.Contains(GetMatchIndex(i, j))) {
                    if (!preload_keys)
                        image_data[j].LoadKeys(false);

                    ImageData &img_j = image_data[j];

                    fprintf(f, "%d %d\n", i, j);
                    models.GetModel(GetMatchIndex(i,j)).
                        WriteWithProjections(f, track_data, i, j, 
                                             img_i, img_j);

                    if (!preload_keys)
                        image_data[j].UnloadKeys();
                }
            }

            if (!preload_keys)
                image_data[i].UnloadKeys();
        }
        
        fclose(f);
    }

    if (!preload_keys) {
        for (int i = 0; i < num_images; i++) {
            image_data[i].UnloadKeys();
        }
    }
}
#endif

#if 0
static void FixPEdges(int num_images, ModelMap &models, PEdgeMap &p_edges)
{
    const double MIN_ANGLE = 1.5;
    
    for (int i = 0; i < num_images; i++) {
        for (int j = i+1; j < num_images; j++) {
            int idx = i * num_images + j;
            if (models.find(idx) != models.end()) {
                if (models[idx].m_angle < MIN_ANGLE) {
                    printf("[FixPEdges] Angle (%d,%d) too small [%0.3f]\n",
                           i, j, models[idx].m_angle);

                    if (models[idx].m_num_points >= 64) {
                        printf("[FixPEdges] Replacing with p-edge\n");
                        p_edges[idx] = true;
                    }

                    models.erase(idx);
                }
            }
        }
    }
}
#endif

static void FixScaffoldEdges(int num_images, ModelMap &models) 
{
    const double MIN_ANGLE = 1.5;
    
    // FIXME LOOP ?
    for (int i = 0; i < num_images; i++) {
        for (int j = i+1; j < num_images; j++) {
            MatchIndex idx = GetMatchIndex(i,j); // i * num_images + j;
            if (models.Contains(idx)) {
                TwoFrameModel &m = models.GetModel(idx);

                if (m.m_angle < MIN_ANGLE) {
                    printf("[FixScaffoldEdges] "
                           "Angle (%d,%d) too small [%0.3f]\n",
                           i, j, m.m_angle);

#if 0
                    if (m.m_num_points >= 64) {
                        printf("[FixScaffoldEdges] "
                               "Replacing with connector edge\n");
                        m.TurnOffScaffold();
                        // p_edges[idx] = true;
                    } else {
                        models.erase(idx);
                    }
#else
                    models.RemoveModel(idx);
#endif
                }
            }
        }
    }
}

#define MATCH_THRESHOLD 28   // 16

static void PermuteMatches(std::vector<KeypointMatch> &matches) 
{
    int num_matches = (int) matches.size();
    int *perm = new int[num_matches];
    generate_permutation(num_matches, perm);
    std::vector<KeypointMatch> matches_new;
    matches_new.resize(num_matches);

    for (int i = 0; i < num_matches; i++) {
        matches_new[i] = matches[perm[i]];
    }

    matches = matches_new;

    delete [] perm;
}

static void ClearKeys(ImageData &data)
{
    /* Clear keys */
    std::vector<Keypoint>::iterator iter;
    for (iter = data.m_keys.begin(); iter != data.m_keys.end(); iter++) {
        iter->m_extra = -1;
    }

    for (iter = data.m_keys.begin(); iter != data.m_keys.end(); iter++) {
        iter->m_extra = -1;
    }    
}

double BundlerApp::RunSFMNecker(int i1, int i2, 
                                camera_params_t *cameras, 
                                int num_points, v3_t *points, v3_t *colors,
                                std::vector<ImageKeyVector> &pt_views,
                                camera_params_t *cameras_new,
                                v3_t *points_new,
                                double threshold)
{
    memcpy(points_new, points, sizeof(v3_t) * num_points);
    memcpy(cameras_new, cameras, sizeof(camera_params_t) * 2);

    /* Swap cameras */
    camera_params_t tmp = cameras_new[0];
    memcpy(cameras_new[0].R, cameras_new[1].R, sizeof(double) * 9);
    memcpy(cameras_new[0].t, cameras_new[1].t, sizeof(double) * 3);

    memcpy(cameras_new[1].R, tmp.R, sizeof(double) * 9);
    memcpy(cameras_new[1].t, tmp.t, sizeof(double) * 3);	

    for (int i = 0; i < num_points; i++) {
        if (pt_views[i].size() == 0)
            continue;

        int k1 = pt_views[i][0].second;
        int k2 = pt_views[i][1].second;

        double proj1[2] = { GetKey(i1,k1).m_x, GetKey(i1,k1).m_y };
        double proj2[2] = { GetKey(i2,k2).m_x, GetKey(i2,k2).m_y };

        if (m_optimize_for_fisheye) {
            double x1 = proj1[0];
            double y1 = proj1[1];
                    
            double x2 = proj2[0];
            double y2 = proj2[1];
                    
            m_image_data[i1].UndistortPoint(x1, y1, proj1[0], proj1[1]);
            m_image_data[i2].UndistortPoint(x2, y2, proj2[0], proj2[1]);
        }
	    
        v2_t p = v2_new(proj1[0], proj1[1]);
        v2_t q = v2_new(proj2[0], proj2[1]);

        /* Triangulate the point */
        bool in_front = true;
        double angle = 0.0;
        double error = 0.0;
        points_new[i] = Triangulate(p, q, cameras_new[0], cameras_new[1], 
                                    error, in_front, angle, 
                                    true);
    }

    int added_order[2] = { i1, i2 };
    double error1;
    error1 = RunSFM_SBA(num_points, 2, 0, false,
                        cameras_new, points_new, added_order, 
                        colors, pt_views,
                        threshold, NULL, NULL, NULL, NULL, true);

    return error1;
}



bool BundlerApp::BundleTwoFrame(int i1, int i2, TwoFrameModel *model,
                                double &angle_out, int &num_pts_out, 
                                bool bundle_from_tracks) 
{
    const double TERM_THRESH = 1.0e-12;

    // assert(!m_estimate_distortion && !m_fixed_focal_length);
    assert(!m_fixed_focal_length);

    if (!m_image_data[i1].m_has_init_focal || 
        !m_image_data[i2].m_has_init_focal) {

        printf("[BundleTwoFrame] "
               "Error: two frames must have focal length estimates\n");
        return false;
    }

    camera_params_t cameras[2];

    /* Load the keys for the images */
    if (!m_image_data[i1].m_keys_loaded)
        m_image_data[i1].LoadKeys(false, !m_optimize_for_fisheye);
    if (!m_image_data[i2].m_keys_loaded)
        m_image_data[i2].LoadKeys(false, !m_optimize_for_fisheye);

    // #define USE_COLORS
#ifdef USE_COLORS
    m_image_data[i1].ReadKeyColors();
    m_image_data[i2].ReadKeyColors();
#endif /* USE_COLORS */

    if (bundle_from_tracks)
        SetMatchesFromTracks(i1, i2);

    MatchIndex list_idx;

    if (i1 < i2) 
        list_idx = GetMatchIndex(i1, i2); // i1 * num_images + i2;
    else
        list_idx = GetMatchIndex(i2, i1); // i2 * num_images + i1;

    unsigned int num_matches = m_matches.GetNumMatches(list_idx);

    if (num_matches > m_image_data[i1].m_keys.size() - 5 && 
        num_matches > m_image_data[i2].m_keys.size() - 5) {

        printf("[BundleTwoFrame] Identical images!\n");
        angle_out = 0.0;
        num_pts_out = (int) num_matches;
        return true;
    }

    if (m_keypoint_border_width > 0) {
        RemoveMatchesNearBorder(i1, i2, m_keypoint_border_width);

        // int num_matches = (int) m_match_lists[list_idx].size();
        int num_matches = m_matches.GetNumMatches(list_idx);

        if (num_matches < MATCH_THRESHOLD) {
            printf("[BundleTwoFrame] Removed too many matches\n");
            return false;
        }
    }

    InitializeCameraParams(m_image_data[i1], cameras[0]);
    InitializeCameraParams(m_image_data[i2], cameras[1]);

    std::vector<int> tracks;
    std::vector<ImageKeyVector> pt_views;
    // int num_init_cams = 0;

    /* Clear keys */
    std::vector<Keypoint>::iterator iter;
    for (iter = m_image_data[i1].m_keys.begin(); 
         iter != m_image_data[i1].m_keys.end(); 
         iter++) {
        
        iter->m_extra = -1;
    }

    for (iter = m_image_data[i2].m_keys.begin(); 
         iter != m_image_data[i2].m_keys.end(); 
         iter++) {
        
        iter->m_extra = -1;
    }

    /* Put first camera at origin */
    cameras[0].R[0] = 1.0;  cameras[0].R[1] = 0.0;  cameras[0].R[2] = 0.0;
    cameras[0].R[3] = 0.0;  cameras[0].R[4] = 1.0;  cameras[0].R[5] = 0.0;
    cameras[0].R[6] = 0.0;  cameras[0].R[7] = 0.0;  cameras[0].R[8] = 1.0;

    /* Initialize the positions of the cameras (using constraints,
     * if provided) */
    if (m_image_data[i1].m_camera.m_constrained[0])
	cameras[0].t[0] = m_image_data[i1].m_camera.m_constraints[0];
    else
	cameras[0].t[0] = 0.0;

    if (m_image_data[i1].m_camera.m_constrained[1])
	cameras[0].t[1] = m_image_data[i1].m_camera.m_constraints[1];
    else
	cameras[0].t[1] = 0.0;

    if (m_image_data[i1].m_camera.m_constrained[2])
	cameras[0].t[2] = m_image_data[i1].m_camera.m_constraints[2];
    else
	cameras[0].t[2] = 0.0;

    if (m_image_data[i1].m_has_init_focal)
	cameras[0].f = m_image_data[i1].m_init_focal;
    else 
	cameras[0].f = m_init_focal_length;

    if (m_image_data[i2].m_has_init_focal)
        cameras[1].f = m_image_data[i2].m_init_focal;
    else
        cameras[1].f = m_init_focal_length;

    SetCameraConstraints(i1, cameras + 0);
    SetCameraConstraints(i2, cameras + 1);

    if (m_constrain_focal) {
	SetFocalConstraint(m_image_data[i1], cameras + 0);
	SetFocalConstraint(m_image_data[i2], cameras + 1);
    }

    printf("[Sifter::BundleTwoFrame] Estimating relative pose...\n");
    fflush(stdout);
    bool success = EstimateRelativePose2(i1, i2, cameras[0], cameras[1]);
    fflush(stdout);

    if (!success) 
        return false;

    // unsigned int num_matches = (int) m_match_lists[list_idx].size();
    // unsigned int num_matches = m_matches.GetNumMatches(list_idx);

    /* **** Set up the initial 3D points **** */
    printf("[BundleTwoFrame] Adding initial matches...\n");
    fflush(stdout);

    int pt_count = 0;

    v3_t *points = new v3_t[num_matches];
#ifdef USE_COLORS
    v3_t *colors = new v3_t[num_matches];
#endif

    double angle_sum = 0.0;
    int num_in_back = 0, num_skipped = 0;

    // std::vector<KeypointMatch> matches = m_match_lists[list_idx];
    std::vector<KeypointMatch> matches = m_matches.GetMatchList(list_idx);

    PermuteMatches(matches);

    // #define OUTPUT_POINT_STATUS
#ifdef OUTPUT_POINT_STATUS
    m_image_data[i1].LoadImage();
    m_image_data[i2].LoadImage();
    img_t *pt_img1 = img_scale(m_image_data[i1].m_img, 2);
    img_t *pt_img2 = img_scale(m_image_data[i2].m_img, 2);
    m_image_data[i1].UnloadImage();
    m_image_data[i2].UnloadImage();
#endif

    for (unsigned int i = 0; i < num_matches; i++) {
	int key_idx1 = matches[i].m_idx1;
	int key_idx2 = matches[i].m_idx2;

#if 0
	printf("  Adding match %d ==> %d [%d]\n", 
	       key_idx1, key_idx2, pt_count);
#endif

        /* Triangulate the point */
        double xp1, yp1, xp2, yp2;
        double x_proj1, y_proj1, x_proj2, y_proj2;

        xp1 = x_proj1 = GetKey(i1,key_idx1).m_x;
        yp1 = y_proj1 = GetKey(i1,key_idx1).m_y;
        xp2 = x_proj2 = GetKey(i2,key_idx2).m_x;
        yp2 = y_proj2 = GetKey(i2,key_idx2).m_y;
        
        if (m_optimize_for_fisheye) {
            m_image_data[i1].UndistortPoint(x_proj1, y_proj1, 
                                            x_proj1, y_proj1);
            m_image_data[i2].UndistortPoint(x_proj2, y_proj2, 
                                            x_proj2, y_proj2);
        }

        double error;
            
        v2_t p = v2_new(x_proj1, y_proj1);
        v2_t q = v2_new(x_proj2, y_proj2);

        bool in_front = true;
        double angle = 0.0;
        points[pt_count] = Triangulate(p, q, cameras[0], cameras[1], 
                                       error, in_front, angle, 
                                       true);

        
        if (m_optimize_for_fisheye) {
            /* Project the point */
            double tmp[3], tmp2[3];
            matrix_diff(3, 1, 3, 1, points[pt_count].p, cameras[0].t, tmp);
            matrix_product331(cameras[0].R, tmp, tmp2);
            double px1 = -cameras[0].f * tmp2[0] / tmp2[2];
            double py1 = -cameras[0].f * tmp2[1] / tmp2[2];

            matrix_diff(3, 1, 3, 1, points[pt_count].p, cameras[1].t, tmp);
            matrix_product331(cameras[1].R, tmp, tmp2);
            double px2 = -cameras[1].f * tmp2[0] / tmp2[2];
            double py2 = -cameras[1].f * tmp2[1] / tmp2[2];

            m_image_data[i1].DistortPoint(px1, py1, px1, py1);
            m_image_data[i2].DistortPoint(px2, py2, px2, py2);

            double dx1 = px1 - xp1;
            double dy1 = py1 - yp1;
            double dx2 = px2 - xp2;
            double dy2 = py2 - yp2;
            
            error = 0.5 * (sqrt(dx1 * dx1 + dy1 * dy1) +
                           sqrt(dx2 * dx2 + dy2 * dy2));
        }

#if 0
        printf(" tri.error[%d] = %0.3f, %0.3f, %d\n", 
               i, error, RAD2DEG(angle), in_front ? 1 : 0);
#endif

#ifdef OUTPUT_POINT_STATUS
        double x_img1 = 0.5 * (xp1 + 0.5 * m_image_data[i1].GetWidth());
        double y_img1 = 0.5 * (yp1 + 0.5 * m_image_data[i1].GetHeight());
        double x_img2 = 0.5 * (xp2 + 0.5 * m_image_data[i2].GetWidth());
        double y_img2 = 0.5 * (yp2 + 0.5 * m_image_data[i2].GetHeight());
#endif

        if (error > 10.0) {
#ifdef OUTPUT_POINT_STATUS
            img_draw_pt(pt_img1, iround(x_img1), iround(y_img1), 4,
                        0xff, 0x0, 0x0);
            img_draw_pt(pt_img2, iround(x_img2), iround(y_img2), 4,
                        0xff, 0x0, 0x0);
#endif

            // printf(" skipping point\n");
            num_skipped++;
            continue;
        }

        if (!in_front) {
#ifdef OUTPUT_POINT_STATUS
            img_draw_pt(pt_img1, iround(x_img1), iround(y_img1), 4,
                        0x0, 0x0, 0xff);
            img_draw_pt(pt_img2, iround(x_img2), iround(y_img2), 4,
                        0x0, 0x0, 0xff);
#endif

            num_in_back++;
            continue;
        }

#ifdef OUTPUT_POINT_STATUS
        img_draw_pt(pt_img1, iround(x_img1), iround(y_img1), 4,
                    0x0, 0xff, 0x0);
        img_draw_pt(pt_img2, iround(x_img2), iround(y_img2), 4,
                    0x0, 0xff, 0x0);
#endif

        angle_sum += angle;

#ifdef USE_COLORS
        /* Get the color of the point */
        unsigned char r = GetKey(i1,key_idx1).m_r;
        unsigned char g = GetKey(i1,key_idx1).m_g;
        unsigned char b = GetKey(i1,key_idx1).m_b;
        colors[pt_count] = v3_new((double) r, (double) g, (double) b);
#endif     
   
        GetKey(i1,key_idx1).m_extra = pt_count;
        GetKey(i2,key_idx2).m_extra = pt_count;

        if (bundle_from_tracks) {
            int track_idx = GetKey(i1,key_idx1).m_track;
            m_track_data[track_idx].m_extra = pt_count;
            tracks.push_back(track_idx);
        } 

        ImageKeyVector views;
        views.push_back(ImageKey(0, key_idx1));
        views.push_back(ImageKey(1, key_idx2));
        pt_views.push_back(views);
        
        pt_count++;
    } /* end loop through matches */

#ifdef OUTPUT_POINT_STATUS
    char ptbuf[256]; 
    sprintf(ptbuf, "pt%03d-%03d.bmp", i1, i2);
    img_write_bmp_file(pt_img1, ptbuf);
    sprintf(ptbuf, "pt%03d-%03d.bmp", i2, i1);
    img_write_bmp_file(pt_img2, ptbuf);

    img_free(pt_img1);
    img_free(pt_img2);
#endif

    double angle_avg = angle_sum / pt_count;
    printf("  Average angle: %0.3f\n", RAD2DEG(angle_avg));
    printf("  In-back pct  : %d / %d (%0.3f%%)\n", num_in_back, num_matches,
           100.0 * num_in_back / num_matches);
    printf("  Skipped pct  : %d / %d (%0.3f%%)\n", num_skipped, num_matches,
           100.0 * num_skipped / num_matches);
    fflush(stdout);

    angle_out = RAD2DEG(angle_avg);
    num_pts_out = pt_count;

    if (RAD2DEG(angle_avg) < 0.5 /*1.5*/) { /* Thresh on triangulation angle */
        printf("[BundleTwoFrame] Average tri.angle too small, "
               "aborting!\n");
        
        ClearKeys(m_image_data[i1]);
        ClearKeys(m_image_data[i2]);

        return false;
    }

    if ((double) num_in_back / num_matches > 0.30 /*0.035*/) {
        /* Too many points are in back of the cameras, abort */
        printf("[BundleTwoFrame] Too many points [%d / %d / %0.3f] "
               "in back of the cameras, aborting!\n", 
               num_in_back, num_matches, 100.0 * num_in_back / num_matches);

        ClearKeys(m_image_data[i1]);
        ClearKeys(m_image_data[i2]);
        return false;
    }

    if (pt_count < 20) {   
        ClearKeys(m_image_data[i1]);
        ClearKeys(m_image_data[i2]);
        return false;
    }

    // m_match_lists[list_idx].clear();
    m_matches.ClearMatch(list_idx);
    matches.clear();

    /* Add constraints to camera 0 to fix position and rotation */
    cameras[0].constrained[0] = true;
    cameras[0].constrained[1] = true;
    cameras[0].constrained[2] = true;
    cameras[0].constrained[3] = true;
    cameras[0].constrained[4] = true;
    cameras[0].constrained[5] = true;

    cameras[0].constraints[0] = 0.0;
    cameras[0].constraints[1] = 0.0;
    cameras[0].constraints[2] = 0.0;
    cameras[0].constraints[3] = 0.0;
    cameras[0].constraints[4] = 0.0;
    cameras[0].constraints[5] = 0.0;

    cameras[0].weights[0] = 1.0e6;
    cameras[0].weights[1] = 1.0e6;
    cameras[0].weights[2] = 1.0e6;
    cameras[0].weights[3] = 1.0e6;
    cameras[0].weights[4] = 1.0e6;
    cameras[0].weights[5] = 1.0e6;

    /* Constraint radial distortion parameters */
    cameras[0].constrained[7] = true;
    cameras[0].constrained[8] = true;
    cameras[0].constraints[7] = 0.0;
    cameras[0].constraints[8] = 0.0;
    cameras[0].weights[7] = 1.0e2;
    cameras[0].weights[8] = 1.0e2;

    cameras[1].constrained[7] = true;
    cameras[1].constrained[8] = true;
    cameras[1].constraints[7] = 0.0;
    cameras[1].constraints[8] = 0.0;
    cameras[1].weights[7] = 1.0e2;
    cameras[1].weights[8] = 1.0e2;

    int added_order[2] = { i1, i2 };

    /* ********** Bundle adjust the two views ********** */
    // std::vector<ImageKeyVector> pt_views_new = pt_views;

#if 1
#ifdef USE_COLORS
    double error0 = RunSFM_SBA(pt_count, 2, 0, false,
                               cameras, points, added_order, 
                               colors, pt_views,
                               TERM_THRESH, NULL, NULL, NULL, NULL, true);
#else
    double error0 = RunSFM_SBA(pt_count, 2, 0, false,
                               cameras, points, added_order, NULL, pt_views,
                               TERM_THRESH, NULL, NULL, NULL, NULL, true);
#endif
#endif

#if 0
    camera_params_t cameras_new[2];
    v3_t *points_new = new v3_t[pt_count];
    double error1 = RunSFMNecker(i1, i2, cameras, pt_count, points, colors,
                                 pt_views_new, cameras_new, points_new, 
                                 TERM_THRESH);
    
    double error = MIN(error0, error1);

    if (error1 < error0) {
        printf("  Switching to reflected solution (%0.3f < %0.3f)\n",
               error1, error0);
        memcpy(points, points_new, pt_count * sizeof(v3_t));
        memcpy(cameras, cameras_new, 2 * sizeof(camera_params_t));
        pt_views = pt_views_new;
    } else {
        printf("  Keeping initial solution (%0.3f < %0.3f)\n",
               error0, error1);
    }

    delete [] points_new;
    pt_views_new.clear();
#endif

    if (error0 > 0.5) {
        printf("[BundleTwoFrame] "
               "Error of %0.3f [%d,%d, %d pts] is too high!\n",
               error0, i1, i2, pt_count);

#ifdef USE_COLORS
        char buf[256];
        camera_params_t cameras_tmp[2] = 
            { model->m_camera0, model->m_camera1 };
        sprintf(buf, "model-%03d-%03d.ply", i1, i2);
        DumpPointsToPly(m_output_directory, buf, pt_count, 2, 
                        model->m_points, colors, cameras_tmp, false);
#endif

        delete [] points;

        ClearKeys(m_image_data[i1]);
        ClearKeys(m_image_data[i2]);

        return false;
    }

    /* Remove outliers and points that are outside the 90% distance
     * threshold */
    double *dists = new double[pt_count];
    for (int i = 0; i < pt_count; i++) {
        dists[i] = v3_magsq(points[i]);
    }

    qsort_ascending();
    double dist_threshold = 
        kth_element_copy(pt_count, iround(0.90 * pt_count), dists);

    dist_threshold = MAX(10000.0 /*500.0*/, dist_threshold);

    if (pt_count < 30)
        dist_threshold = DBL_MAX;

    /* Compute the angle between the rays used to triangulate the
     * point */

    double *angles = new double[pt_count];
        
    for (int i = 0; i < pt_count; i++) {
        if ((int) pt_views[i].size() > 0) {
            double *pos = points[i].p;
            double ray1[3], ray2[3];
            
            matrix_diff(3, 1, 3, 1, pos, cameras[0].t, ray1);
            matrix_diff(3, 1, 3, 1, pos, cameras[1].t, ray2);

            double dot;
            matrix_product(1, 3, 3, 1, ray1, ray2, &dot);
            
            double norm = matrix_norm(3, 1, ray1) * matrix_norm(3, 1, ray2);
            dot /= norm;

            double angle = acos(CLAMP(dot, -1.0 + 1.0e-8, 1.0 - 1.0e-8));
            angles[i] = RAD2DEG(angle);
        } else {
            /* Hack to make sure outliers aren't considered */
            angles[i] = 30.0;
        }
    }

    qsort_descending();
    double angle_threshold = 
        kth_element_copy(pt_count, iround(0.90 * pt_count), angles);
    
    /* Angle threshold should be at most 1.0 -- scratch that, 0.15 */
    angle_threshold = MIN(0.15 /*0.5*/, angle_threshold);

    if (pt_count < 30)
        angle_threshold = 0.0;

    printf("  Using angle threshold: %0.3f\n", angle_threshold);
    fflush(stdout);

    std::vector<ImageKeyVector> pt_views_new;
    int inlier_count = 0;
    for (int i = 0; i < pt_count; i++) {

        if ((int) pt_views[i].size() > 0 && 
            (dists[i] <= dist_threshold && angles[i] >= angle_threshold)) {

            int k1 = pt_views[i][0].second;
            int k2 = pt_views[i][1].second;
            
            m_image_data[i1].m_keys[k1].m_extra = inlier_count;
            m_image_data[i2].m_keys[k2].m_extra = inlier_count;

            pt_views_new.push_back(pt_views[i]);
            points[inlier_count] = points[i];

            if (bundle_from_tracks)
                tracks[inlier_count] = tracks[i];

            inlier_count++;
        } else if (dists[i] > dist_threshold || angles[i] < angle_threshold) {
            int k1 = pt_views[i][0].second;
            int k2 = pt_views[i][1].second;

            double angle = angles[i];

            if (angle < 0.0) {
                printf("Error: angle < 0.0\n");
            }

            m_image_data[i1].m_keys[k1].m_extra = -1;
            m_image_data[i2].m_keys[k2].m_extra = -1;

            printf("Threw out point [%d] with angle %0.3f [dist: %0.3f]\n", 
                   i, angles[i], dists[i]);
            fflush(stdout);
        }
    }

    delete [] dists;
    delete [] angles;

    printf("  pt_count old: %d, new %d\n", pt_count, inlier_count);
    fflush(stdout);

    pt_views = pt_views_new;
    pt_count = inlier_count;

    num_pts_out = pt_count;

    if (pt_count < 20) {
        printf("  Too few points remain, exiting!\n");
        ClearKeys(m_image_data[i1]);
        ClearKeys(m_image_data[i2]);

        return false;
    }

    printf("  focal1, focal2: %0.3f (%0.3f), %0.3f (%0.3f); "
           "%0.3f, %0.3f; %0.3f, %0.3f\n",
           cameras[0].f, cameras[0].constraints[6],
           cameras[1].f, cameras[1].constraints[6],
           cameras[0].k[0], cameras[0].k[1], cameras[1].k[0], cameras[1].k[1]);
    fflush(stdout);

    /* ********** Transform the scene to canonical form ********** */

    double diff[3];
    
    double eye1[3];

    memcpy(eye1, cameras[0].t, sizeof(double) * 3);
    matrix_diff(3, 1, 3, 1, cameras[1].t, cameras[0].t, diff);

    double dist = matrix_norm(3, 1, diff);

    double eye2[3];
    matrix_product(3, 3, 3, 1, cameras[0].R, diff, eye2);
    matrix_scale(3, 1, eye2, 1.0 / dist, eye2);

    double tmp[9];
    matrix_transpose_product2(3, 3, 3, 3, cameras[1].R, cameras[0].R, tmp);
    memcpy(cameras[1].R, tmp, sizeof(double) * 9);

    /* Transform all the scene points */
    for (int i = 0; i < pt_count; i++) {
        double p_tmp[3];
        matrix_diff(3, 1, 3, 1, points[i].p, eye1, p_tmp);
        matrix_product(3, 3, 3, 1, cameras[0].R, p_tmp, points[i].p);
        matrix_scale(3, 1, points[i].p, 1.0 / dist, points[i].p);
    }

    matrix_ident(3, cameras[0].R);
    cameras[0].t[0] = cameras[0].t[1] = cameras[0].t[2] = 0.0;

    memcpy(cameras[1].t, eye2, sizeof(double) * 3);

    /* Now center on the mean of the selected points */
#define CHOOSE_PCT 0.20
    int selected = MIN(MAX(32, iround(CHOOSE_PCT * pt_count)), pt_count);
    selected = MIN(selected, 512);

    v3_t mean = v3_mean(selected, points);

    matrix_diff(3, 1, 3, 1, cameras[0].t, mean.p, cameras[0].t);
    matrix_diff(3, 1, 3, 1, cameras[1].t, mean.p, cameras[1].t);

    for (int i = 0; i < pt_count; i++)
        matrix_diff(3, 1, 3, 1, points[i].p, mean.p, points[i].p);
    
    cameras[0].constraints[0] = cameras[0].t[0];
    cameras[0].constraints[1] = cameras[0].t[1];
    cameras[0].constraints[2] = cameras[0].t[2];

    num_in_back = 0;

    /* Recompute angles */
    int count = 0;
    angle_sum = 0.0;
    for (unsigned int i = 0; i < num_matches; i++) {
	int key_idx1 = matches[i].m_idx1;
	int key_idx2 = matches[i].m_idx2;

        /* Triangulate the point */
        double x_proj1 = GetKey(i1,key_idx1).m_x;
        double y_proj1 = GetKey(i1,key_idx1).m_y;
        double x_proj2 = GetKey(i2,key_idx2).m_x;
        double y_proj2 = GetKey(i2,key_idx2).m_y;

        if (m_optimize_for_fisheye) {
            m_image_data[i1].UndistortPoint(x_proj1, y_proj1, 
                                            x_proj1, y_proj1);
            m_image_data[i2].UndistortPoint(x_proj2, y_proj2, 
                                            x_proj2, y_proj2);
        }
        
        double error;
            
        v2_t p = v2_new(x_proj1, y_proj1);
        v2_t q = v2_new(x_proj2, y_proj2);

        bool in_front = true;
        double angle = 0.0;
        Triangulate(p, q, cameras[0], cameras[1], error, in_front, angle, true);
            
        if (error > 4.0) {
            continue;
        }

        if (!in_front) {
            num_in_back++;
            continue;
        }

        angle_sum += angle;
        count++;
    }
    
    printf("  Average angle [after BA]: %0.3f\n", RAD2DEG(angle_sum / count));
    fflush(stdout);

    angle_out = RAD2DEG(angle_sum / count);

    if (num_in_back > 0.5 * num_matches) {
        printf("  Too many points in back (%d / %d), exiting!\n", 
               num_in_back, num_matches);
        ClearKeys(m_image_data[i1]);
        ClearKeys(m_image_data[i2]);

        return false;
    }

    int cnp = m_estimate_distortion ? 9 : 7;
    double *S = new double[4 * cnp * cnp];
    double *U = new double[2 * cnp * cnp];
    double *Uf = new double[4 * cnp * cnp];
    double *V = new double[pt_count * 3 * 3];
    double *W = new double[pt_count * 2 * 3 * cnp];

#if 1
    /* Set constraints on focal length */
    cameras[0].constrained[6] = true;
    cameras[0].constraints[6] = cameras[0].f;
    cameras[0].weights[6] = 1.0e3;

    cameras[0].constraints[7] = cameras[0].k[0];
    cameras[0].constraints[8] = cameras[0].k[1];
    cameras[0].weights[7] = 1.0e6;
    cameras[0].weights[8] = 1.0e6;

    cameras[1].constrained[6] = true;
    cameras[1].constraints[6] = cameras[1].f;
    cameras[1].weights[6] = 1.0e3;
    cameras[1].constraints[7] = cameras[1].k[0];
    cameras[1].constraints[8] = cameras[1].k[1];
    cameras[1].weights[7] = 1.0e6;
    cameras[1].weights[8] = 1.0e6;
#endif

#if 1
    error0 = RunSFM_SBA(pt_count, 2, 0, false,
                        cameras, points, added_order, NULL, pt_views,
                        TERM_THRESH, S, U, V, W, false);
#endif

    memset(Uf, 0, sizeof(double) * 4 * cnp * cnp);
    for (int i = 0; i < cnp; i++) {
        memcpy(Uf + i * 2 * cnp, U + i * cnp, cnp * sizeof(double));
        memcpy(Uf + (i+cnp) * 2 * cnp + cnp, U + (i+cnp) * cnp, 
               cnp * sizeof(double));
    }

    // double *Vfull = new double[pt_count * pt_count * 3 * 3];
    // double *Vinv = new double[pt_count * pt_count * 3 * 3];
    // memset(Vinv, 0, pt_count * pt_count * 3 * 3 * sizeof(double));
    // memset(Vfull, 0, pt_count * pt_count * 3  * 3 * sizeof(double));

    // #define OLD_GAUGE

#ifdef OLD_GAUGE
    for (int i = 0; i < pt_count; i++) {
        double Vsubinv[9];
        matrix_invert(3, V + 9 * i, Vsubinv);
        
        memcpy(Vinv + (3 * i + 0) * 3 * pt_count + 3 * i, Vsubinv + 0,
               sizeof(double) * 3);
        memcpy(Vinv + (3 * i + 1) * 3 * pt_count + 3 * i, Vsubinv + 3,
               sizeof(double) * 3);
        memcpy(Vinv + (3 * i + 2) * 3 * pt_count + 3 * i, Vsubinv + 6,
               sizeof(double) * 3);
    }
#else

    double *VfA = new double[3 * selected * 3 * selected];
    memset(VfA, 0, selected * selected * 3  * 3 * sizeof(double));

    for (int i = 0; i < selected; i++) {
        memcpy(VfA + (3 * i + 0) * 3 * selected + 3 * i, V + 9 * i + 0,
               sizeof(double) * 3);
        memcpy(VfA + (3 * i + 1) * 3 * selected + 3 * i, V + 9 * i + 3,
               sizeof(double) * 3);
        memcpy(VfA + (3 * i + 2) * 3 * selected + 3 * i, V + 9 * i + 6,
               sizeof(double) * 3);
    }

    /* Augment Vf */
    // mean = v3_mean(selected, points);
    double norm = sqrt(v3_variance_zm(selected, points)); 
        // matrix_norm(3, 1, mean.p);
        // norm = norm * norm;
    printf("  norm2: %0.3e\n", norm);
    fflush(stdout);

    for (int i = 0; i < selected; i++) {
        double dx1 = 1.0e2 * 2.0 * Vx(points[i]) / (norm * selected);
        double dy1 = 1.0e2 * 2.0 * Vy(points[i]) / (norm * selected);
        double dz1 = 1.0e2 * 2.0 * Vz(points[i]) / (norm * selected);

        for (int j = 0; j < selected; j++) {
            double dx2 = 1.0e2 * 2.0 * Vx(points[j]) / (norm * selected);
            double dy2 = 1.0e2 * 2.0 * Vy(points[j]) / (norm * selected);
            double dz2 = 1.0e2 * 2.0 * Vz(points[j]) / (norm * selected);

            VfA[(i * 3 + 0) * 3 * selected + j * 3 + 0] += dx1 * dx2;
            VfA[(i * 3 + 0) * 3 * selected + j * 3 + 1] += dx1 * dy2;
            VfA[(i * 3 + 0) * 3 * selected + j * 3 + 2] += dx1 * dz2;

            VfA[(i * 3 + 1) * 3 * selected + j * 3 + 0] += dy1 * dx2;
            VfA[(i * 3 + 1) * 3 * selected + j * 3 + 1] += dy1 * dy2;
            VfA[(i * 3 + 1) * 3 * selected + j * 3 + 2] += dy1 * dz2;

            VfA[(i * 3 + 2) * 3 * selected + j * 3 + 0] += dz1 * dx2;
            VfA[(i * 3 + 2) * 3 * selected + j * 3 + 1] += dz1 * dy2;
            VfA[(i * 3 + 2) * 3 * selected + j * 3 + 2] += dz1 * dz2;
        }
    }

    double *VinvA = new double[3 * selected * 3 * selected];
    matrix_invert(3 * selected, VfA, VinvA);

    double *Y = new double[pt_count * 2 * 3 * cnp];

    cblas_dgemm_driver_x(2 * cnp, 3 * selected, 3 * selected,
                         3 * pt_count, 3 * selected, 3 * pt_count,
                         W, VinvA, Y);

    double VfB[9], VinvB[9];
    for (int i = selected; i < pt_count; i++) {
        memcpy(VfB + 0, V + 9 * i + 0, 3 * sizeof(double));
        memcpy(VfB + 3, V + 9 * i + 3, 3 * sizeof(double));
        memcpy(VfB + 6, V + 9 * i + 6, 3 * sizeof(double));

        matrix_invert(3, VfB, VinvB);
        
        cblas_dgemm_driver_x(2 * cnp, 3, 3, 
                             3 * pt_count, 3, 3 * pt_count,
                             W + 3 * i, VinvB, Y + 3 * i);
    }
#endif

    double *WViWT = new double[4 * cnp * cnp];
    matrix_transpose_product2(2 * cnp, 3 * pt_count, 2 * cnp, 3 * pt_count,
                              Y, W, WViWT);
    double *Stest = new double[4 * cnp * cnp];
    matrix_diff(2*cnp, 2*cnp, 2*cnp, 2*cnp, Uf, WViWT, Stest);

    /* Add the scale constraint to the hessian S */
    memcpy(eye2, cameras[1].t, sizeof(double) * 3);
    
#ifdef OLD_GAUGE
    double dx = 1.0e4 * 2.0 * eye2[0];
    double dy = 1.0e4 * 2.0 * eye2[1];
    double dz = 1.0e4 * 2.0 * eye2[2];

    S[105] += dx * dx;  S[106] += dx * dy;  S[107] += dx * dz;
    S[119] += dx * dy;  S[120] += dy * dy;  S[121] += dy * dz;
    S[133] += dx * dz;  S[134] += dy * dz;  S[135] += dz * dz;
#else
    memcpy(S, Stest, 4 * cnp * cnp * sizeof(double));
#endif

    double *Ufull = new double[4 * cnp * cnp];
    double *Sfull = new double[2 * cnp];
    double *VTfull = new double[4 * cnp * cnp];

    dgesvd_driver(2*cnp, 2*cnp, S, Ufull, Sfull, VTfull);
    printf("S-values (full):\n");
    for (int i = 0; i < 2*cnp; i++) {
        printf("  [%02d] %0.3e\n", i, Sfull[i]);
    }    

    double *Sinv = new double[4*cnp*cnp];
    matrix_invert(2*cnp, S, Sinv);

    int row1 = 2 * cnp * (cnp+0) + cnp;
    int row2 = 2 * cnp * (cnp+1) + cnp;
    int row3 = 2 * cnp * (cnp+2) + cnp;
    
    double C2[9] = { Sinv[row1+0], Sinv[row1+1], Sinv[row1+2],
                     Sinv[row2+0], Sinv[row2+1], Sinv[row2+2],
                     Sinv[row3+0], Sinv[row3+1], Sinv[row3+2] };

    printf("C2:\n");
    matrix_print(3, 3, C2);
    fflush(stdout);

    double U2[9], S2[3], VT2[9];
    dgesvd_driver(3, 3, C2, U2, S2, VT2);

    bool sym_error = false;
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            double ij = C2[i * 3 + j];
            double ji = C2[j * 3 + i];
            
            if (fabs(ij - ji) > 1.0e-2) {
                printf("C2: Symmetry error (%d,%d)!\n", i, j);
                sym_error = true;
            }
        }
    }

    printf("Singular values (C2): %0.3e, %0.3e, %0.3e\n", S2[0], S2[1], S2[2]);
    fflush(stdout);

    if (model != NULL) {
        model->m_error = error0;
        memcpy(model->m_C1, C2, sizeof(double) * 9);
        model->m_points = new v3_t[pt_count];
        memcpy(model->m_points, points, sizeof(v3_t) * pt_count);
        model->m_camera0 = cameras[0];
        model->m_camera1 = cameras[1];

        /* Reflect the model */
        for (int i = 0; i < pt_count; i++) {
            // Vz(model->m_points[i]) = -Vz(model->m_points[i]);
            Vz(model->m_points[i]) = Vz(model->m_points[i]);
        }

#if 0
        model->m_camera0.R[2] = -model->m_camera0.R[2];
        model->m_camera0.R[5] = -model->m_camera0.R[5];
        model->m_camera0.R[6] = -model->m_camera0.R[6];
        model->m_camera0.R[7] = -model->m_camera0.R[7];
        model->m_camera0.t[2] = -model->m_camera0.t[2];

        model->m_camera1.R[2] = -model->m_camera1.R[2];
        model->m_camera1.R[5] = -model->m_camera1.R[5];
        model->m_camera1.R[6] = -model->m_camera1.R[6];
        model->m_camera1.R[7] = -model->m_camera1.R[7];
        model->m_camera1.t[2] = -model->m_camera1.t[2];
#endif
    }

    /* Now clear constraints */

    cameras[0].constrained[0] = false;
    cameras[0].constrained[1] = false;
    cameras[0].constrained[2] = false;
    cameras[0].constrained[3] = false;
    cameras[0].constrained[4] = false;
    cameras[0].constrained[5] = false;

    cameras[1].constrained[0] = true;
    cameras[1].constrained[1] = true;
    cameras[1].constrained[2] = true;
    cameras[1].constrained[3] = true;
    cameras[1].constrained[4] = true;
    cameras[1].constrained[5] = true;

    cameras[1].constraints[0] = 0.0;
    cameras[1].constraints[1] = 0.0;
    cameras[1].constraints[2] = 0.0;
    cameras[1].constraints[3] = 0.0;
    cameras[1].constraints[4] = 0.0;
    cameras[1].constraints[5] = 0.0;

    cameras[1].weights[0] = 1.0e6;
    cameras[1].weights[1] = 1.0e6;
    cameras[1].weights[2] = 1.0e6;
    cameras[1].weights[3] = 1.0e6;
    cameras[1].weights[4] = 1.0e6;
    cameras[1].weights[5] = 1.0e6;


    /* ********** Transform the scene to canonical form ********** */

    /* Find the scale factor that places the two views a unit apart */
    memcpy(eye2, cameras[1].t, sizeof(double) * 3);
    matrix_diff(3, 1, 3, 1, cameras[0].t, cameras[1].t, diff);

#if 0
    dist = matrix_norm(3, 1, diff);

    matrix_product(3, 3, 3, 1, cameras[1].R, diff, eye1);
    matrix_scale(3, 1, eye1, 1.0 / dist, eye1);

    matrix_transpose_product2(3, 3, 3, 3, cameras[0].R, cameras[1].R, tmp);
    memcpy(cameras[0].R, tmp, sizeof(double) * 9);

    /* Transform all the scene points */
    for (int i = 0; i < pt_count; i++) {
        double p_tmp[3];
        matrix_diff(3, 1, 3, 1, points[i].p, eye2, p_tmp);
        matrix_product(3, 3, 3, 1, cameras[1].R, p_tmp, points[i].p);
        matrix_scale(3, 1, points[i].p, 1.0 / dist, points[i].p);
    }

    matrix_ident(3, cameras[1].R);
    cameras[1].t[0] = cameras[1].t[1] = cameras[1].t[2] = 0.0;

    if (m_explicit_camera_centers) {
        memcpy(cameras[0].t, eye1, sizeof(double) * 3);
    } else {
        matrix_product(3, 3, 3, 1, cameras[0].R, eye1, cameras[0].t);
        matrix_scale(3, 1, cameras[0].t, -1.0, cameras[0].t);
    }
#endif

    /* Subtract out point mean */

    mean = v3_mean(selected, points);

    matrix_diff(3, 1, 3, 1, cameras[0].t, mean.p, cameras[0].t);
    matrix_diff(3, 1, 3, 1, cameras[1].t, mean.p, cameras[1].t);

    for (int i = 0; i < pt_count; i++)
        matrix_diff(3, 1, 3, 1, points[i].p, mean.p, points[i].p);

    cameras[1].constraints[0] = cameras[1].t[0];
    cameras[1].constraints[1] = cameras[1].t[1];
    cameras[1].constraints[2] = cameras[1].t[2];

#if 1
#ifdef USE_COLORS
    error0 = RunSFM_SBA(pt_count, 2, 0, false,
                        cameras, points, added_order, colors, pt_views, 
                        TERM_THRESH, S, U, V, W, false);
#else
    error0 = RunSFM_SBA(pt_count, 2, 0, false,
                        cameras, points, added_order, NULL, pt_views,
                        TERM_THRESH, S, U, V, W, false);
#endif
#endif

    memset(Uf, 0, sizeof(double) * 4 * cnp * cnp);
    for (int i = 0; i < cnp; i++) {
        memcpy(Uf + i * 2 * cnp, U + i * cnp, cnp * sizeof(double));
        memcpy(Uf + (i+cnp) * 2 * cnp + cnp, U + (i+cnp) * cnp, 
               cnp * sizeof(double));
    }

    // Vinv = new double[pt_count * pt_count * 3 * 3];
    // memset(Vinv, 0, pt_count * pt_count * 3 * 3 * sizeof(double));

#ifdef OLD_GAUGE
    for (int i = 0; i < pt_count; i++) {
        double Vsubinv[9];
        matrix_invert(3, V + 9 * i, Vsubinv);
        
        memcpy(Vinv + (3 * i + 0) * 3 * pt_count + 3 * i, Vsubinv + 0,
               sizeof(double) * 3);
        memcpy(Vinv + (3 * i + 1) * 3 * pt_count + 3 * i, Vsubinv + 3,
               sizeof(double) * 3);
        memcpy(Vinv + (3 * i + 2) * 3 * pt_count + 3 * i, Vsubinv + 6,
               sizeof(double) * 3);
    }
#else
    memset(VfA, 0, selected * selected * 3  * 3 * sizeof(double));

    for (int i = 0; i < selected; i++) {
        memcpy(VfA + (3 * i + 0) * 3 * selected + 3 * i, V + 9 * i + 0,
               sizeof(double) * 3);
        memcpy(VfA + (3 * i + 1) * 3 * selected + 3 * i, V + 9 * i + 3,
               sizeof(double) * 3);
        memcpy(VfA + (3 * i + 2) * 3 * selected + 3 * i, V + 9 * i + 6,
               sizeof(double) * 3);
    }

    mean = v3_mean(selected, points);

    norm = sqrt(v3_variance_zm(selected, points)); 
    // norm = matrix_norm(3, 1, mean.p);
    // norm = norm * norm;
    printf("  norm1: %0.3e\n", norm);

    /* Augment Vf */
    for (int i = 0; i < selected; i++) {
        double dx1 = 1.0e2 * 2.0 * Vx(points[i]) / (norm * selected);
        double dy1 = 1.0e2 * 2.0 * Vy(points[i]) / (norm * selected);
        double dz1 = 1.0e2 * 2.0 * Vz(points[i]) / (norm * selected);

        for (int j = 0; j < selected; j++) {
            double dx2 = 1.0e2 * 2.0 * Vx(points[j]) / (norm * selected);
            double dy2 = 1.0e2 * 2.0 * Vy(points[j]) / (norm * selected);
            double dz2 = 1.0e2 * 2.0 * Vz(points[j]) / (norm * selected);

            VfA[(i * 3 + 0) * 3 * selected + j * 3 + 0] += dx1 * dx2;
            VfA[(i * 3 + 0) * 3 * selected + j * 3 + 1] += dx1 * dy2;
            VfA[(i * 3 + 0) * 3 * selected + j * 3 + 2] += dx1 * dz2;

            VfA[(i * 3 + 1) * 3 * selected + j * 3 + 0] += dy1 * dx2;
            VfA[(i * 3 + 1) * 3 * selected + j * 3 + 1] += dy1 * dy2;
            VfA[(i * 3 + 1) * 3 * selected + j * 3 + 2] += dy1 * dz2;

            VfA[(i * 3 + 2) * 3 * selected + j * 3 + 0] += dz1 * dx2;
            VfA[(i * 3 + 2) * 3 * selected + j * 3 + 1] += dz1 * dy2;
            VfA[(i * 3 + 2) * 3 * selected + j * 3 + 2] += dz1 * dz2;
        }
    }

    matrix_invert(3 * selected, VfA, VinvA);

    cblas_dgemm_driver_x(2 * cnp, 3 * selected, 3 * selected,
                         3 * pt_count, 3 * selected, 3 * pt_count,
                         W, VinvA, Y);

    for (int i = selected; i < pt_count; i++) {
        memcpy(VfB + 0, V + 9 * i + 0, 3 * sizeof(double));
        memcpy(VfB + 3, V + 9 * i + 3, 3 * sizeof(double));
        memcpy(VfB + 6, V + 9 * i + 6, 3 * sizeof(double));

        matrix_invert(3, VfB, VinvB);
        
        cblas_dgemm_driver_x(2 * cnp, 3, 3, 
                             3 * pt_count, 3, 3 * pt_count,
                             W + 3 * i, VinvB, Y + 3 * i);
    }
#endif

    matrix_transpose_product2(2 * cnp, 3 * pt_count, 2 * cnp, 3 * pt_count,
                              Y, W, WViWT);
    matrix_diff(2*cnp, 2*cnp, 2*cnp, 2*cnp, Uf, WViWT, Stest);

    /* Add the scale constraint to the hessian S */
    memcpy(eye1, cameras[0].t, sizeof(double) * 3);
    
#ifdef OLD_GAUGE
    double dx = 1.0e4 * 2.0 * eye1[0];
    double dy = 1.0e4 * 2.0 * eye1[1];
    double dz = 1.0e4 * 2.0 * eye1[2];

    S[0] += dx * dx;   S[1] += dx * dy;   S[2] += dx * dz;
    S[14] += dx * dy;  S[15] += dy * dy;  S[16] += dy * dz;
    S[28] += dx * dz;  S[29] += dy * dz;  S[30] += dz * dz;
#else
    memcpy(S, Stest, 4 * cnp * cnp * sizeof(double));
#endif

    dgesvd_driver(2*cnp, 2*cnp, S, Ufull, Sfull, VTfull);
    printf("Singular values (full):\n");
    for (int i = 0; i < 2*cnp; i++) {
        printf("  [%02d] %0.3e\n", i, Sfull[i]);
    }    

    matrix_invert(2*cnp, S, Sinv);

    row1 = 0;
    row2 = 2 * cnp;
    row3 = 4 * cnp;

    double C1[9] = { Sinv[row1+0], Sinv[row1+1], Sinv[row1+2],
                     Sinv[row2+0], Sinv[row2+1], Sinv[row2+2],
                     Sinv[row3+0], Sinv[row3+1], Sinv[row3+2] };

    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            double ij = C1[i * 3 + j];
            double ji = C1[j * 3 + i];
            
            if (fabs(ij - ji) > 1.0e-2) {
                printf("C1: Symmetry error (%d,%d)!\n", i, j);
                sym_error = true;
            }
        }
    }

#if 0
    v2_t *projs1 = new v2_t[pt_count];
    v2_t *projs2 = new v2_t[pt_count];
    int num_vars = 8 + 3 * pt_count;
    double *H = new double[num_vars * num_vars];

    for (int i = 0; i < pt_count; i++) {
        Vx(projs1[i]) = m_image_data[i1].m_keys[pt_views[i][0].second].m_x;
        Vy(projs1[i]) = m_image_data[i1].m_keys[pt_views[i][0].second].m_y;
        Vx(projs2[i]) = m_image_data[i2].m_keys[pt_views[i][1].second].m_x;
        Vy(projs2[i]) = m_image_data[i2].m_keys[pt_views[i][1].second].m_y;
    }
    
    camera_refine_fix_free(pt_count, points, projs2, projs1,
                           cameras + 1, cameras + 0, H);

    double *Hinv = new double[num_vars * num_vars];
    matrix_invert(num_vars, H, Hinv);

    memcpy(C1 + 0, Hinv + 1 * num_vars + 1, 3 * sizeof(double));
    memcpy(C1 + 3, Hinv + 2 * num_vars + 1, 3 * sizeof(double));
    memcpy(C1 + 6, Hinv + 3 * num_vars + 1, 3 * sizeof(double));
#endif

    printf("C1:\n");
    matrix_print(3, 3, C1);
    fflush(stdout);

    double U1[9], S1[3], VT1[9];
    dgesvd_driver(3, 3, C1, U1, S1, VT1);

    printf("Singular values (C1): %0.3e, %0.3e, %0.3e\n", 
           S1[0], S1[1], S1[2]);
    fflush(stdout);

    if (model != NULL) {
        memcpy(model->m_C0, C1, sizeof(double) * 9);

        if (bundle_from_tracks) {
            double *tracks_tmp = new double[pt_count];
            int *perm = new int[pt_count];
            v3_t *points_tmp = new v3_t[pt_count];

            for (int i = 0; i < pt_count; i++) {
                tracks_tmp[i] = (double) tracks[i];
                points_tmp[i] = model->m_points[i];
            }

            qsort_ascending();
            qsort_perm(pt_count, tracks_tmp, perm);

            model->m_tracks = new int[pt_count];
            model->m_keys1 = model->m_keys2 = NULL;

            for (int i = 0; i < pt_count; i++) {
                model->m_tracks[i] = iround(tracks_tmp[i]);
                model->m_points[i] = points_tmp[perm[i]];
            }

            delete [] tracks_tmp;
            delete [] perm;
            delete [] points_tmp;
        } else {
            model->m_tracks = NULL;
            model->m_keys1 = new int[pt_count];
            model->m_keys2 = new int[pt_count];
            
            for (int i = 0; i < pt_count; i++) {
                int k1 = pt_views[i][0].second;
                int k2 = pt_views[i][1].second;
                model->m_keys1[i] = k1;
                model->m_keys2[i] = k2;
            }
        }

        model->m_num_points = pt_count;
        model->m_angle = RAD2DEG(angle_sum / count);
    }

#ifdef USE_COLORS
    char buf[256];
    camera_params_t cameras_tmp[2] = 
        { model->m_camera0, model->m_camera1 };
    sprintf(buf, "model-%03d-%03d.ply", i1, i2);
    DumpPointsToPly(m_output_directory, buf, pt_count, 2, 
                    model->m_points, colors, cameras_tmp, false);
#endif

#if 0
    delete [] H;
    delete [] Hinv;

    delete [] projs1;
    delete [] projs2;
#endif

    delete [] points;

    delete [] S;
    delete [] U;
    delete [] Uf;
    delete [] V;
    delete [] W;
    delete [] Y;
    delete [] Sinv;

    delete [] Ufull;
    delete [] Sfull;
    delete [] VTfull;

    delete [] WViWT;
    delete [] Stest;

    // delete [] Vinv;

#ifndef OLD_GAUGE
    // delete [] Vf;
    delete [] VfA;
    delete [] VinvA;
#endif

    ClearKeys(m_image_data[i1]);
    ClearKeys(m_image_data[i2]);

    if (sym_error)
        return false;
 
    return true;
}

void BundlerApp::ComputeCameraCovariance()
{
    int num_images = GetNumImages();

    int cnp = (m_estimate_distortion) ? 9 : 7;
    int num_vars = cnp * num_images;
    double *S = new double[num_vars * num_vars];

    /* Add constraints */
    m_use_point_constraints = true;
    m_point_constraint_weight = 1000.0;

    int num_points = (int) m_point_data.size();
    m_point_constraints = new v3_t[num_points];
    
    for (int i = 0; i < num_points; i++) {
        PointData &p = m_point_data[i];
        m_point_constraints[i] = v3_new(p.m_pos[0], p.m_pos[1], p.m_pos[2]);
    }

    ReRunSFM(S);

#if 0
    /* Debugging stuff */
    for (int i = 1; i < num_images + 1; i++) {
        int num_vars_sub = i * cnp;
        double *Ssub = new double[num_vars_sub * num_vars_sub];

        for (int j = 0; j < num_vars_sub; j++) {
            memcpy(Ssub + j * num_vars_sub, S + j * num_vars, 
                   sizeof(double) * num_vars_sub);
        }

        double *U, *S, *VT;
        U = new double[num_vars_sub * num_vars_sub];
        S = new double[num_vars_sub];
        VT = new double[num_vars_sub * num_vars_sub];
        
        dgesvd_driver(num_vars_sub, num_vars_sub, Ssub, U, S, VT);
        printf("Sing. values [round %d]\n", i);
        for (int j = 0; j < num_vars_sub; j++) {
            printf("   S[%d] = %0.5e\n", j, S[j]);
        }

        printf("Smallest eigenvector [round %d]:\n", i);
        matrix_print(i, cnp, VT + (num_vars_sub-1) * num_vars_sub);

        delete [] U;
        delete [] S;
        delete [] VT;

        double *Sinv_sub = new double[num_vars_sub * num_vars_sub];
        matrix_invert(num_vars_sub, Ssub, Sinv_sub);

        for (int j = 0; j < i; j++) {
            int row0 = (cnp * j + 0) * num_vars_sub;
            int row1 = (cnp * j + 1) * num_vars_sub;
            int row2 = (cnp * j + 2) * num_vars_sub;
            int off = cnp * j;

            double C[9] = { Sinv_sub[row0+off+0], Sinv_sub[row0+off+1], Sinv_sub[row0+off+2],
                            Sinv_sub[row1+off+0], Sinv_sub[row1+off+1], Sinv_sub[row1+off+2],
                            Sinv_sub[row2+off+0], Sinv_sub[row2+off+1], Sinv_sub[row2+off+2] };
            
            printf("Covariance[%d] [Image %d]:\n", i, j);
            matrix_print(3, 3, C);
        }

        delete [] Ssub;
        delete [] Sinv_sub;
    }
#endif

    /* Invert S */
    double *Sinv = new double[num_vars * num_vars];
    matrix_invert(num_vars, S, Sinv);
    
    FILE *f = fopen("covariance.txt", "w");
    if (f == NULL) {
        printf("[ComputeCameraCovariance] Error opening file %s for writing\n",
               "covariance.txt");
        return;
    }

    int count = 0;
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted) {
            int row0 = (cnp * count + 0) * num_vars;
            int row1 = (cnp * count + 1) * num_vars;
            int row2 = (cnp * count + 2) * num_vars;
            int off = cnp * count;

            double C[9] = 
                { Sinv[row0+off+0], Sinv[row0+off+1], Sinv[row0+off+2],
                  Sinv[row1+off+0], Sinv[row1+off+1], Sinv[row1+off+2],
                  Sinv[row2+off+0], Sinv[row2+off+1], Sinv[row2+off+2] };

            fprintf(f, "%d\n", i);
            fprintf(f, "%0.6e %0.6e %0.6e "
                    "%0.6e %0.6e %0.6e "
                    "%0.6e %0.6e %0.6e\n", 
                    C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8]);
            fprintf(f, "%0.6e\n", C[0] + C[4] + C[8]);

            printf("Covariance [Image %d]:\n", i);
            matrix_print(3, 3, C);

            count++;
        }
    }
}

#if 0
void BundlerApp::ComputeCameraCovariance()
{
    if (m_covariance_fix1 == -1 || m_covariance_fix2 == -1) {
        printf("[ComputeCameraCovariance] "
               "Two images need to be specified.\n");
    }
    
    int num_images = GetNumImages();

    int num_vars = 7 * num_images;
    double *S = new double[num_vars * num_vars];

    /* Transform the scene to a canonical form */
    int c1 = m_covariance_fix1;
    int c2 = m_covariance_fix2;

    // RepositionScene();
    TransformSceneCanonical(c1, c2);

    /* Add constraints */
    m_image_data[c1].m_camera.m_constrained[0] = true;
    m_image_data[c1].m_camera.m_constrained[1] = true;
    m_image_data[c1].m_camera.m_constrained[2] = true;
    m_image_data[c1].m_camera.m_constrained[3] = true;
    m_image_data[c1].m_camera.m_constrained[4] = true;
    m_image_data[c1].m_camera.m_constrained[5] = true;

    m_image_data[c1].m_camera.m_constraints[0] = 0.0;
    m_image_data[c1].m_camera.m_constraints[1] = 0.0;
    m_image_data[c1].m_camera.m_constraints[2] = 0.0;
    m_image_data[c1].m_camera.m_constraints[3] = 0.0;
    m_image_data[c1].m_camera.m_constraints[4] = 0.0;
    m_image_data[c1].m_camera.m_constraints[5] = 0.0;

    m_image_data[c1].m_camera.m_constraint_weights[0] = 1.0e6;
    m_image_data[c1].m_camera.m_constraint_weights[1] = 1.0e6;
    m_image_data[c1].m_camera.m_constraint_weights[2] = 1.0e6;
    m_image_data[c1].m_camera.m_constraint_weights[3] = 1.0e6;
    m_image_data[c1].m_camera.m_constraint_weights[4] = 1.0e6;
    m_image_data[c1].m_camera.m_constraint_weights[5] = 1.0e6;

    ReRunSFM(S);
    printf("  focal1, focal2: %0.3f (%0.3f), %0.3f (%0.3f)\n",
           m_image_data[c1].m_camera.m_focal,m_image_data[c1].m_init_focal,
           m_image_data[c2].m_camera.m_focal,m_image_data[c2].m_init_focal);

    /* Add constraints to camera 2 */
    double eye[3];
    m_image_data[c2].m_camera.GetPosition(eye);
    eye[2] *= -1.0;

    double dx = 1.0e4 * 2.0 * eye[0];
    double dy = 1.0e4 * 2.0 * eye[1];
    double dz = 1.0e4 * 2.0 * eye[2];

    int row0 = (7 * c2 + 0) * num_vars;
    int row1 = (7 * c2 + 1) * num_vars;
    int row2 = (7 * c2 + 2) * num_vars;
    int off = 7 * c2;

    S[row0 + off + 0] += dx * dx;   
    S[row0 + off + 1] += dx * dy;
    S[row0 + off + 2] += dx * dz;

    S[row1 + off + 0] += dx * dy;  
    S[row1 + off + 1] += dy * dy;  
    S[row1 + off + 2] += dy * dz;

    S[row2 + off + 0] += dx * dz;  
    S[row2 + off + 1] += dy * dz;  
    S[row2 + off + 2] += dz * dz;

#if 0
    /* Debugging stuff */
    for (int i = 1; i < num_images + 1; i++) {
        int num_vars_sub = i * 7;
        double *Ssub = new double[num_vars_sub * num_vars_sub];

        for (int j = 0; j < num_vars_sub; j++) {
            memcpy(Ssub + j * num_vars_sub, S + j * num_vars, 
                   sizeof(double) * num_vars_sub);
        }

        double *U, *S, *VT;
        U = new double[num_vars_sub * num_vars_sub];
        S = new double[num_vars_sub];
        VT = new double[num_vars_sub * num_vars_sub];
        
        dgesvd_driver(num_vars_sub, num_vars_sub, Ssub, U, S, VT);
        printf("Sing. values [round %d]\n", i);
        for (int j = 0; j < num_vars_sub; j++) {
            printf("   S[%d] = %0.5e\n", j, S[j]);
        }

        printf("Smallest eigenvector [round %d]:\n", i);
        matrix_print(i, 7, VT + (num_vars_sub-1) * num_vars_sub);

        delete [] U;
        delete [] S;
        delete [] VT;

        double *Sinv_sub = new double[num_vars_sub * num_vars_sub];
        matrix_invert(num_vars_sub, Ssub, Sinv_sub);

        for (int j = 0; j < i; j++) {
            int row0 = (7 * j + 0) * num_vars_sub;
            int row1 = (7 * j + 1) * num_vars_sub;
            int row2 = (7 * j + 2) * num_vars_sub;
            int off = 7 * j;

            double C[9] = { Sinv_sub[row0+off+0], Sinv_sub[row0+off+1], Sinv_sub[row0+off+2],
                            Sinv_sub[row1+off+0], Sinv_sub[row1+off+1], Sinv_sub[row1+off+2],
                            Sinv_sub[row2+off+0], Sinv_sub[row2+off+1], Sinv_sub[row2+off+2] };
            
            printf("Covariance[%d] [Image %d]:\n", i, j);
            matrix_print(3, 3, C);
        }

        delete [] Ssub;
        delete [] Sinv_sub;
    }
#endif

    /* Invert S */
    double *Sinv = new double[num_vars * num_vars];
    matrix_invert(num_vars, S, Sinv);
    
    FILE *f = fopen("covariance.txt", "w");
    if (f == NULL) {
        printf("[ComputeCameraCovariance] Error opening file %s for writing\n",
               "covariance.txt");
        return;
    }

    int count = 0;
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted) {
            int row0 = (7 * count + 0) * num_vars;
            int row1 = (7 * count + 1) * num_vars;
            int row2 = (7 * count + 2) * num_vars;
            int off = 7 * count;

            double C[9] = 
                { Sinv[row0+off+0], Sinv[row0+off+1], Sinv[row0+off+2],
                  Sinv[row1+off+0], Sinv[row1+off+1], Sinv[row1+off+2],
                  Sinv[row2+off+0], Sinv[row2+off+1], Sinv[row2+off+2] };

            printf("Covariance [Image %d]:\n", i);
            matrix_print(3, 3, C);

            fprintf(f, "%d\n", i);
            fprintf(f, "%0.6f %0.6f %0.6f "
                    "%0.6f %0.6f %0.6f "
                    "%0.6f %0.6f %0.6f\n", 
                    C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8]);
            fprintf(f, "%0.6f\n", C[0] + C[4] + C[8]);

            count++;
        }
    }

    fclose(f);
}
#endif
