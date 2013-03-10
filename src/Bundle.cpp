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

/* Bundle.cpp */
/* Bundle adjustment routines */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <vector>

#include "defines.h"
#include "horn.h"
#include "image.h"
#include "matrix.h"
#include "qsort.h"
#include "resample.h"
#include "sfm.h"
#include "triangulate.h"
#include "util.h"

#include "BundlerApp.h"
// #include "SifterGraph.h"
#include "BundleAdd.h"
#include "Epipolar.h"
#include "Distortion.h"

/* Use a 180 rotation to fix up the intrinsic matrix */
void FixIntrinsics(double *P, double *K, double *R, double *t) 
{
    /* Check the parity along the diagonal */
    int neg = (K[0] < 0.0) + (K[4] < 0.0) + (K[8] < 0.0);

    /* If odd parity, negate the instrinsic matrix */
    if ((neg % 2) == 1) {
        matrix_scale(3, 3, K, -1.0, K);
        matrix_scale(3, 4, P, -1.0, P);
    }

    /* Now deal with case of even parity */
    double fix[9];
    matrix_ident(3, fix);
    double tmp[9], tmp2[12];

    if (K[0] < 0.0 && K[4] < 0.0) {
        fix[0] = -1.0;
        fix[4] = -1.0;
    } else if (K[0] < 0.0) {
        fix[0] = -1.0;
        fix[8] = -1.0;
    } else if (K[4] < 0.0) {
        fix[4] = -1.0;
        fix[8] = -1.0;
    } else {
        /* No change needed */
    }

    matrix_product(3, 3, 3, 3, K, fix, tmp);
    memcpy(K, tmp, sizeof(double) * 3 * 3);

    double Kinv[9];
    matrix_invert(3, K, Kinv);

    matrix_product(3, 3, 3, 4, Kinv, P, tmp2);

    memcpy(R + 0, tmp2 + 0, sizeof(double) * 3);
    memcpy(R + 3, tmp2 + 4, sizeof(double) * 3);
    memcpy(R + 6, tmp2 + 8, sizeof(double) * 3);

    t[0] = tmp2[3];
    t[1] = tmp2[7];
    t[2] = tmp2[11];
}

void GetIntrinsics(const camera_params_t &camera, double *K) {
    if (!camera.known_intrinsics) {
        K[0] = camera.f;  K[1] = 0.0;       K[2] = 0.0;
        K[3] = 0.0;       K[4] = camera.f;  K[5] = 0.0;
        K[6] = 0.0;       K[7] = 0.0;       K[8] = 1.0;    
    } else {
        memcpy(K, camera.K_known, 9 * sizeof(double));
    }
}

/* Compute the angle between two rays */
double ComputeRayAngle(v2_t p, v2_t q, 
                       const camera_params_t &cam1, 
                       const camera_params_t &cam2)
{
    double K1[9], K2[9];
    GetIntrinsics(cam1, K1);
    GetIntrinsics(cam2, K2);

    double K1_inv[9], K2_inv[9];
    matrix_invert(3, K1, K1_inv);
    matrix_invert(3, K2, K2_inv);

    double p3[3] = { Vx(p), Vy(p), 1.0 };
    double q3[3] = { Vx(q), Vy(q), 1.0 };

    double p3_norm[3], q3_norm[3];
    matrix_product331(K1_inv, p3, p3_norm);
    matrix_product331(K2_inv, q3, q3_norm);

    v2_t p_norm = v2_new(p3_norm[0] / p3_norm[2], p3_norm[1] / p3_norm[2]);
    v2_t q_norm = v2_new(q3_norm[0] / q3_norm[2], q3_norm[1] / q3_norm[2]);

    double R1_inv[9], R2_inv[9];
    matrix_transpose(3, 3, (double *) cam1.R, R1_inv);
    matrix_transpose(3, 3, (double *) cam2.R, R2_inv);

    double p_w[3], q_w[3];

    double pv[3] = { Vx(p_norm), Vy(p_norm), -1.0 };
    double qv[3] = { Vx(q_norm), Vy(q_norm), -1.0 };

    double Rpv[3], Rqv[3];

    matrix_product331(R1_inv, pv, Rpv);
    matrix_product331(R2_inv, qv, Rqv);

    matrix_sum(3, 1, 3, 1, Rpv, (double *) cam1.t, p_w);
    matrix_sum(3, 1, 3, 1, Rqv, (double *) cam2.t, q_w);

    /* Subtract out the camera center */
    double p_vec[3], q_vec[3];
    matrix_diff(3, 1, 3, 1, p_w, (double *) cam1.t, p_vec);
    matrix_diff(3, 1, 3, 1, q_w, (double *) cam2.t, q_vec);

    /* Compute the angle between the rays */
    double dot;
    matrix_product(1, 3, 3, 1, p_vec, q_vec, &dot);

    double mag = matrix_norm(3, 1, p_vec) * matrix_norm(3, 1, q_vec);

    return acos(CLAMP(dot / mag, -1.0 + 1.0e-8, 1.0 - 1.0e-8));
}

/* Utility routine for updating an image key vector */
static void RemoveImageKey(ImageKeyVector &vec, int view) 
{
    int size = (int) vec.size();
    int found = 0;

    for (int i = 0; i < size; i++) {
        if (vec[i].first == view) {
            vec.erase(vec.begin() + i);
            i--;
            size--;
            found++;
        }
    }

    if (found == 0)
        printf("[RemoveImageKey] Error! Couldn't find view %d\n", view);
    else if (found > 1) 
        printf("[RemoveImageKey] Error! Found too many views %d\n", view);
}

/* Check cheirality for a camera and a point */
bool CheckCheirality(v3_t p, const camera_params_t &camera) 
{
    double pt[3] = { Vx(p), Vy(p), Vz(p) };
    double cam[3];

    pt[0] -= camera.t[0];
    pt[1] -= camera.t[1];
    pt[2] -= camera.t[2];
    matrix_product(3, 3, 3, 1, (double *) camera.R, pt, cam);

    // EDIT!!!
    if (cam[2] > 0.0)
        return false;
    else
        return true;
}

double GetCameraDistance(camera_params_t *c1, camera_params_t *c2) 
{
    double center1[3]; 
    double Rinv1[9];
    matrix_invert(3, c1->R, Rinv1);

    memcpy(center1, c1->t, 3 * sizeof(double));

    double center2[3];
    double Rinv2[9];
    matrix_invert(3, c2->R, Rinv2);

    memcpy(center2, c2->t, 3 * sizeof(double));

    double dx = center1[0] - center2[0];
    double dy = center1[1] - center2[1];
    double dz = center1[2] - center2[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}

void InitializeCameraParams(const ImageData &data, camera_params_t &camera)
{
    matrix_ident(3, camera.R);
    camera.t[0] = camera.t[1] = camera.t[2] = 0.0;
    camera.f = 0.0;
    camera.k[0] = camera.k[1] = 0.0;

    camera.k_inv[0] = camera.k_inv[2] = camera.k_inv[3] = 0.0;
    camera.k_inv[4] = camera.k_inv[5] = 0.0;
    camera.k_inv[1] = 1.0;

    camera.f_scale = 1.0;
    camera.k_scale = 1.0;

    for (int i = 0; i < NUM_CAMERA_PARAMS; i++) {
        camera.constrained[i] = 0;
        camera.constraints[i] = 0.0;
        camera.weights[i] = 0.0;
    }

    camera.fisheye = data.m_fisheye;
    camera.f_cx = data.m_fCx;
    camera.f_cy = data.m_fCy;
    camera.f_rad = data.m_fRad;
    camera.f_angle = data.m_fAngle;
    camera.f_focal = data.m_fFocal;

    if (data.m_known_intrinsics) {
        camera.known_intrinsics = 1;
        memcpy(camera.K_known, data.m_K, 9 * sizeof(double));
        memcpy(camera.k_known, data.m_k, 5 * sizeof(double));
    } else {
        camera.known_intrinsics = 0;
    }
}

void BundlerApp::
CheckPointKeyConsistency(const std::vector<ImageKeyVector> pt_views,
                         int *added_order)
{
    int num_points = (int) pt_views.size();
    int errors = 0;

    for (int i = 0; i < num_points; i++) {
        int num_views = pt_views[i].size();

        for (int j = 0; j < num_views; j++) {
            /* check consistency */
            ImageKey ik = pt_views[i][j];

            int img = added_order[ik.first];
            int key = ik.second;

            if (m_image_data[img].m_keys[key].m_extra != i) {
                printf("[CheckPointKeyConsistency] Error: (%d,%d).m_extra "
                    "should be %d (m_extra is %d)\n", img, key, i,
                    m_image_data[img].m_keys[key].m_extra);
                errors++;
            }
        }
    }

    printf("[CheckPointKeyConsistency] There were %d errors\n", errors);
}

void BundlerApp::ReRunSFM(double *S, double *U, double *V, double *W) 
{
    // #define RERUN_ADD_POINTS
#ifdef RERUN_ADD_POINTS
    /* Compute initial image information */
    ComputeGeometricConstraints();
#endif

    int num_pts = (int) m_point_data.size();
    int num_images = GetNumImages();

    /* Initialize all keypoints to have not been matched */
    printf("[SifterApp::ReRunSFM] Initializing keypoints...\n");
    for (int i = 0; i < num_images; i++) {
        int num_keys = GetNumKeys(i);
        for (int j = 0; j < num_keys; j++) {
            GetKey(i,j).m_extra = -1;
        }
    }

    /* Set up the cameras */
    int num_init_cams = 0;
    int *added_order = new int[num_images];
    int *added_order_inv = new int[num_images];
    std::vector<ImageKeyVector> pt_views;

    camera_params_t *cameras = new camera_params_t[num_images];

    printf("[SifterApp::ReRunSFM] Setting up cameras\n");
    for (int i = 0; i < num_images; i++) {
        printf(".");
        fflush(stdout);

        if (m_image_data[i].m_camera.m_adjusted) {
            m_image_data[i].LoadKeys(false, !m_optimize_for_fisheye);

#ifdef RERUN_ADD_POINTS
            m_image_data[i].ReadKeyColors();
            SetTracks(i);
#endif

            added_order[num_init_cams] = i;
            added_order_inv[i] = num_init_cams;

            InitializeCameraParams(m_image_data[i], cameras[num_init_cams]);

            /* Restore the camera parameters */
            memcpy(cameras[num_init_cams].R, m_image_data[i].m_camera.m_R, 
                sizeof(double) * 9);

            matrix_transpose_product(3, 3, 3, 1, 
                m_image_data[i].m_camera.m_R,
                m_image_data[i].m_camera.m_t,
                cameras[num_init_cams].t);

            cameras[num_init_cams].t[0] *= -1.0;
            cameras[num_init_cams].t[1] *= -1.0;
            cameras[num_init_cams].t[2] *= -1.0;

            cameras[num_init_cams].f = m_image_data[i].m_camera.m_focal;
            cameras[num_init_cams].k[0] = m_image_data[i].m_camera.m_k[0];
            cameras[num_init_cams].k[1] = m_image_data[i].m_camera.m_k[1];

            /* Set the camera constraints */
            SetCameraConstraints(i, cameras + num_init_cams);

            if (m_constrain_focal) {
                /* Bad hack... */
                if (m_image_data[i].m_has_init_focal) {
                    double diff = cameras[num_init_cams].f - 
                        m_image_data[i].m_init_focal;

                    if (fabs(diff) / m_image_data[i].m_init_focal < 0.4) {
                        printf("[ReRunSFM] Constraining focal "
                            "length for camera %d\n", i);

                        /* Setup the focal length constraints */
                        SetFocalConstraint(m_image_data[i], 
                            cameras + num_init_cams);
                    }
                }
            }

            num_init_cams++;
        } else {
            added_order_inv[i] = -1;
        }
    }

    printf("\n");

    /* Set up the points, visibility mask and projections */
    printf("[ReRunSFM] Setting up views...\n");

#ifndef RERUN_ADD_POINTS
    v3_t *init_pts = new v3_t[num_pts];
#else
    v3_t *init_pts = new v3_t[m_track_data.size()];
#endif

    v3_t *colors = new v3_t[num_pts];

    for (int i = 0; i < num_pts; i++) {
        PointData &pt = m_point_data[i];
        int num_views = pt.m_views.size();

        // init_pts[i] = v3_new(pt.m_pos[0], pt.m_pos[1], -pt.m_pos[2]);
        init_pts[i] = v3_new(pt.m_pos[0], pt.m_pos[1], pt.m_pos[2]);
        colors[i] = v3_new(pt.m_color[0], pt.m_color[1], pt.m_color[2]);

        ImageKeyVector views;
        double *views_arr = new double[num_views];
        int *perm = new int[num_views];

        for (int j = 0; j < num_views; j++) {
            ImageKey ik = pt.m_views[j];

            int v = ik.first;
            int k = ik.second;

            if (m_image_data[v].m_keys[k].m_extra != -1) {
                printf("Error! Already assigned this key "
                    "[%d,%d] <- %d != %d!\n", v, k, 
                    m_image_data[v].m_keys[k].m_extra, i);
            }

            m_image_data[v].m_keys[k].m_extra = i;
            ik.first = added_order_inv[ik.first];
            views_arr[j] = (double) v;

            views.push_back(ik);
        }

        /* Sort the views */
        qsort_ascending();
        qsort_perm(num_views, views_arr, perm);

        ImageKeyVector views_sorted;
        for (int j = 0; j < num_views; j++) {
            views_sorted.push_back(views[perm[j]]);
        }

        pt_views.push_back(views_sorted);

        delete [] views_arr;
        delete [] perm;
    }

    CheckPointKeyConsistency(pt_views, added_order);

#ifdef RERUN_ADD_POINTS
    BundleAdjustAddAllNewPoints(num_pts, num_init_cams,
        added_order, cameras, init_pts, colors,
        0.0, pt_views, 16.0, 2);
#endif 

    DumpOutputFile(m_output_directory, m_bundle_output_file, 
        num_images, num_init_cams, num_pts,
        added_order, cameras, init_pts, colors, pt_views);

    RunSFM(num_pts, num_init_cams, 0, false, cameras, init_pts, added_order, 
        colors, pt_views, 0.0 /*eps2*/, S, U, V, W);

    /* Save the camera parameters and points */

    /* Cameras */

    for (int i = 0; i < num_images; i++) {
        m_image_data[i].m_camera.m_adjusted = false;
    }

    printf("Focal lengths:\n");
    for (int i = 0; i < num_init_cams; i++) {
        int img = added_order[i];

        m_image_data[img].m_camera.m_adjusted = true;
        memcpy(m_image_data[img].m_camera.m_R, cameras[i].R, 
            9 * sizeof(double));

        matrix_product(3, 3, 3, 1, 
            cameras[i].R, cameras[i].t,
            m_image_data[img].m_camera.m_t);

        matrix_scale(3, 1, 
            m_image_data[img].m_camera.m_t, -1.0, 
            m_image_data[img].m_camera.m_t);

        m_image_data[img].m_camera.m_focal = cameras[i].f;

        printf("  [%d]: %0.3f\n", img, cameras[i].f);

        m_image_data[img].m_camera.Finalize();
    }

    fflush(stdout);

    /* Points */
    m_point_data.clear();
    for (int i = 0; i < num_pts; i++) {
        /* Check if the point is visible in any view */
        if ((int) pt_views[i].size() == 0) 
            continue; /* Invisible */

        PointData pdata;
        pdata.m_pos[0] = Vx(init_pts[i]);
        pdata.m_pos[1] = Vy(init_pts[i]);
        pdata.m_pos[2] = Vz(init_pts[i]);

        pdata.m_color[0] = (float) Vx(colors[i]);
        pdata.m_color[1] = (float) Vy(colors[i]);
        pdata.m_color[2] = (float) Vz(colors[i]);

#if 1
        for (int j = 0; j < (int) pt_views[i].size(); j++) {
            int v = pt_views[i][j].first;
            int vnew = added_order[v];
            pdata.m_views.push_back(ImageKey(vnew, pt_views[i][j].second));
        }
#else
        pdata.m_views = pt_views[i];
#endif

        m_point_data.push_back(pdata);
    }

    /* Save the output file */
    DumpPointsToPly(m_output_directory, "points_readjusted.ply", 
        num_pts, num_init_cams, init_pts, colors, cameras);

    /* Dump output */
    if (m_bundle_output_file != NULL) {
        DumpOutputFile(m_output_directory, m_bundle_output_file, 
            num_images, num_init_cams, num_pts,
            added_order, cameras, init_pts, colors, pt_views);
    }

    delete [] cameras;
    delete [] init_pts;

    delete [] added_order;
    delete [] added_order_inv;
}

static int compare_doubles(const void *d1, const void *d2)
{
    double a = *(double *) d1;
    double b = *(double *) d2;

    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}


double BundlerApp::RunSFM(int num_pts, int num_cameras, int start_camera,
                          bool fix_points, camera_params_t *init_camera_params,
                          v3_t *init_pts, int *added_order, v3_t *colors,
                          std::vector<ImageKeyVector> &pt_views, double eps2, 
                          double *S, double *U, double *V, double *W,
                          bool remove_outliers)
{
#define MIN_POINTS 20
    int num_outliers = 0;
    int total_outliers = 0;
    double dist_total = 0.0;
    int num_dists = 0;

    int *remap = new int [num_pts];
    v3_t *nz_pts = new v3_t[num_pts];

    do {
        if (num_pts - total_outliers < MIN_POINTS) {
            printf("[RunSFM] Too few points remaining, exiting!\n");
            fflush(stdout);

            dist_total = DBL_MAX;
            break;
        }

        /* Set up the vmask and projections */
        char *vmask = NULL;
        double *projections = NULL;

        int num_projections = 0;
        for (int i = 0; i < num_pts; i++) {
            num_projections += (int) pt_views[i].size();
        }

        vmask = new char[num_pts * num_cameras];
        projections = new double[2 * num_projections];

        for (int i = 0; i < num_pts * num_cameras; i++)
            vmask[i] = 0;

        int arr_idx = 0;
        int nz_count = 0;
        for (int i = 0; i < num_pts; i++) {
            int num_views = (int) pt_views[i].size();

            if (num_views > 0) {
                for (int j = 0; j < num_views; j++) {
                    int c = pt_views[i][j].first;
                    int v = added_order[c];
                    int k = pt_views[i][j].second;

                    vmask[nz_count * num_cameras + c] = 1;

                    projections[2 * arr_idx + 0] = GetKey(v,k).m_x;
                    projections[2 * arr_idx + 1] = GetKey(v,k).m_y;

                    arr_idx++;
                }

                remap[i] = nz_count;
                nz_pts[nz_count] = init_pts[i];
                nz_count++;
            } else {
                remap[i] = -1;
            }
        }

        dist_total = 0.0;
        num_dists = 0;

        bool fixed_focal = m_fixed_focal_length;
        clock_t start = clock();

        run_sfm(nz_count, num_cameras, start_camera, vmask, projections, 
            fixed_focal ? 0 : 1, 0,
            m_estimate_distortion ? 1 : 0, 1,
            init_camera_params, nz_pts, 
            (m_use_constraints || m_constrain_focal) ? 1 : 0,
            (m_use_point_constraints) ? 1 : 0,
            m_point_constraints, m_point_constraint_weight,
            fix_points ? 1 : 0, m_optimize_for_fisheye, eps2, V, S, U, W);

        clock_t end = clock();

        printf("[RunSFM] run_sfm took %0.3fs\n",
            (double) (end - start) / (double) CLOCKS_PER_SEC);

        /* Check for outliers */

        start = clock();

        std::vector<int> outliers;
        std::vector<double> reproj_errors;

        for (int i = 0; i < num_cameras; i++) {
            ImageData &data = m_image_data[added_order[i]];

            double K[9] = { init_camera_params[i].f, 0.0, 0.0, 
                0.0, init_camera_params[i].f, 0.0,
                0.0, 0.0, 1.0 };
            // double w[3] = { 0.0, 0.0, 0.0 };
            double dt[3] = { init_camera_params[i].t[0],
                init_camera_params[i].t[1],
                init_camera_params[i].t[2] };		    

            /* Compute inverse distortion parameters */
            if (m_estimate_distortion) {
                double *k = init_camera_params[i].k;
                double k_dist[6] = { 0.0, 1.0, 0.0, k[0], 0.0, k[1] };
                double w_2 = 0.5 * data.GetWidth();
                double h_2 = 0.5 * data.GetHeight();
                double max_radius = 
                    sqrt(w_2 * w_2 + h_2 * h_2) / init_camera_params[i].f;

                InvertDistortion(6, 6, 0.0, max_radius, 
                    k_dist, init_camera_params[i].k_inv);
            }

            if (data.m_known_intrinsics) {
                double *k = init_camera_params[i].k_known;
                double k_dist[8] = 
                { 0.0, 1.0, 0.0, k[0], 0.0, k[1], 0.0, k[4]};
                double w_2 = 0.5 * data.GetWidth();
                double h_2 = 0.5 * data.GetHeight();
                double max_radius = 
                    sqrt(w_2 * w_2 + h_2 * h_2) /
                    init_camera_params[i].K_known[0];

                InvertDistortion(8, 6, 0.0, max_radius, k_dist, 
                    init_camera_params[i].k_inv);
            }

            int num_keys = GetNumKeys(added_order[i]);

            int num_pts_proj = 0;
            for (int j = 0; j < num_keys; j++) {
                if (GetKey(added_order[i], j).m_extra >= 0) {
                    num_pts_proj++;
                }
            }

            double *dists = new double[num_pts_proj];
            int pt_count = 0;

            std::vector<Keypoint>::iterator iter;
            // for (int j = 0; j < num_keys; j++) {
            for (iter = m_image_data[added_order[i]].m_keys.begin();
                iter != m_image_data[added_order[i]].m_keys.end();
                iter++) {

                    const Keypoint &key = *iter;

                    if (key.m_extra >= 0) {
                        double b[3], pr[2];
                        double dx, dy, dist;
                        int pt_idx = key.m_extra;

                        b[0] = Vx(nz_pts[remap[pt_idx]]);
                        b[1] = Vy(nz_pts[remap[pt_idx]]);
                        b[2] = Vz(nz_pts[remap[pt_idx]]);

                        sfm_project_rd(&(init_camera_params[i]), K, 
                            init_camera_params[i].k,
                            init_camera_params[i].R, dt, b, pr, 
                            m_estimate_distortion, true);

                        if (m_optimize_for_fisheye) {
                            /* Distort the points */
                            double x = pr[0], y = pr[1];
                            m_image_data[added_order[i]].
                                DistortPoint(x, y, pr[0], pr[1]);
                        }

                        dx = pr[0] - key.m_x;
                        dy = pr[1] - key.m_y;

                        dist = sqrt(dx * dx + dy * dy);
                        dist_total += dist;
                        num_dists++;

                        dists[pt_count] = dist;

                        pt_count++;
                    }
            }

            /* Estimate the median of the distances */
            double med = kth_element_copy(num_pts_proj, 
                iround(0.8 /* 0.9 */ * num_pts_proj),
                dists);

            median_copy(num_pts_proj, dists);

#define NUM_STDDEV 2.0 // 3.0 // 6.0
            double thresh = 1.2 * NUM_STDDEV * med; /* k * stddev */
            thresh = CLAMP(thresh, m_min_proj_error_threshold, 
                m_max_proj_error_threshold);  

            /* Compute the average reprojection error for this
            * camera */

            double sum = 0.0;
            for (int j = 0; j < num_pts_proj; j++) {
                sum += dists[j];
            }

            double avg = sum / num_pts_proj;
            printf("[RunSFM] Mean error cam %d[%d] [%d pts]: %0.3e "
                "[med: %0.3e, %0.3e]\n",
                i, added_order[i], num_pts_proj, avg, 
                kth_element_copy(num_pts_proj, 
                iround(0.5 * num_pts_proj), dists), 
                thresh);

            // printf("Outlier threshold is %0.3f\n", thresh);

            pt_count = 0;
            for (int j = 0; j < num_keys; j++) {
                int pt_idx = GetKey(added_order[i],j).m_extra;

                if (pt_idx < 0)
                    continue;

                /* Don't remove constrained points */
                if (m_use_point_constraints && 
                    Vx(m_point_constraints[pt_idx]) != 0.0) {
                        pt_count++;
                        continue;
                }

                if (dists[pt_count] > thresh) {
                    /* Remove this point from consideration */
                    bool found = false;
                    for (int k = 0; k < (int) outliers.size(); k++) {
                        if (outliers[k] == pt_idx) {
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        outliers.push_back(pt_idx);
                        reproj_errors.push_back(dists[pt_count]);
                    }
                }
                pt_count++;
            }

#define OUTPUT_VERBOSE_STATS
#ifdef OUTPUT_VERBOSE_STATS
#define NUM_ERROR_BINS 10
            qsort(dists, num_pts_proj, sizeof(double), compare_doubles);

            double pr_min = dists[0];
            double pr_max = dists[num_pts_proj-1];
            double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

            /* Break histogram into 10 bins */
            int idx_count = 0;
            for (int i = 0; i < NUM_ERROR_BINS; i++) {
                double max = pr_min + (i+1) * pr_step;
                int start = idx_count;

                while (idx_count < num_pts_proj && dists[idx_count] <= max)
                    idx_count++;

                int bin_size = idx_count - start;
                printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n", 
                    max - pr_step, max, bin_size, 
                    bin_size / (double) num_pts_proj);
            }
#endif

            delete [] dists;
        }

        /* Remove outlying points */
        if (remove_outliers) {
            for (int i = 0; i < (int) outliers.size(); i++) {
                int idx = outliers[i];

                printf("[RunSFM] Removing outlier %d "
                    "(reproj error: %0.3f)\n", idx, reproj_errors[i]);

                if (colors != NULL) {
                    Vx(colors[idx]) = 0x0;
                    Vy(colors[idx]) = 0x0;
                    Vz(colors[idx]) = 0xff;
                }

                int num_views = (int) pt_views[idx].size();

                for (int j = 0; j < num_views; j++) {
                    int v = pt_views[idx][j].first;
                    int k = pt_views[idx][j].second;

                    vmask[idx * num_cameras + v] = 0;

                    /* Sanity check */
                    if (GetKey(added_order[v], k).m_extra != idx)
                        printf("Error!  Entry for (%d,%d) "
                        "should be %d, but is %d\n",
                        added_order[v], k,
                        idx, GetKey(added_order[v], k).m_extra);

                    GetKey(added_order[v], k).m_extra = -2;
                }

                pt_views[idx].clear();
            }

            num_outliers = outliers.size();
            total_outliers += num_outliers;

            end = clock();
            printf("[RunSFM] outlier removal took %0.3fs\n",
                (double) (end - start) / (double) CLOCKS_PER_SEC);

            printf("[RunSFM] Removing %d outliers\n", num_outliers);
        }

        delete [] vmask;
        delete [] projections;

        for (int i = 0; i < num_pts; i++) {
            if (remap[i] != -1) {
                init_pts[i] = nz_pts[remap[i]];
            }
        }

        if (!remove_outliers) break;

    } while (num_outliers > 0);

    delete [] remap;
    delete [] nz_pts;

    return dist_total / num_dists;
}

void BundlerApp::ClearCameraConstraints(camera_params_t *params) 
{
    for (int i = 0; i < NUM_CAMERA_PARAMS; i++) {
        params->constrained[i] = false;
        params->constraints[i] = 0.0;
        params->weights[i] = 0.0;
    }
}

void BundlerApp::SetCameraConstraints(int cam_idx, camera_params_t *params) 
{
    const CameraInfo &cam = m_image_data[cam_idx].m_camera;

    params->constrained[0] = cam.m_constrained[0];
    params->constrained[1] = cam.m_constrained[1];
    params->constrained[2] = cam.m_constrained[2];
    params->constrained[3] = cam.m_constrained[3];
    params->constrained[4] = cam.m_constrained[4];
    params->constrained[5] = cam.m_constrained[5];
    params->constrained[6] = cam.m_constrained[6];

    if (m_estimate_distortion) {
        params->constrained[7] = true;
        params->constrained[8] = true;
    } else {
        params->constrained[7] = false;
        params->constrained[8] = false;
    }

    params->constraints[0] = cam.m_constraints[0];
    params->constraints[1] = cam.m_constraints[1];
    params->constraints[2] = cam.m_constraints[2];
    params->constraints[3] = cam.m_constraints[3];
    params->constraints[4] = cam.m_constraints[4];
    params->constraints[5] = cam.m_constraints[5];
    params->constraints[6] = cam.m_constraints[6];
    params->constraints[7] = 0.0;
    params->constraints[8] = 0.0;

    params->weights[0] = cam.m_constraint_weights[0];
    params->weights[1] = cam.m_constraint_weights[1];
    params->weights[2] = cam.m_constraint_weights[2];
    params->weights[3] = cam.m_constraint_weights[3];
    params->weights[4] = cam.m_constraint_weights[4];
    params->weights[5] = cam.m_constraint_weights[5];
    params->weights[6] = cam.m_constraint_weights[6];

    if (m_estimate_distortion) {
        params->weights[7] = m_distortion_weight;
        params->weights[8] = m_distortion_weight;
    } else {
        params->weights[7] = 0.0;
        params->weights[8] = 0.0;
    }
}

void BundlerApp::SetFocalConstraint(const ImageData &data, 
                                    camera_params_t *params)
{
    if (data.m_has_init_focal) {
        params->constrained[6] = true;
        params->constraints[6] = data.m_init_focal;
        params->weights[6] = m_constrain_focal_weight;
    }
}

/* Initialize the bundle adjustment procedure (loading an existing
* model if one exists) */
void BundlerApp::InitializeBundleAdjust(int &num_init_cams,
                                        int *added_order,
                                        int *added_order_inv,
                                        camera_params_t *cameras,
                                        v3_t *points, v3_t *colors,
                                        std::vector<ImageKeyVector> &pt_views,
                                        bool use_constraints)
{
    int num_images = GetNumImages();

    /* Initialize all keypoints to have not been matched */
    for (int i = 0; i < num_images; i++) {
        std::vector<Keypoint>::iterator iter;
        for (iter = m_image_data[i].m_keys.begin(); 
            iter != m_image_data[i].m_keys.end(); 
            iter++) {

                iter->m_extra = -1;
        }
    }

    /* Initialize the bundle adjustment with the existing model (if
    * there is one) */
    /* Cameras */
    num_init_cams = 0;

    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted) {
            printf("[InitializeBundleAdjust] Loading keys for "
                "image %d\n", i);
            m_image_data[i].LoadKeys(false, !m_optimize_for_fisheye);
            m_image_data[i].ReadKeyColors();
            SetTracks(i);

            added_order[num_init_cams] = i;
            added_order_inv[i] = num_init_cams;

            InitializeCameraParams(m_image_data[i], cameras[num_init_cams]);

            /* Restore the camera parameters */
            memcpy(cameras[num_init_cams].R, m_image_data[i].m_camera.m_R, 
                sizeof(double) * 9);

            matrix_transpose_product(3, 3, 3, 1, 
                m_image_data[i].m_camera.m_R,
                m_image_data[i].m_camera.m_t,
                cameras[num_init_cams].t);

            cameras[num_init_cams].t[0] *= -1.0;
            cameras[num_init_cams].t[1] *= -1.0;
            cameras[num_init_cams].t[2] *= -1.0;

            cameras[num_init_cams].f = m_image_data[i].m_camera.m_focal;

            if (m_estimate_distortion) {
                double k[2];

                k[0] = cameras[num_init_cams].k[0] = 
                    m_image_data[i].m_camera.m_k[0];
                k[1] = cameras[num_init_cams].k[1] = 
                    m_image_data[i].m_camera.m_k[1];

                double k_dist[6] = { 0.0, 1.0, 0.0, k[0], 0.0, k[1] };
                double w_2 = 0.5 * m_image_data[i].GetWidth();
                double h_2 = 0.5 * m_image_data[i].GetHeight();
                double max_radius = 
                    sqrt(w_2 * w_2 + h_2 * h_2) / cameras[num_init_cams].f;

                InvertDistortion(6, 6, 0.0, max_radius, 
                    k_dist, cameras[num_init_cams].k_inv);
            }

            if (m_image_data[i].m_known_intrinsics) {
                memcpy(cameras[num_init_cams].K_known,
                    m_image_data[i].m_K, 9 * sizeof(double));
                memcpy(cameras[num_init_cams].k_known,
                    m_image_data[i].m_k, 5 * sizeof(double));

                double *k = cameras[num_init_cams].k_known;
                double k_dist[8] = 
                { 0.0, 1.0, 0.0, k[0], 0.0, k[1], 0.0, k[4]};
                double w_2 = 0.5 * m_image_data[i].GetWidth();
                double h_2 = 0.5 * m_image_data[i].GetHeight();
                double max_radius = 
                    sqrt(w_2 * w_2 + h_2 * h_2) / 
                    cameras[num_init_cams].K_known[0];
                // cameras[num_init_cams].f;

                InvertDistortion(8, 6, 0.0, max_radius, k_dist, 
                    cameras[num_init_cams].k_inv);
            }

            ClearCameraConstraints(cameras + num_init_cams);
            if (use_constraints) {
                /* Setup the camera constraints */
                SetCameraConstraints(i, cameras + num_init_cams);
            }

            if (m_constrain_focal) {
                /* Setup the focal length constraints */
                /* Bad hack... */
                if (m_image_data[i].m_has_init_focal) {
                    double diff = cameras[num_init_cams].f - 
                        m_image_data[i].m_init_focal;

                    if (fabs(diff) / m_image_data[i].m_init_focal < 0.1) {
                        printf("[ReRunSFM] Constraining focal "
                            "length for camera %d\n", i);

                        /* Setup the focal length constraints */
                        SetFocalConstraint(m_image_data[i], 
                            cameras + num_init_cams);
                    }
                }
            }

            num_init_cams++;
        } else {
            if (added_order_inv != NULL)
                added_order_inv[i] = -1;
        }
    }

    /* Points */
    int n = 0;
    double error = 0.0;
    for (int i = 0; i < (int) m_point_data.size(); i++) {
        points[i] = v3_new(m_point_data[i].m_pos[0],
            m_point_data[i].m_pos[1],
            m_point_data[i].m_pos[2]);
        // -m_point_data[i].m_pos[2]);

        colors[i] = v3_new(m_point_data[i].m_color[0],
            m_point_data[i].m_color[1],
            m_point_data[i].m_color[2]);

        ImageKeyVector views;
        int num_views = (int) m_point_data[i].m_views.size();
        double *views_arr = new double[num_views];
        int *perm = new int[num_views];

        for (int j = 0; j < num_views; j++) {
            ImageKey p = m_point_data[i].m_views[j];

            if (added_order_inv[p.first] == -1) {
                printf("Error: added_order_inv[%d] doesn't exist!\n", p.first);
            }

            if (p.first < 0 || p.first >= num_images) {
                printf("Error: p.first[%d] out of range\n", p.first);
            }

            if (p.second < 0 || 
                p.second >= (int) m_image_data[p.first].m_keys.size()) {

                    printf("Error: p.second[%d,%d] out of range\n", p.second,
                        (int) m_image_data[p.first].m_keys.size());  
            }

            int image_idx = p.first;
            int key_idx = p.second;

            m_image_data[p.first].m_keys[p.second].m_extra = i;
            p.first = added_order_inv[p.first];

            views_arr[j] = (double) p.first;
            views.push_back(p);

            int track = m_image_data[image_idx].m_keys[key_idx].m_track;

            if (track != -1) { // assert(track != -1);
                m_track_data[track].m_extra = i;
            }
        }

        /* Sort the views */
        qsort_ascending();
        qsort_perm(num_views, views_arr, perm);

        ImageKeyVector views_sorted;
        for (int j = 0; j < num_views; j++) {
            views_sorted.push_back(views[perm[j]]);
        }

        pt_views.push_back(views_sorted);

        delete [] views_arr;
        delete [] perm;

        // pt_views.push_back(views);
    }

    printf("  Avg. proj error [%d projections] = %0.3e\n", n,
        sqrt(error / n));
}

/* Set up the matrix of projections and the visibility mask */
void BundlerApp::SetupProjections(int num_cameras, int num_points, 
                                  int *added_order,
                                  v2_t *projections, char *vmask) 
{
    for (int i = 0; i < num_cameras * num_points; i++)
        vmask[i] = 0;

    for (int i = 0; i < num_cameras; i++) {
        int idx = added_order[i];
        int num_keys = GetNumKeys(idx);

        for (int j = 0; j < num_keys; j++) {
            int pidx = GetKey(idx,j).m_extra;
            if (pidx >= 0) {
                /* This is a good point */
                vmask[pidx * num_cameras + i] = 1;
                projections[pidx * num_cameras + i] = 
                    v2_new(GetKey(idx,j).m_x, GetKey(idx,j).m_y);
            }
        }
    }
}

int BundlerApp::FindCameraWithMostConnectivity(int num_cameras, int num_points,
                                               int *added_order,
                                               int &parent_idx, 
                                               int &max_matches)
{
    int num_images = GetNumImages();

    /* Find the current "frontier" */
    bool *frontier = new bool[num_images];

    for (int i = 0; i < num_images; i++)
        frontier[i] = false;

    for (int i = 0; i < num_images; i++) {
        /* Check if we added this image already */
        bool added = false;
        for (int j = 0; j < num_cameras; j++) {
            if (added_order[j] == i) {
                added = true;
                break;
            }
        }

        if (added)
            frontier[i] = true;	
        else
            continue;

        /* Find the adjacent nodes in the graph */
        for (int j = 0; j < num_images; j++) {
            SetMatchesFromTracks(i, j);

            MatchIndex base = GetMatchIndex(i, j);
            unsigned int num_matches = // (int) m_match_lists[base].size();
                m_matches.GetNumMatches(base); 

            if (num_matches > 32) {
                frontier[j] = true;
            }

            // m_match_lists[base].clear();
            m_matches.ClearMatch(base);
        }
    }

    max_matches = 0;

    int i_best = -1;
    double top_score = 0.0;

    parent_idx = -1;

    int *frontier_scores = new int[num_images];
    int *seen_scores = new int[num_images];
    int *parents = new int[num_images];

    for (int i = 0; i < num_images; i++) 
        frontier_scores[i] = seen_scores[i] = 0;

    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_ignore_in_bundle)
            continue;

        if (m_only_bundle_init_focal && !m_image_data[i].m_has_init_focal)
            continue;

        /* Check if we added this image already */
        bool added = false;
        for (int j = 0; j < num_cameras; j++) {
            if (added_order[j] == i) {
                added = true;
                break;
            }
        }

        if (added)
            continue;

        /* Count the number of points that have already been
        * placed in the scene */
        int *saw = new int[num_points];

        for (int j = 0; j < num_points; j++)
            saw[j] = 0;

        int num_existing_matches = 0;
        unsigned int max_matches_curr = 0;
        int parent_idx_best = -1;

        for (int j = 0; j < num_cameras; j++) {
            int camera_idx = added_order[j];
            MatchIndex base = GetMatchIndex(i, camera_idx);

            SetMatchesFromTracks(i, camera_idx);
            // int num_matches = (int) m_match_lists[base].size();
            unsigned int num_matches = m_matches.GetNumMatches(base);

            int num_existing_matches_curr = 0;

            const std::vector<KeypointMatch> &list = 
                m_matches.GetMatchList(base);

            for (unsigned int k = 0; k < num_matches; k++) {
                // int idx2 = m_match_lists[base][k].m_idx2;
                int idx2 = list[k].m_idx2;

                if (GetKey(camera_idx, idx2).m_extra >= 0) {
                    int pidx = GetKey(camera_idx, idx2).m_extra;

                    if (saw[pidx] == 0) {
                        num_existing_matches++;
                        saw[pidx] = 1;
                    }

                    num_existing_matches_curr++;
                }
            }

            if (num_matches > max_matches_curr) {
                parent_idx_best = j;
                max_matches_curr = num_matches;
            }

            // m_match_lists[base].clear();
            m_matches.ClearMatch(base);
        }

        if (num_existing_matches > 0)
            printf("  existing_matches[%d] = %d\n", i, num_existing_matches);

        /* Now, count how many images could potentially be added to
        * the frontier after this image */
        for (int j = 0; j < num_images; j++) {
            if (frontier[j]) continue;

            SetMatchesFromTracks(i, j);

            MatchIndex base = GetMatchIndex(i, j);
            // int num_matches = (int) m_match_lists[base].size();
            unsigned int num_matches = m_matches.GetNumMatches(base);

            if (num_matches > 32)
                frontier_scores[i]++;

            // m_match_lists[base].clear();
        }

        seen_scores[i] = num_existing_matches;
        parents[i] = parent_idx_best;

        delete [] saw;
    }

    /* Find the image with the max seen_score */
    int max_seen_score = 0;

    for (int i = 0; i < num_images; i++) {
        if (seen_scores[i] > max_seen_score) 
            max_seen_score = seen_scores[i];	
    }

    if (max_seen_score == 0) {
        max_matches = 0;
        return -1;
    }

    /* Now, find the image with the highest score */
    for (int i = 0; i < num_images; i++) {
        if (seen_scores[i] == 0)
            continue;

        printf("  score[%d] = %d / %d\n", 
            i, seen_scores[i], frontier_scores[i]);
        fflush(stdout);

        /* Only accept images with at least 20% of the current max
        * matched image */
        if (seen_scores[i] < 0.20 * max_seen_score || seen_scores[i] < 32)
            continue;

        double score = frontier_scores[i];

        if (i_best == -1) {
            i_best = i;
            top_score = score;
        } else if (score == top_score) {
            if (seen_scores[i] > seen_scores[i_best]) {
                i_best = i;
                top_score = score;
            }
        } else if (score > top_score) {
            i_best = i;
            top_score = score;
        }
    }

    if (i_best == -1) {
        parent_idx = -1;
        max_matches = 0;

        delete [] frontier;
        delete [] frontier_scores;
        delete [] seen_scores;
        delete [] parents;

        return -1;
    }

    parent_idx = parents[i_best];
    max_matches = seen_scores[i_best];

    if (parent_idx == -1) {
        printf("Error: parent not found\n");
    }

    printf("  accepting image %d (%d / %d)\n", i_best, 
        seen_scores[i_best], frontier_scores[i_best]);
    fflush(stdout);

    delete [] frontier;
    delete [] frontier_scores;
    delete [] seen_scores;
    delete [] parents;

    return i_best;    
}

/* Find the camera with the most matches to existing points */
int BundlerApp::FindCameraWithMostMatches(int num_cameras, int num_points,
                                          int *added_order,
                                          int &parent_idx, int &max_matches,
                                          const std::vector<ImageKeyVector> &pt_views) 
{
    max_matches = 0;

    int i_best = -1;
    double top_score = 0.0;

    parent_idx = -1;

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_ignore_in_bundle)
            continue;

        if (m_only_bundle_init_focal && !m_image_data[i].m_has_init_focal)
            continue;

        /* Check if we added this image already */
        bool added = false;
        for (int j = 0; j < num_cameras; j++) {
            if (added_order[j] == i) {
                added = true;
                break;
            }
        }

        if (added)
            continue;

        int num_existing_matches = 0;
        int parent_idx_best = -1;

        /* Find the tracks seen by this image */
        const std::vector<int> &tracks = m_image_data[i].m_visible_points;
        int num_tracks = (int) tracks.size();

        for (int j = 0; j < num_tracks; j++) {
            int tr = tracks[j];
            if (m_track_data[tr].m_extra < 0)
                continue;

            /* This tracks corresponds to a point */
            int pt = m_track_data[tr].m_extra;
            if ((int) pt_views[pt].size() == 0)
                continue;

            num_existing_matches++;
        }

        if (num_existing_matches > 0)
            printf("  existing_matches[%d] = %d\n", i, num_existing_matches);

        double score = num_existing_matches; 

        if (score > 0.0)
            printf("  score[%d]   = %0.3f\n", i, score);

        if (score > top_score) {
            i_best = i;
            parent_idx = parent_idx_best;
            max_matches = num_existing_matches;
            top_score = score;
        }

        // delete [] saw;
    }

    if (parent_idx == -1) {
        printf("Error: parent not found\n");
    }

    return i_best;
}

/* Find all cameras with at least N matches to existing points */
std::vector<ImagePair> BundlerApp::FindCamerasWithNMatches(int n, 
                                                           int num_cameras, 
                                                           int num_points,
                                                           int *added_order,
                                                           const std::vector<ImageKeyVector> &pt_views) 
{
    std::vector<ImagePair> image_pairs;

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_ignore_in_bundle)
            continue;

        if (m_only_bundle_init_focal && !m_image_data[i].m_has_init_focal)
            continue;

        /* Check if we added this image already */
        bool added = false;
        for (int j = 0; j < num_cameras; j++) {
            if (added_order[j] == i) {
                added = true;
                break;
            }
        }

        if (added)
            continue;

        int num_existing_matches = 0;
        int parent_idx_best = -1;

        /* Find the tracks seen by this image */
        const std::vector<int> &tracks = m_image_data[i].m_visible_points;
        int num_tracks = (int) tracks.size();

        for (int j = 0; j < num_tracks; j++) {
            int tr = tracks[j];
            if (m_track_data[tr].m_extra < 0)
                continue;

            /* This tracks corresponds to a point */
            int pt = m_track_data[tr].m_extra;
            if ((int) pt_views[pt].size() == 0)
                continue;

            num_existing_matches++;
        }

        if (num_existing_matches >= n)
            image_pairs.push_back(ImagePair(i, parent_idx_best));

        // delete [] saw;
    }

    return image_pairs;
}

#define MIN_INLIERS_EST_PROJECTION 6 /* 7 */ /* 30 */ /* This constant needs
* adjustment */
#define INIT_REPROJECTION_ERROR 16.0 /* 6.0 */ /* 8.0 */

/* Pick a good initial pair of cameras to bootstrap the bundle
 * adjustment */
void BundlerApp::BundlePickInitialPair(int &i_best, int &j_best, 
                                       bool use_init_focal_only)
{
    /* Compute the match matrix */
    int num_images = GetNumImages();
    int max_matches = 0;
    double max_score = 0.0;

    int max_matches_2 = 0;
    double max_score_2 = 0.0;

    i_best = j_best = -1;
    int i_best_2 = -1;
    int j_best_2 = -1;

    int i_best_3 = -1;
    int j_best_3 = -1;

    if (m_initial_pair[0] != -1 && m_initial_pair[1] != -1) {
        i_best = m_initial_pair[0];
        j_best = m_initial_pair[1];

        printf("[BundleAdjust] Setting initial pair to "
            "%d and %d\n", i_best, j_best);

        return;
    }

    double SCORE_THRESHOLD;

    if (m_use_angular_score) {
        SCORE_THRESHOLD = 0.2;
    } else {
        SCORE_THRESHOLD = 2.0; // 1.0 // 2.0 // 1.0
    }

    /* Compute score for each image pair */
    int max_pts = 0;
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_ignore_in_bundle)
            continue;

        if (use_init_focal_only && m_use_focal_estimate && 
            !m_image_data[i].m_has_init_focal)
            continue;

        for (int j = i+1; j < num_images; j++) {
            if (m_image_data[j].m_ignore_in_bundle)
                continue;

            if (use_init_focal_only && m_use_focal_estimate && 
                !m_image_data[j].m_has_init_focal)
                continue;

            MatchIndex idx = GetMatchIndex(i, j);

            int num_matches = GetNumTrackMatches(i, j);
            max_pts += num_matches;

#define MATCH_THRESHOLD 32
#define MIN_SCORE 1.0e-1
#define MIN_MATCHES 80

            if (num_matches <= MATCH_THRESHOLD) {
                continue;
            }

            double score = 0.0;
            if (m_use_angular_score) {

            } else {
                double ratio = m_transforms[idx].m_inlier_ratio;

                if (ratio == 0.0) {
                    score = MIN_SCORE;
                } else {
                    score = 1.0 / m_transforms[idx].m_inlier_ratio;
                }
            }

            /* Compute the primary score */
            if (num_matches > max_matches && score > SCORE_THRESHOLD) {
                max_matches = num_matches;
                max_score = score;

                i_best = i;
                j_best = j;
            }

            /* Compute the backup score */
            if (num_matches > MIN_MATCHES && score > max_score_2) {
                max_matches_2 = num_matches;
                max_score_2 = score;
                i_best_2 = i;
                j_best_2 = j;
            }
        }
    }

    /* Set track pointers to -1 (GetNumTrackMatches alters these
     * values) */
    for (int i = 0; i < (int) m_track_data.size(); i++) {
        m_track_data[i].m_extra = -1;
    }

    if (i_best == -1 && j_best == -1) {
        if (i_best_2 == -1 && j_best_2 == -1) {
            printf("[BundleAdjust] Error: no good camera pairs found!\n");

            if (use_init_focal_only) {
                printf("[BundleAdjust] Trying a backup approach...\n");
                BundlePickInitialPair(i_best, j_best, false);
            } else {
                printf("[BundleAdjust] Picking first two cameras...\n");

                i_best = 0;
                j_best = 1;
            }
        } else {
            i_best = i_best_2;
            j_best = j_best_2;
        }
    }
}

/* Setup the initial camera pair for bundle adjustment */
int BundlerApp::SetupInitialCameraPair(int i_best, int j_best,
                                       double &init_focal_length_0,
                                       double &init_focal_length_1,
                                       camera_params_t *cameras,
                                       v3_t *points, v3_t *colors,
                                       std::vector<ImageKeyVector> &pt_views)
{
    /* Load the keys for the images */
    m_image_data[i_best].LoadKeys(false, !m_optimize_for_fisheye);
    m_image_data[j_best].LoadKeys(false, !m_optimize_for_fisheye);
    m_image_data[i_best].ReadKeyColors();
    m_image_data[j_best].ReadKeyColors();

    SetMatchesFromTracks(i_best, j_best);

    SetTracks(i_best);
    SetTracks(j_best);

    InitializeCameraParams(m_image_data[i_best], cameras[0]);
    InitializeCameraParams(m_image_data[j_best], cameras[1]);
    SetCameraConstraints(i_best, cameras + 0);
    SetCameraConstraints(j_best, cameras + 1);

    /* Put first camera at origin */
    cameras[0].R[0] = 1.0;  cameras[0].R[1] = 0.0;  cameras[0].R[2] = 0.0;
    cameras[0].R[3] = 0.0;  cameras[0].R[4] = 1.0;  cameras[0].R[5] = 0.0;
    cameras[0].R[6] = 0.0;  cameras[0].R[7] = 0.0;  cameras[0].R[8] = 1.0;

    /* Initialize the positions of the cameras (using constraints,
     * if provided) */
    if (m_image_data[i_best].m_camera.m_constrained[0])
        cameras[0].t[0] = m_image_data[i_best].m_camera.m_constraints[0];
    else
        cameras[0].t[0] = 0.0;

    if (m_image_data[i_best].m_camera.m_constrained[1])
        cameras[0].t[1] = m_image_data[i_best].m_camera.m_constraints[1];
    else
        cameras[0].t[1] = 0.0;

    if (m_image_data[i_best].m_camera.m_constrained[2])
        cameras[0].t[2] = m_image_data[i_best].m_camera.m_constraints[2];
    else
        cameras[0].t[2] = 0.0;

    if (m_image_data[i_best].m_has_init_focal)
        init_focal_length_0 = cameras[0].f = 
        m_image_data[i_best].m_init_focal;
    else 
        init_focal_length_0 = cameras[0].f = 
        m_init_focal_length; // INITIAL_FOCAL_LENGTH;

    if (m_image_data[j_best].m_has_init_focal)
        init_focal_length_1 = cameras[1].f = 
        m_image_data[j_best].m_init_focal;
    else
        init_focal_length_1 = cameras[1].f = 
        m_init_focal_length; // INITIAL_FOCAL_LENGTH;

    bool solved_for_extrinsics = false;
    if (m_factor_essential && m_image_data[i_best].m_has_init_focal && 
        m_image_data[j_best].m_has_init_focal && 
        !m_use_constraints) {

            /* Solve for the initial locations */
            if (EstimateRelativePose2(i_best, j_best, cameras[0], cameras[1])) {
                solved_for_extrinsics = true;
            }        
    } else {
#define INITIAL_DEPTH 3.0 // 1000.0 // 3.0 /* was 3.0, fixme! */

        /* Put second camera at origin too */
        cameras[1].R[0] = 1.0;  cameras[1].R[1] = 0.0;  cameras[1].R[2] = 0.0;
        cameras[1].R[3] = 0.0;  cameras[1].R[4] = 1.0;  cameras[1].R[5] = 0.0;
        cameras[1].R[6] = 0.0;  cameras[1].R[7] = 0.0;  cameras[1].R[8] = 1.0;

        if (m_image_data[j_best].m_camera.m_constrained[0])
            cameras[1].t[0] = m_image_data[j_best].m_camera.m_constraints[0];
        else
            cameras[1].t[0] = 0.0;

        if (m_image_data[j_best].m_camera.m_constrained[1])
            cameras[1].t[1] = m_image_data[j_best].m_camera.m_constraints[1];
        else
            cameras[1].t[1] = 0.0;

        if (m_image_data[j_best].m_camera.m_constrained[2])
            cameras[1].t[2] = m_image_data[j_best].m_camera.m_constraints[2];
        else
            cameras[1].t[2] = 0.0;        
    }

    if (m_constrain_focal) {
        SetFocalConstraint(m_image_data[i_best], cameras + 0);
        SetFocalConstraint(m_image_data[j_best], cameras + 1);
    }

    /* **** Set up the initial 3D points **** */
    printf("[BundleAdjust] Adding initial matches...\n");

    int pt_count = 0;
    MatchIndex list_idx = GetMatchIndex(i_best, j_best);
    std::vector<KeypointMatch> &list = m_matches.GetMatchList(list_idx);

    // int num_matches = (int) m_match_lists[list_idx].size();

    unsigned int num_matches = list.size();

    for (unsigned int i = 0; i < num_matches; i++) {
        int key_idx1 = list[i].m_idx1;
        int key_idx2 = list[i].m_idx2;


        printf("  Adding match %d ==> %d [%d]\n", 
            key_idx1, key_idx2, pt_count);

        double x_proj = GetKey(i_best,key_idx1).m_x;
        double y_proj = GetKey(i_best,key_idx1).m_y;

        /* Back project the point to a constant depth */
        if (!solved_for_extrinsics) {
            double x_pt = (x_proj / m_init_focal_length) * INITIAL_DEPTH;
            double y_pt = (y_proj / m_init_focal_length) * INITIAL_DEPTH;
            double z_pt = INITIAL_DEPTH + cameras[0].t[2];

            points[pt_count] = v3_new(x_pt, y_pt, z_pt);
        } else {
            double x_proj1 = GetKey(i_best,key_idx1).m_x;
            double y_proj1 = GetKey(i_best,key_idx1).m_y;
            double x_proj2 = GetKey(j_best,key_idx2).m_x;
            double y_proj2 = GetKey(j_best,key_idx2).m_y;

            double error;

            v2_t p = v2_new(x_proj1, y_proj1);
            v2_t q = v2_new(x_proj2, y_proj2);

            bool in_front = true;
            double angle = 0.0;
            points[pt_count] = Triangulate(p, q, cameras[0], cameras[1], 
                error, in_front, angle, true);

            printf(" tri.error[%d] = %0.3f\n", i, error);

            if (error > /*4.0*/ m_projection_estimation_threshold) {
                printf(" skipping point\n");
                continue;
            }
        }

        /* Get the color of the point */
        unsigned char r = GetKey(i_best,key_idx1).m_r;
        unsigned char g = GetKey(i_best,key_idx1).m_g;
        unsigned char b = GetKey(i_best,key_idx1).m_b;
        colors[pt_count] = v3_new((double) r, (double) g, (double) b);

        GetKey(i_best,key_idx1).m_extra = pt_count;
        GetKey(j_best,key_idx2).m_extra = pt_count;

        int track_idx = GetKey(i_best,key_idx1).m_track;
        m_track_data[track_idx].m_extra = pt_count;

        ImageKeyVector views;
        views.push_back(ImageKey(0, key_idx1));
        views.push_back(ImageKey(1, key_idx2));
        pt_views.push_back(views);

        pt_count++;
    }

    // m_match_lists[list_idx].clear();

    m_matches.ClearMatch(list_idx);

    return pt_count;
}


void BundlerApp::EstimateIgnoredCameras(int &curr_num_cameras,
                                        camera_params_t *cameras,
                                        int *added_order,
                                        int &curr_num_pts,
                                        v3_t *points, 
                                        v3_t *colors,
                                        std::vector<ImageKeyVector> &pt_views)
{
    int num_images = GetNumImages();

    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_ignore_in_bundle) {
            printf("[BundleAdjust] Adding camera %d\n", i);
            added_order[curr_num_cameras] = i;

            bool init_success;
            camera_params_t camera = 
                BundleInitializeImage(m_image_data[i], i,
                curr_num_cameras,
                curr_num_cameras, curr_num_pts,
                added_order, points, NULL,
                cameras, pt_views, &init_success, false);

            if (init_success) {
                cameras[curr_num_cameras] = camera;
                curr_num_cameras++;
            }
        }
    }

#if 1
    /* TEST */
    RunSFM(curr_num_pts, curr_num_cameras, 0, true,
        cameras, points, added_order, colors, pt_views, 1.0e-20 /*eps*/);
#endif

    int pt_count = 
        BundleAdjustAddAllNewPoints(curr_num_pts, curr_num_cameras,
        added_order, cameras, points, colors,
        0.0, pt_views, 16.0, 2);

    curr_num_pts = pt_count;

    /* Do another round of pose estimation to grab any remaining views */
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_ignore_in_bundle && 
            !m_image_data[i].m_camera.m_adjusted) {
                printf("[BundleAdjust] Adding camera %d\n", i);
                added_order[curr_num_cameras] = i;

                bool init_success;
                camera_params_t camera = 
                    BundleInitializeImage(m_image_data[i], i,
                    curr_num_cameras,
                    curr_num_cameras, curr_num_pts,
                    added_order, points, NULL,
                    cameras, pt_views, &init_success, false);

                if (init_success) {
                    cameras[curr_num_cameras] = camera;
                    curr_num_cameras++;
                }
        }
    }

    /* Dump output before bundle adjustment */
    if (m_bundle_output_file != NULL) {
        char buf[256];
        sprintf(buf, "%sadd.out", m_bundle_output_base);
        DumpOutputFile(m_output_directory, buf, 
            num_images, curr_num_cameras, pt_count,
            added_order, cameras, points, colors, pt_views);

#if 0
        if (m_estimate_distortion) {
            sprintf(buf, "%sadd.rd.out", m_bundle_output_base);
            DumpOutputFile(m_output_directory, buf, num_images, 
                curr_num_cameras, pt_count,
                added_order, cameras, points, colors, pt_views,
                true);
        }
#endif
    }

#if 1
    /* TEST */
    RunSFM(curr_num_pts, curr_num_cameras, 0, true,
        cameras, points, added_order, colors, pt_views, 1.0e-20);
#endif

    pt_count = 
        BundleAdjustAddAllNewPoints(curr_num_pts, curr_num_cameras,
        added_order, cameras, points, colors,
        0.0, pt_views, 16.0, 2);

#if 0
    /* Do one more round of polishing for each camera */
    for (int i = 0; i < curr_num_cameras; i++) {
        clock_t start = clock();

        int image = added_order[i];
        int num_keys = GetNumKeys(image);
        int num_visible = 0;

        printf("[EstimateIgnoredCameras] Polishing image %d (%d)\n", image, i);

        v3_t *points_visible = new v3_t[num_keys];
        v2_t *projs = new v2_t[num_keys];
        int *pt_idxs = new int[num_keys];

        int count = 0;
        for (int j = 0; j < num_keys; j++) {
            if (m_image_data[image].m_keys[j].m_extra < 0)
                continue;

            const Keypoint &key = m_image_data[image].m_keys[j];

            int track = key.m_extra;
            points_visible[count] = points_visible[track];
            projs[count] = v2_new(key.m_x, key.m_y);
            pt_idxs[count] = track;

            count++;
        }

        RefineCameraAndPoints(m_image_data[image], num_visible,
            points_visible, projs, pt_idxs, 
            cameras, added_order, pt_views,
            cameras + i, false);

        delete [] points_visible;
        delete [] projs;
        delete [] pt_idxs;

        clock_t end = clock();

        printf("[BundleAdjust] Polishing took %0.3fs\n",
            (double) (end - start) / CLOCKS_PER_SEC);
    }
#endif

    /* Dump output before bundle adjustment */
    if (m_bundle_output_file != NULL) {
        char buf[256];
        sprintf(buf, "%spre.out", m_bundle_output_base);
        DumpOutputFile(m_output_directory, buf, 
            num_images, curr_num_cameras, pt_count,
            added_order, cameras, points, colors, pt_views);

#if 0
        if (m_estimate_distortion) {
            sprintf(buf, "%spre.rd.out", m_bundle_output_base);
            DumpOutputFile(m_output_directory, buf, num_images, 
                curr_num_cameras, pt_count,
                added_order, cameras, points, colors, pt_views,
                true);
        }
#endif
    }

    /* Run a final bundle */
    RunSFM(pt_count, curr_num_cameras, 0, false,
        cameras, points, added_order, colors, pt_views, 1.0e-20);

    curr_num_pts = pt_count;

    printf("  focal lengths:\n");

    for (int i = 0; i < curr_num_cameras; i++) {
        if(m_image_data[added_order[i]].m_has_init_focal) {
            printf("   [%03d] %0.3f (%0.3f) %s %d; %0.3e, %0.3e\n", 
                i, cameras[i].f, 
                m_image_data[added_order[i]].m_init_focal,
                m_image_data[added_order[i]].m_name,
                added_order[i], cameras[i].k[0], cameras[i].k[1]); 
        } else {
            printf("   [%03d] %0.3f %s %d; %0.3e %0.3e\n", 
                i, cameras[i].f, m_image_data[added_order[i]].m_name,
                added_order[i], cameras[i].k[0], cameras[i].k[1]);
        }
    }

    DumpPointsToPly(m_output_directory, "final.ply", 
        curr_num_pts, curr_num_cameras, points, colors, cameras);
}

/* Compute pose of all cameras */
void BundlerApp::BundleAdjust() 
{
    clock_t start = clock();

    /* Compute initial image information */
    ComputeGeometricConstraints();

#if 0
    ComputeImageGraphLayout();

    CreateImageGraph();
    CreateWorkingImageGraphLargestComponent();

    std::vector<int> interior;
    ImageGraph graph = ComputeMSTWorkingGraph(interior);
    // WriteGraphMETIS(graph, "graph.met");
    PartitionGraph(graph, interior);

    OutputImageGraph("graph.partition.dot", false);
#endif

    /* For now, assume all images form one connected component */
    int num_images = GetNumImages();
    int *added_order = new int[num_images];
    int *added_order_inv = new int[num_images];

    /* Set track pointers to -1 */
    for (int i = 0; i < (int) m_track_data.size(); i++) {
        m_track_data[i].m_extra = -1;
    }

    /* **** Run bundle adjustment! **** */

    camera_params_t *cameras = new camera_params_t[num_images];

    int max_pts = (int) m_track_data.size();
    v3_t *points = new v3_t[max_pts];
    v3_t *colors = new v3_t[max_pts];
    std::vector<ImageKeyVector> pt_views;

    // int good_pair_1 = -1;
    // int good_pair_2 = -1;

    /* Initialize the bundle adjustment */
    int num_init_cams = 0;
    InitializeBundleAdjust(num_init_cams, added_order, added_order_inv,
        cameras, points, colors, pt_views, 
        m_use_constraints);

    int i_best = -1, j_best = -1, max_matches = 0;
    double max_score = 0.0;
    int curr_num_cameras, curr_num_pts;
    int pt_count;

    if (num_init_cams == 0) {
        BundlePickInitialPair(i_best, j_best, true);

        added_order[0] = i_best;
        added_order[1] = j_best;

        // good_pair_1 = 0; // i_best;
        // good_pair_2 = 1; // j_best;

        printf("[BundleAdjust] Adjusting cameras "
            "%d and %d (score = %0.3f)\n", 
            i_best, j_best, max_score);

        /* **** Set up the initial cameras **** */
        double init_focal_length_0 = 0.0, init_focal_length_1 = 0.0;
        pt_count = curr_num_pts = 
            SetupInitialCameraPair(i_best, j_best, 
            init_focal_length_0, init_focal_length_1,
            cameras, points, colors, pt_views);

        DumpOutputFile(m_output_directory, "bundle.init.out",
            num_images, 2, curr_num_pts,
            added_order, cameras, points, colors, pt_views);

        int cnp = GetNumCameraParameters();

        /* Run sfm for the first time */
        double *S = new double[2 * 2 * cnp * cnp];
        double error0;
        error0 = RunSFM(curr_num_pts, 2, 0, false,
            cameras, points, added_order, colors, pt_views, 0.0,
            S, NULL, NULL, NULL, !m_fix_necker);

        delete [] S;

        printf("  focal lengths: %0.3f, %0.3f\n", cameras[0].f, cameras[1].f);

        if (m_fix_necker) {
            /* Swap the cameras and flip the depths to deal with Necker
            * reversal */

            camera_params_t cameras_old[2];
            v3_t *points_old;

            points_old = new v3_t [curr_num_pts];

            memcpy(points_old, points, sizeof(v3_t) * curr_num_pts);
            memcpy(cameras_old, cameras, sizeof(camera_params_t) * 2);

            camera_params_t tmp = cameras[0];
            memcpy(cameras[0].R, cameras[1].R, sizeof(double) * 9);
            memcpy(cameras[0].t, cameras[1].t, sizeof(double) * 3);
            cameras[0].f = init_focal_length_0;
            cameras[0].k[0] = cameras[0].k[1] = 0.0;

            memcpy(cameras[1].R, tmp.R, sizeof(double) * 9);
            memcpy(cameras[1].t, tmp.t, sizeof(double) * 3);
            cameras[1].f = init_focal_length_1;
            cameras[1].k[0] = cameras[1].k[1] = 0.0;

            double K1inv[9] = 
            { 1.0 / cameras[0].f, 0.0, 0.0,
            0.0, 1.0 / cameras[0].f, 0.0,
            0.0, 0.0, 1.0 };

            double K2inv[9] = 
            { 1.0 / cameras[1].f, 0.0, 0.0,
            0.0, 1.0 / cameras[1].f, 0.0,
            0.0, 0.0, 1.0 };

            for (int i = 0; i < curr_num_pts; i++) {
                int k1 = pt_views[i][0].second;
                int k2 = pt_views[i][1].second;

                double proj1[3] = { GetKey(added_order[0],k1).m_x,
                    GetKey(added_order[0],k1).m_y,
                    1.0 };

                double proj2[3] = { GetKey(added_order[1],k2).m_x,
                    GetKey(added_order[1],k2).m_y,
                    1.0 };

                if (m_optimize_for_fisheye) {
                    double x1 = proj1[0];
                    double y1 = proj1[1];

                    double x2 = proj2[0];
                    double y2 = proj2[1];

                    m_image_data[added_order[0]].
                        UndistortPoint(x1, y1, proj1[0], proj1[1]);
                    m_image_data[added_order[1]].
                        UndistortPoint(x2, y2, proj2[0], proj2[1]);
                }


                double proj1_norm[3], proj2_norm[3];

                matrix_product(3, 3, 3, 1, K1inv, proj1, proj1_norm);
                matrix_product(3, 3, 3, 1, K2inv, proj2, proj2_norm);

                v2_t p = v2_new(proj1_norm[0] / proj1_norm[2],
                    proj1_norm[1] / proj1_norm[2]);

                v2_t q = v2_new(proj2_norm[0] / proj2_norm[2],
                    proj2_norm[1] / proj2_norm[2]);

                double proj_error;

                double t1[3], t2[3];

                /* Put the translation in standard form */
                matrix_product(3, 3, 3, 1, cameras[0].R, cameras[0].t, t1);
                matrix_scale(3, 1, t1, -1.0, t1);
                matrix_product(3, 3, 3, 1, cameras[1].R, cameras[1].t, t2);
                matrix_scale(3, 1, t2, -1.0, t2);

                points[i] = triangulate(p, q, 
                    cameras[0].R, t1, cameras[1].R, t2,
                    &proj_error);
            }

            double error1;
            error1 = RunSFM(curr_num_pts, 2, 0, false,
                cameras, points, added_order, colors, pt_views);

            printf("  focal lengths: %0.3f, %0.3f\n", 
                cameras[0].f, cameras[1].f);

#if 0
            if (error0 < error1) {
                /* Swap back */
                printf("Restoring pre-Necker configuration\n");

                memcpy(points, points_old, sizeof(v3_t) * curr_num_pts);
                memcpy(cameras, cameras_old, sizeof(camera_params_t) * 2);
            }

            delete [] points_old;
#endif
        }

        DumpPointsToPly(m_output_directory, "points001.ply", 
            curr_num_pts, 2, points, colors, cameras);

        if (m_bundle_output_base != NULL) {
            char buf[256];
            sprintf(buf, "%s%03d.out", m_bundle_output_base, 1);
            DumpOutputFile(m_output_directory, buf, num_images, 2, curr_num_pts,
                added_order, cameras, points, colors, pt_views);

#if 0
            if (m_estimate_distortion) {
                sprintf(buf, "%s%03d.rd.out", m_bundle_output_base, 1);
                DumpOutputFile(m_output_directory, buf, 
                    num_images, 2, curr_num_pts,
                    added_order, cameras, points, colors, pt_views, 
                    true);
            }
#endif
        }

        curr_num_cameras = 2;
    } else {
#if 0
        if (m_initial_pair[0] == -1 || m_initial_pair[1] == -1) {
            printf("[BundleAdjust] Error: initial good pair "
                "not provided!\n");
            printf("[BundleAdjust] Please specify a pair of "
                "cameras with medium baseline using\n"
                "  --init_pair1 <img1> and --init_pair2 <img2>\n");
            exit(1);
        }

        good_pair_1 = added_order_inv[m_initial_pair[0]];
        good_pair_2 = added_order_inv[m_initial_pair[1]];

        if (good_pair_1 == -1 || good_pair_2 == -1) {
            printf("[BundleAdjust] Error: initial pair haven't "
                "been adjusted!\n");
            printf("[BundleAdjust] Please specify another pair!\n");
            exit(0);
        }
#endif

        curr_num_cameras = num_init_cams;
        pt_count = curr_num_pts = (int) m_point_data.size();
    }

    for (int round = curr_num_cameras; 
        round < num_images; 
        round++, curr_num_cameras++) {

            int parent_idx = -1;
            int next_idx;

            if (m_construct_max_connectivity)
                next_idx = FindCameraWithMostConnectivity(round, curr_num_pts, 
                added_order, parent_idx, 
                max_matches);
            else
                next_idx = FindCameraWithMostMatches(round, curr_num_pts, 
                added_order, parent_idx, 
                max_matches, pt_views);

            printf("[BundleAdjust] max_matches = %d\n", max_matches);

            if (max_matches < 16)
                break; /* No more connections */

            /* Now, throw the new camera into the mix and redo bundle
            * adjustment */
            added_order[round] = next_idx;

            printf("[BundleAdjust[%d]] Adjusting camera %d "
                "(parent = %d, matches = %d)\n", 
                round, next_idx, 
                (parent_idx == -1 ? -1 : added_order[parent_idx]), max_matches);


            /* **** Set up the new camera **** */
#if 0
            cameras[round] = 
                BundleInitializeImage(m_image_data[next_idx], 
                next_idx, round, curr_num_pts,
                added_order, points, cameras + parent_idx,
                cameras, pt_views);
#else
            BundleInitializeImageFullBundle(next_idx, parent_idx, 
                round, curr_num_pts,
                added_order, cameras, points, colors,
                pt_views);
#endif

            /* Compute the distance between the first pair of cameras */
#if 0
            double dist0 = GetCameraDistance(cameras + good_pair_1, 
                cameras + good_pair_2, 
                m_explicit_camera_centers);
#else
            double dist0 = 0.0;
#endif

            printf("[BundleAdjust] Adding new matches\n");

            if (!m_skip_add_points) {
                pt_count = 
                    BundleAdjustAddAllNewPoints(curr_num_pts, curr_num_cameras + 1,
                    added_order, cameras, 
                    points, colors,
                    dist0, pt_views);
            }

            curr_num_pts = pt_count;
            printf("[BundleAdjust] Number of points = %d\n", pt_count);
            fflush(stdout);

            /* Run sfm again to update parameters */
            RunSFM(curr_num_pts, round + 1, 0, false,
                cameras, points, added_order, colors, pt_views);

            /* Remove bad points and cameras */
            RemoveBadPointsAndCameras(curr_num_pts, curr_num_cameras + 1, 
                added_order, cameras, points, colors, 
                pt_views);

            printf("  focal lengths:\n");

            for (int i = 0; i <= round; i++) {
                if(m_image_data[added_order[i]].m_has_init_focal) {
                    printf("   [%03d] %0.3f (%0.3f) %s %d; %0.3e, %0.3e\n", 
                        i, cameras[i].f, 
                        m_image_data[added_order[i]].m_init_focal,
                        m_image_data[added_order[i]].m_name,
                        added_order[i], cameras[i].k[0], cameras[i].k[1]); 
                } else {
                    printf("   [%03d] %0.3f %s %d; %0.3e %0.3e\n", 
                        i, cameras[i].f, m_image_data[added_order[i]].m_name,
                        added_order[i], cameras[i].k[0], cameras[i].k[1]);
                }
            }

            /* Dump output for this round */
            char buf[256];
            sprintf(buf, "points%03d.ply", round);

            DumpPointsToPly(m_output_directory, buf, 
                curr_num_pts, round+1, points, colors, cameras);

            if (m_bundle_output_base != NULL) {
                sprintf(buf, "%s%03d.out", m_bundle_output_base, round);
                DumpOutputFile(m_output_directory, buf, num_images, 
                    curr_num_cameras + 1, curr_num_pts,
                    added_order, cameras, points, colors, pt_views);

#if 0
                if (m_estimate_distortion) {
                    sprintf(buf, "%s%03d.rd.out", m_bundle_output_base, round);
                    DumpOutputFile(m_output_directory, buf, num_images, 
                        curr_num_cameras + 1, curr_num_pts,
                        added_order, cameras, points, colors, pt_views,
                        true);
                }
#endif
            }
    }

    clock_t end = clock();

    printf("[BundleAdjust] Bundle adjustment took %0.3fs\n",
        (end - start) / ((double) CLOCKS_PER_SEC));

    if (m_estimate_ignored) {
        EstimateIgnoredCameras(curr_num_cameras,
            cameras, added_order,
            curr_num_pts, points, colors, pt_views);
    }

    /* Dump output */
    if (m_bundle_output_file != NULL) {
        DumpOutputFile(m_output_directory, m_bundle_output_file, 
            num_images, curr_num_cameras, curr_num_pts,
            added_order, cameras, points, colors, pt_views);

#if 0
        if (m_estimate_distortion) {
            char buf[256];
            sprintf(buf, "%s.rd.out", m_bundle_output_file);
            DumpOutputFile(m_output_directory, buf, num_images, 
                curr_num_cameras, curr_num_pts,
                added_order, cameras, points, colors, pt_views,
                true);
        }
#endif
    }

    /* Save the camera parameters and points */

    /* Cameras */
    for (int i = 0; i < num_images; i++) {
        m_image_data[i].m_camera.m_adjusted = false;
    }

    for (int i = 0; i < curr_num_cameras; i++) {
        int img = added_order[i];

        m_image_data[img].m_camera.m_adjusted = true;
        memcpy(m_image_data[img].m_camera.m_R, cameras[i].R, 
            9 * sizeof(double));

        matrix_product(3, 3, 3, 1, 
            cameras[i].R, cameras[i].t,
            m_image_data[img].m_camera.m_t);

        matrix_scale(3, 1, 
            m_image_data[img].m_camera.m_t, -1.0, 
            m_image_data[img].m_camera.m_t);	    

        m_image_data[img].m_camera.m_focal = cameras[i].f;

        m_image_data[img].m_camera.Finalize();
    }

    /* Points */
    for (int i = 0; i < curr_num_pts; i++) {
        /* Check if the point is visible in any view */
        if ((int) pt_views[i].size() == 0) 
            continue; /* Invisible */

        PointData pdata;
        pdata.m_pos[0] = Vx(points[i]);
        pdata.m_pos[1] = Vy(points[i]);
        pdata.m_pos[2] = Vz(points[i]);

        pdata.m_color[0] = (float) Vx(colors[i]);
        pdata.m_color[1] = (float) Vy(colors[i]);
        pdata.m_color[2] = (float) Vz(colors[i]);

#if 1
        for (int j = 0; j < (int) pt_views[i].size(); j++) {
            int v = pt_views[i][j].first;
            int vnew = added_order[v];
            pdata.m_views.push_back(ImageKey(vnew, pt_views[i][j].second));
        }
#else
        pdata.m_views = pt_views[i];
#endif
        // pdata.m_views = pt_views[i];

        m_point_data.push_back(pdata);
    }

    delete [] added_order;
    delete [] added_order_inv;

    SetMatchesFromPoints();

#if 0
    bool *image_mask = new bool[num_images];

    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted)
            image_mask[i] = true;
        else 
            image_mask[i] = false;
    }
#endif
}




std::vector<int> RefineCameraParameters(const ImageData &data,
                                        int num_points, 
                                        v3_t *points, v2_t *projs, 
                                        int *pt_idxs, camera_params_t *camera,
                                        double *error_out, 
                                        bool adjust_focal,
                                        bool remove_outliers,
                                        bool optimize_for_fisheye,
                                        bool estimate_distortion,
                                        double min_proj_error_threshold,
                                        double max_proj_error_threshold)
{
    int num_points_curr = num_points;
    v3_t *points_curr = new v3_t[num_points];
    v2_t *projs_curr = new v2_t[num_points];

    memcpy(points_curr, points, num_points * sizeof(v3_t));
    memcpy(projs_curr, projs, num_points * sizeof(v2_t));

    std::vector<int> inliers;

    for (int i = 0; i < num_points; i++)
        inliers.push_back(i);

    int round = 0;

    /* First refine with the focal length fixed */
    camera_refine(num_points_curr, points_curr, projs_curr, camera, 0, 0);    

    while (1) {
        printf("[RefineCameraParameters] Calling with %d points\n", 
            num_points_curr);

        camera_refine(num_points_curr, points_curr, projs_curr, camera, 
            adjust_focal ? 1 : 0, estimate_distortion ? 1 : 0);

        if (!remove_outliers)
            break;

        v3_t *points_next = new v3_t[num_points];
        v2_t *projs_next = new v2_t[num_points];

        int count = 0;
        double error = 0.0;
        std::vector<int> inliers_next;

        double *errors = new double[num_points_curr];

        for (int i = 0; i < num_points_curr; i++) {
            v2_t pr = sfm_project_final(camera, points_curr[i], 1,
                estimate_distortion ? 1 : 0);

            if (optimize_for_fisheye) {
                /* Distort pr */

                double x = Vx(pr);
                double y = Vy(pr);
                data.DistortPoint(x, y, Vx(pr), Vy(pr));
            }

            double dx = Vx(pr) - Vx(projs_curr[i]);
            double dy = Vy(pr) - Vy(projs_curr[i]);
            double diff = sqrt(dx * dx + dy * dy);

            errors[i] = diff;
            error += diff;
        }

        printf("[RefineCameraParameters] Error: %0.3f\n", 
            error / num_points_curr);

        /* Sort and histogram errors */
        double med = kth_element_copy(num_points_curr, 
            iround(0.95 * num_points_curr),
            errors);

        /* We will tolerate any match with projection error < 8.0 */
        double threshold = 1.2 * NUM_STDDEV * med; /* k * stddev */
        threshold = CLAMP(threshold, min_proj_error_threshold, 
            max_proj_error_threshold);  

        // double threshold = min_proj_error_threshold;

        // double threshold = MAX(8.0, med);

        printf("[RefineCameraParameters] Threshold = %0.3f\n", threshold);
        for (int i = 0; i < num_points_curr; i++) {
            if (errors[i] < threshold) {
                inliers_next.push_back(inliers[i]);

                points_next[count] = points_curr[i];
                projs_next[count] = projs_curr[i];

                count++;
            } else {
                if (pt_idxs != NULL) {
                    printf("[RefineCameraParameters] Removing point [%d] with "
                        "reprojection error %0.3f\n", pt_idxs[i], 
                        errors[i]);
                } else {
                    printf("[RefineCameraParameters] Removing point with "
                        "reprojection error %0.3f\n", errors[i]);
                }
            }
        }

#if 1
        qsort(errors, num_points_curr, sizeof(double), compare_doubles);

        double pr_min = errors[0];
        double pr_max = errors[num_points_curr-1];
        double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

        /* Break histogram into 10 bins */
        int idx_count = 0;
        for (int i = 0; i < NUM_ERROR_BINS; i++) {
            double max = pr_min + (i+1) * pr_step;
            int start = idx_count;

            while (idx_count < num_points_curr && errors[idx_count] <= max)
                idx_count++;

            int bin_size = idx_count - start;
            printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n", 
                max - pr_step, max, bin_size, 
                bin_size / (double) num_points_curr);
        }
#endif

        delete [] points_curr;
        delete [] projs_curr;
        delete [] errors;

        points_curr = points_next;
        projs_curr = projs_next;

        if (count == num_points_curr)
            break;  /* We're done */

        num_points_curr = count;

        inliers = inliers_next;

        if (count == 0) /* Out of measurements */
            break;

        round++;

        if (error_out != NULL) {
            *error_out = error;
        }
    }

    printf("[RefineCameraParameters] Exiting after %d rounds "
        "with %d / %d points\n", round + 1, num_points_curr, num_points);

    delete [] points_curr;
    delete [] projs_curr;

    return inliers;
}


double BundlerApp::RefinePoints(int num_points, v3_t *points, v2_t *projs,
                                int *pt_idxs, camera_params_t *cameras,
                                int *added_order,
                                const std::vector<ImageKeyVector> &pt_views,
                                camera_params_t *camera_out)
{
    double error = 0.0;

    /* Triangulate each of the points */
    for (int i = 0; i < num_points; i++) {
        int pt_idx = pt_idxs[i];

        int num_views = (int) pt_views[pt_idx].size() + 1;

        if (num_views < 2) continue;

        v2_t *pv = new v2_t[num_views];
        double *Rs = new double[9 * num_views];
        double *ts = new double[3 * num_views];

        for (int j = 0; j < num_views; j++) {
            camera_params_t *cam = NULL;

            if (j < num_views - 1) {
                int camera_idx = pt_views[pt_idx][j].first;
                int image_idx = added_order[camera_idx];
                int key_idx = pt_views[pt_idx][j].second;
                Keypoint &key = GetKey(image_idx, key_idx);

                double p3[3] = { key.m_x, key.m_y, 1.0 };
                double K[9], Kinv[9];
                GetIntrinsics(cameras[camera_idx], K);
                matrix_invert(3, K, Kinv);

                double p_n[3];
                matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

                pv[j] = v2_new(p_n[0], p_n[1]);

                cam = cameras + camera_idx;
            } else {
                double p3[3] = { Vx(projs[i]), Vy(projs[i]), 1.0 };
                double K[9], Kinv[9];
                GetIntrinsics(*camera_out, K);
                matrix_invert(3, K, Kinv);

                double p_n[3];
                matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

                pv[j] = v2_new(p_n[0], p_n[1]);
                cam = camera_out;
            }

            memcpy(Rs + 9 * j, cam->R, 9 * sizeof(double));

            matrix_product(3, 3, 3, 1, cam->R, cam->t, ts + 3 * j);
            matrix_scale(3, 1, ts + 3 * j, -1.0, ts + 3 * j);
        }

        // points[i] = triangulate_n(num_views, pv, Rs, ts, &error_curr);
        double error_curr = 0.0;
        points[i] = triangulate_n_refine(points[i], num_views, pv, Rs, ts, 
            &error_curr);

        v2_t pr = sfm_project_final(camera_out, points[i], 1,
            m_estimate_distortion ? 1 : 0);

        double dx = Vx(pr) - Vx(projs[i]);
        double dy = Vy(pr) - Vy(projs[i]);

        error += dx * dx + dy * dy;

        delete [] pv;
        delete [] Rs;
        delete [] ts;
    }

    return sqrt(error / num_points);
}

std::vector<int> 
BundlerApp::RefineCameraAndPoints(const ImageData &data, int num_points,
                                  v3_t *points, v2_t *projs,
                                  int *pt_idxs, 
                                  camera_params_t *cameras,
                                  int *added_order,
                                  const std::vector<ImageKeyVector> &pt_views,
                                  camera_params_t *camera_out,
                                  bool remove_outliers)
{
    // double error_thresh = 1.0e-6;
    double error_old = DBL_MAX;
    double derror;

    bool removed;
    std::vector<int> inliers_out;

    for (int i = 0; i < num_points; i++)
        inliers_out.push_back(i);

    int num_points_curr = num_points;
    v3_t *points_curr = new v3_t[num_points];
    v2_t *projs_curr = new v2_t[num_points];
    int *pt_idxs_curr = new int[num_points];

    memcpy(points_curr, points, sizeof(v3_t) * num_points);
    memcpy(projs_curr, projs, sizeof(v2_t) * num_points);
    memcpy(pt_idxs_curr, pt_idxs, sizeof(int) * num_points);

    do {
        removed = false;
        // do {
        {
            double error = 0.0;

            /* Refine the camera */
            RefineCameraParameters(data, num_points_curr, 
                points_curr, projs_curr, 
                pt_idxs_curr, camera_out, &error, 
                !m_fixed_focal_length, false, // true,
                m_optimize_for_fisheye,
                m_estimate_distortion,
                m_min_proj_error_threshold,
                m_max_proj_error_threshold);


#if 1 /* Change to 0 to skip point polishing */
            /* Refine the points */
            error = RefinePoints(num_points_curr, points_curr, projs_curr, 
                pt_idxs_curr, cameras, added_order, 
                pt_views, camera_out);

            printf("[RefineCameraAndPoints] "
                "Error (after point polishing): %0.3f\n", error);

            derror = error_old - error;
            error_old = error;
#endif
        }
        // } while (derror > error_thresh);

        if (remove_outliers) {
            /* Refine the camera once more and remove outliers */
            std::vector<int> inliers;
            inliers = RefineCameraParameters(data, num_points_curr, 
                points_curr, projs_curr,
                pt_idxs_curr, camera_out, 
                NULL, !m_fixed_focal_length, true,
                m_optimize_for_fisheye,
                m_estimate_distortion,
                m_min_proj_error_threshold,
                m_max_proj_error_threshold);

            int num_inliers = inliers.size();

            std::vector<int> inliers_out_next;
            if (num_inliers < num_points_curr) {
                removed = true;

                for (int i = 0; i < num_inliers; i++) {
                    points_curr[i] = points_curr[inliers[i]];
                    projs_curr[i] = projs_curr[inliers[i]];
                    pt_idxs_curr[i] = pt_idxs_curr[inliers[i]];
                    inliers_out_next.push_back(inliers_out[i]);
                }

                num_points_curr = num_inliers;
                inliers_out = inliers_out_next;
            }
        } else {
            RefineCameraParameters(data, num_points_curr, 
                points_curr, projs_curr,
                pt_idxs_curr, camera_out, 
                NULL, !m_fixed_focal_length, false,
                m_optimize_for_fisheye,
                m_estimate_distortion,
                m_min_proj_error_threshold,
                m_max_proj_error_threshold);
        }

    } while (removed);

    delete [] points_curr;
    delete [] projs_curr;
    delete [] pt_idxs_curr;

    return inliers_out;
}


bool FindAndVerifyCamera(int num_points, v3_t *points_solve, v2_t *projs_solve,
                         int *idxs_solve,
                         double *K, double *R, double *t, 
                         double proj_estimation_threshold,
                         double proj_estimation_threshold_weak,
                         std::vector<int> &inliers,
                         std::vector<int> &inliers_weak,
                         std::vector<int> &outliers)
{
    /* First, find the projection matrix */
    double P[12];
    int r = -1;

    if (num_points >= 9) {
        r = find_projection_3x4_ransac(num_points, 
            points_solve, projs_solve, 
            P, /* 2048 */ 4096 /* 100000 */, 
            proj_estimation_threshold);
    }

    if (r == -1) {
        printf("[FindAndVerifyCamera] Couldn't find projection matrix\n");
        return false;
    }

    /* If number of inliers is too low, fail */
    if (r <= MIN_INLIERS_EST_PROJECTION) {
        printf("[FindAndVerifyCamera] Too few inliers to use "
            "projection matrix\n");
        return false;
    }

    double KRinit[9], Kinit[9], Rinit[9], tinit[3];
    memcpy(KRinit + 0, P + 0, 3 * sizeof(double));
    memcpy(KRinit + 3, P + 4, 3 * sizeof(double));
    memcpy(KRinit + 6, P + 8, 3 * sizeof(double));

    dgerqf_driver(3, 3, KRinit, Kinit, Rinit);	    

    /* We want our intrinsics to have a certain form */
    FixIntrinsics(P, Kinit, Rinit, tinit);
    matrix_scale(3, 3, Kinit, 1.0 / Kinit[8], Kinit);

    printf("[FindAndVerifyCamera] Estimated intrinsics:\n");
    matrix_print(3, 3, Kinit);
    printf("[FindAndVerifyCamera] Estimated extrinsics:\n");
    matrix_print(3, 3, Rinit);
    matrix_print(1, 3, tinit);
    fflush(stdout);

    /* Check cheirality constraint */
    printf("[FindAndVerifyCamera] Checking consistency...\n");

    double Rigid[12] = 
    { Rinit[0], Rinit[1], Rinit[2], tinit[0],
    Rinit[3], Rinit[4], Rinit[5], tinit[1],
    Rinit[6], Rinit[7], Rinit[8], tinit[2] };

    int num_behind = 0;
    for (int j = 0; j < num_points; j++) {
        double p[4] = { Vx(points_solve[j]), 
            Vy(points_solve[j]),
            Vz(points_solve[j]), 1.0 };
        double q[3], q2[3];

        matrix_product(3, 4, 4, 1, Rigid, p, q);
        matrix_product331(Kinit, q, q2);

        double pimg[2] = { -q2[0] / q2[2], -q2[1] / q2[2] };
        double diff = 
            (pimg[0] - Vx(projs_solve[j])) * 
            (pimg[0] - Vx(projs_solve[j])) + 
            (pimg[1] - Vy(projs_solve[j])) * 
            (pimg[1] - Vy(projs_solve[j]));

        diff = sqrt(diff);

        if (diff < proj_estimation_threshold)
            inliers.push_back(j);

        if (diff < proj_estimation_threshold_weak) {
            inliers_weak.push_back(j);
        } else {
            printf("[FindAndVerifyCamera] Removing point [%d] "
                "(reproj. error = %0.3f)\n", idxs_solve[j], diff);
            outliers.push_back(j);
        }

        // EDIT!!!
        if (q[2] > 0.0)
            num_behind++;  /* Cheirality constraint violated */
    }

    if (num_behind >= 0.9 * num_points) {
        printf("[FindAndVerifyCamera] Error: camera is pointing "
            "away from scene\n");
        return false;
    }

    memcpy(K, Kinit, sizeof(double) * 9);
    memcpy(R, Rinit, sizeof(double) * 9);
    memcpy(t, tinit, sizeof(double) * 3);

    // #define COLIN_HACK
#ifdef COLIN_HACK
    matrix_ident(3, R);
    t[0] = t[1] = t[2] = 0.0;
#endif

    return true;
}



camera_params_t 
BundlerApp::BundleInitializeImage(ImageData &data, 
                                  int image_idx, int camera_idx,
                                  int num_cameras, int num_points,
                                  int *added_order,
                                  v3_t *points,
                                  camera_params_t *parent,
                                  camera_params_t *cameras,
                                  std::vector<ImageKeyVector> &pt_views,
                                  bool *success_out,
                                  bool refine_cameras_and_points)
{
    clock_t start = clock();

    if (success_out != NULL)
        *success_out = true;

    /* Load the keys */
    data.LoadKeys(false, !m_optimize_for_fisheye);
    SetTracks(image_idx);

    /* **** Connect the new camera to any existing points **** */
    int num_pts_solve = 0;
    int num_keys = (int) data.m_keys.size();
    v3_t *points_solve = new v3_t[num_keys];
    v2_t *projs_solve = new v2_t[num_keys];
    v2_t *projs_solve_orig = new v2_t[num_keys];
    int *idxs_solve = new int[num_keys];
    int *keys_solve = new int[num_keys];

    printf("[BundleInitializeImage] "
        "Connecting existing matches...\n");

    /* Find the tracks seen by this image */
    std::vector<int> &tracks = data.m_visible_points;
    int num_tracks = (int) tracks.size();

    for (int i = 0; i < num_tracks; i++) {
        int tr = tracks[i];
        if (m_track_data[tr].m_extra < 0)
            continue;

        /* This tracks corresponds to a point */
        int pt = m_track_data[tr].m_extra;
        if ((int) pt_views[pt].size() == 0)
            continue;

        int key =  data.m_visible_keys[i];

        // printf("  Connecting existing point [%d] (cam: %d)\n", 
        //        pt, image_idx);

        /* Add the point to the set we'll use to solve for
        * the camera position */
        points_solve[num_pts_solve] = points[pt];

        if (m_optimize_for_fisheye) {
            double x = data.m_keys[key].m_x;
            double y = data.m_keys[key].m_y;
            double x_u, y_u;
            data.UndistortPoint(x, y, x_u, y_u);
            projs_solve[num_pts_solve] = v2_new(x_u, y_u);
            projs_solve_orig[num_pts_solve] = v2_new(x, y);
        } else {
            projs_solve[num_pts_solve] = 
                v2_new(data.m_keys[key].m_x, data.m_keys[key].m_y);
        }

        idxs_solve[num_pts_solve] = pt;
        keys_solve[num_pts_solve] = key;

        num_pts_solve++;
    }

    if (num_pts_solve < m_min_max_matches) {
        printf("[BundleInitializeImage] Couldn't initialize\n");

        if (success_out != NULL)
            *success_out = false;

        camera_params_t dummy;

        delete [] points_solve;
        delete [] projs_solve;
        delete [] projs_solve_orig;
        delete [] idxs_solve;
        delete [] keys_solve;

        m_image_data[image_idx].UnloadKeys();

        return dummy;
    }

    /* **** Solve for the camera position **** */
    printf("[BundleInitializeImage] Initializing camera...\n");
    fflush(stdout);
    double Kinit[9], Rinit[9], tinit[3];
    std::vector<int> inliers, inliers_weak, outliers;
    bool success = 
        FindAndVerifyCamera(num_pts_solve, points_solve, projs_solve,
        idxs_solve, Kinit, Rinit, tinit, 
        m_projection_estimation_threshold, 
        16.0 * m_projection_estimation_threshold, /*4.0*/
        inliers, inliers_weak, outliers);

    if (!success) {
        printf("[BundleInitializeImage] Couldn't initialize\n");

        if (success_out != NULL)
            *success_out = false;

        camera_params_t dummy;

        delete [] points_solve;
        delete [] projs_solve;
        delete [] projs_solve_orig;
        delete [] idxs_solve;
        delete [] keys_solve;

        m_image_data[image_idx].UnloadKeys();

        return dummy;
    }

    camera_params_t camera_new;

    InitializeCameraParams(data, camera_new);

    /* Start with the new camera at same place as the best match */

    if (success) {
        /* Set up the new camera */
        memcpy(camera_new.R, Rinit, 9 * sizeof(double));

        matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, camera_new.t);
        matrix_scale(3, 1, camera_new.t, -1.0, camera_new.t);

        /* Set up the new focal length */
        SetCameraConstraints(added_order[num_cameras], &camera_new);

        if (m_fixed_focal_length) {
            camera_new.f = m_init_focal_length;
        } else {
            if (m_use_focal_estimate) {
                if (data.m_has_init_focal) {
                    double ratio;
                    double init = data.m_init_focal;
                    double obs = 0.5 * (Kinit[0] + Kinit[4]);

                    printf("[BundleInitializeImage] "
                        "Camera has initial focal length of %0.3f\n", init);

                    if (init > obs) ratio = init / obs;
                    else            ratio = obs / init;

                    if (ratio < 1.4 || m_trust_focal_estimate) {
                        camera_new.f = data.m_init_focal;
                        if (m_constrain_focal)
                            SetFocalConstraint(m_image_data[image_idx], 
                            &camera_new);
                    } else {
                        printf("[BundleInitializeImage] "
                            "Estimated focal length of %0.3f "
                            "is too different\n", obs);
                        camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
                    }
                } else {
                    camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
                }
                // } else if (data.m_has_init_focal) {
                //    camera_new.f = data.m_init_focal;
            } else {
                if (parent != NULL)
                    camera_new.f = parent->f;
                else
                    camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
                // 1.2 * MAX(data.GetWidth(), data.GetHeight());
            }
        }
    } else {
        printf("[BundleInitializeImage] Error!  "
            "Pose estimation failed!\n");
    }

    /* **** Finally, start the bundle adjustment **** */
    printf("[BundleInitializeImage] Adjusting...\n");
    fflush(stdout);

    int num_inliers = (int) inliers_weak.size();

    v3_t *points_final = new v3_t[num_inliers];
    v2_t *projs_final = new v2_t[num_inliers];
    int *idxs_final = new int[num_inliers];
    int *keys_final = new int[num_inliers];
    int num_points_final = num_inliers;

    for (int i = 0; i < num_inliers; i++) {
        points_final[i] = points_solve[inliers_weak[i]];
        if (m_optimize_for_fisheye)
            projs_final[i] = projs_solve_orig[inliers_weak[i]];
        else
            projs_final[i] = projs_solve[inliers_weak[i]];

        idxs_final[i] = idxs_solve[inliers_weak[i]];
        keys_final[i] = keys_solve[inliers_weak[i]];
    }

    if (refine_cameras_and_points) {
        inliers = 
            RefineCameraAndPoints(data, num_points_final,
            points_final, projs_final, idxs_final, 
            cameras, added_order, pt_views, &camera_new, 
            true);
    } else {
        inliers = RefineCameraParameters(data, num_points_final, 
            points_final, projs_final, 
            idxs_final, &camera_new, 
            NULL, !m_fixed_focal_length, true,
            m_optimize_for_fisheye,
            m_estimate_distortion,
            m_min_proj_error_threshold,
            m_max_proj_error_threshold);
    }

    if ((int) inliers.size() < 8 || camera_new.f < 0.1 * data.GetWidth()) {
        printf("[BundleInitializeImage] Bad camera\n");
        if (success_out)
            *success_out = false;

        delete [] points_final;
        delete [] projs_final;
        delete [] idxs_final;
        delete [] keys_final;

        delete [] points_solve;
        delete [] projs_solve;
        delete [] projs_solve_orig;
        delete [] idxs_solve;
        delete [] keys_solve;

        camera_params_t dummy;
        return dummy;
    }

    /* Point the keys to their corresponding points */
    num_inliers = (int) inliers.size();
    for (int i = 0; i < num_inliers; i++) {
        int inlier_idx = inliers[i];
        // printf("[BundleInitializeImage] Connecting point [%d]\n",
        //        idxs_final[inlier_idx]);
        data.m_keys[keys_final[inlier_idx]].m_extra = idxs_final[inlier_idx];
        pt_views[idxs_final[inlier_idx]].
            push_back(ImageKey(camera_idx, keys_final[inlier_idx]));
    }
    fflush(stdout);

    delete [] points_final;
    delete [] projs_final;
    delete [] idxs_final;
    delete [] keys_final;

    delete [] points_solve;
    delete [] projs_solve;
    delete [] projs_solve_orig;
    delete [] idxs_solve;
    delete [] keys_solve;

    clock_t end = clock();

    printf("[BundleInitializeImage] Initializing took %0.3fs\n",
        (double) (end - start) / CLOCKS_PER_SEC);

    data.ReadKeyColors();
    data.m_camera.m_adjusted = true;

    return camera_new;
}


void BundlerApp::BundleInitializeImageFullBundle(int image_idx, int parent_idx,
                                                 int num_cameras,
                                                 int num_points, 
                                                 int *added_order,
                                                 camera_params_t *cameras,
                                                 v3_t *points, v3_t *colors,
                                                 std::vector<ImageKeyVector> 
                                                 &pt_views) 
{
    m_image_data[image_idx].LoadKeys(false, !m_optimize_for_fisheye);
    m_image_data[image_idx].ReadKeyColors();

    SetTracks(image_idx);

    /* **** Add any new points **** */
    int num_pts_solve = 0;
    int num_keys = GetNumKeys(image_idx);
    v3_t *points_solve = new v3_t[num_keys];
    v2_t *projs_solve = new v2_t[num_keys];
    int *idxs_solve = new int[num_keys];
    int *keys_solve = new int[num_keys];

    int *saw = new int[num_points];
    for (int i = 0; i < num_points; i++)
        saw[i] = 0;

    int *keys_seen = new int[num_keys];
    for (int i = 0; i < num_keys; i++)
        keys_seen[i] = -1;

    printf("[BundleInitializeImageFullBundle] "
        "Connecting existing matches...\n");

    for (int i = 0; i < num_cameras; i++) {
        int other = added_order[i];

        int first = MIN(image_idx, other);
        int second = MAX(image_idx, other);

        MatchIndex idx = GetMatchIndex(first, second);
        std::vector<KeypointMatch> &list = m_matches.GetMatchList(idx);

        SetMatchesFromTracks(first, second);

        printf("  Matches[%d,%d] = %d\n", image_idx, other,
            (int) m_matches.GetNumMatches(idx));
        // (int) m_match_lists[idx].size());

        for (unsigned int j = 0; j < list.size(); j++) {
            int idx1 = list[j].m_idx1;
            int idx2 = list[j].m_idx2;

            int this_idx, other_idx;
            if (image_idx == first) {
                this_idx = idx1;
                other_idx = idx2;
            } else {
                other_idx = idx1;
                this_idx = idx2;
            }

            if (GetKey(other,other_idx).m_extra < 0) {
                if (GetKey(other,other_idx).m_extra == -2) {
                    // printf("Removed point as outlier\n");
                } else if (GetKey(other,other_idx).m_extra == -1) {
                    // printf("Haven't added point yet\n");
                } else {
                    printf("Error!  Index = %d\n", 
                        GetKey(other,other_idx).m_extra);
                }

                continue;
            } else {
                int pt_idx = GetKey(other,other_idx).m_extra;

                if (saw[pt_idx] == 1) {
                    // printf("Already saw point %d\n", pt_idx);
                    continue;
                }

                saw[pt_idx] = 1;

                if (keys_seen[this_idx] != -1) {
                    printf("Error!  Already saw key %d (pt: %d)!\n", 
                        this_idx, pt_idx);
                }

                printf("  Connecting existing point "
                    "%d => %d [%d] (cam: %d)\n", 
                    this_idx, other_idx, pt_idx, other);

                keys_seen[this_idx] = pt_idx;

                /* Add the point to the set we'll use to solve for
                * the camera position */

                points_solve[num_pts_solve] = points[pt_idx];

                if (m_optimize_for_fisheye) {
                    double x = GetKey(image_idx,this_idx).m_x;
                    double y = GetKey(image_idx,this_idx).m_y;
                    double x_u, y_u;
                    m_image_data[image_idx].UndistortPoint(x, y, x_u, y_u);
                    projs_solve[num_pts_solve] = v2_new(x_u, y_u);
                } else {
                    projs_solve[num_pts_solve] = 
                        v2_new(GetKey(image_idx,this_idx).m_x, 
                        GetKey(image_idx,this_idx).m_y);
                }

                idxs_solve[num_pts_solve] = pt_idx;
                keys_solve[num_pts_solve] = this_idx;

                num_pts_solve++;

                if (num_pts_solve > num_keys) {
                    printf("Error! Number of points exceeds "
                        "the number of keys (%d > %d)!\n", 
                        num_pts_solve, num_keys);

                }
            }
        }

        // ClearMatches(idx);
        m_matches.ClearMatch(idx);
        // m_match_lists[idx].clear();
    }

    fflush(stdout);

    delete [] saw;
    delete [] keys_seen;

    /* **** Solve for the camera position **** */
    double Kinit[9], Rinit[9], tinit[3];
    std::vector<int> inliers, inliers_weak, outliers;
    bool success = 
        FindAndVerifyCamera(num_pts_solve, points_solve, projs_solve,
        idxs_solve, Kinit, Rinit, tinit, 
        m_projection_estimation_threshold, 
        16.0 * m_projection_estimation_threshold, /*4.0*/
        inliers, inliers_weak, outliers);

    printf("[BundleInitializeImageFullBundle] "
        "%d(%d) inliers out of %d\n", 
        (int) inliers.size(), (int) inliers_weak.size(), num_pts_solve);

    if (success) {

#if 0
        /* Try to fit another camera to the outliers */
        int num_outliers = (int) outliers.size();

        if (num_outliers >= 32) {
            v3_t *points_solve2 = new v3_t[num_outliers];
            v2_t *projs_solve2 = new v2_t[num_keys];
            int *idxs_solve2 = new int[num_keys];
            int *keys_solve2 = new int[num_keys];

            for (int i = 0; i < num_outliers; i++) {
                points_solve2[i] = points_solve[outliers[i]];
                projs_solve2[i] = projs_solve[outliers[i]];
                idxs_solve2[i] = idxs_solve[outliers[i]];
                keys_solve2[i] = keys_solve[outliers[i]];
            }

            double Kinit2[9], Rinit2[9], tinit2[3];
            std::vector<int> inliers2, inliers2_weak, outliers2;

            FindAndVerifyCamera(num_outliers, points_solve2, projs_solve2,
                idxs_solve2, Kinit2, Rinit2, tinit2, 
                m_projection_estimation_threshold,
                16.0 * m_projection_estimation_threshold, /*4*/
                inliers2, inliers2_weak, outliers2);

            /* If 50% of the outliers agree on a camera matrix, add
            * these inliers to the full set */

            int num_inliers2 = (int) inliers2.size();
            if (num_inliers2 > 0.5 * num_outliers) {
                printf("[BundleInitializeImageFullBundle] "
                    "Adding %d inliers from second camera\n", 
                    (int) inliers2_weak.size());

                for (int i = 0; i < (int) inliers2_weak.size(); i++) {
                    int idx = inliers2_weak[i];
                    int real_idx = outliers[idx];

                    inliers_weak.push_back(real_idx);
                }
            }

            delete [] points_solve2;
            delete [] projs_solve2;
            delete [] keys_solve2;
            delete [] idxs_solve2;
        }
#endif	

        /* Setup the new camera */
        InitializeCameraParams(m_image_data[image_idx], cameras[num_cameras]);

        double Kinit_inv[9];
        matrix_invert(3, Kinit, Kinit_inv);

        memcpy(cameras[num_cameras].R, Rinit, 9 * sizeof(double));

        matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, 
            cameras[num_cameras].t);
        matrix_scale(3, 1, cameras[num_cameras].t, -1.0, 
            cameras[num_cameras].t);

        /* Connect the inliers to the camera */
        int num_inliers = (int) inliers_weak.size();

        for (int i = 0; i < num_inliers; i++) {
            int key_idx = keys_solve[inliers_weak[i]];
            int pt_idx = idxs_solve[inliers_weak[i]];

            GetKey(image_idx, key_idx).m_extra = pt_idx;
            pt_views[pt_idx].push_back(ImageKey(num_cameras, key_idx));
        }
    } else {
        /* Reset the camera */

        if (parent_idx != -1) {
            memcpy(cameras[num_cameras].R, cameras[parent_idx].R, 
                9 * sizeof(double));
            memcpy(cameras[num_cameras].t, cameras[parent_idx].t, 
                3 * sizeof(double));
        } else {
            printf("[BundleInitializeImageFullBundle] "
                "Error: parent_idx = -1\n");
        }
    }

    fflush(stdout);

    /* Set up camera constraints */
    SetCameraConstraints(added_order[num_cameras], cameras + num_cameras);

    /* Initialize the focal length of the new camera */
    if (m_fixed_focal_length) {
        cameras[num_cameras].f = m_init_focal_length;
    } else {
        if (m_use_focal_estimate) {
            if (m_image_data[image_idx].m_has_init_focal) {
                double ratio;
                double init = m_image_data[image_idx].m_init_focal;
                double obs = 0.5 * (Kinit[0] + Kinit[4]);

                printf("[BundleInitializeImageFullBundle] "
                    "Camera has initial focal length of %0.3f\n", init);

                if (init > obs) ratio = init / obs;
                else            ratio = obs / init;

                if (ratio < 1.4 || m_trust_focal_estimate) {
                    cameras[num_cameras].f = init;
                    if (m_constrain_focal)
                        SetFocalConstraint(m_image_data[added_order[num_cameras]], 
                        cameras + num_cameras);
                } else {
                    printf("[BundleInitializeImageFullBundle] "
                        "Estimated focal length of %0.3f "
                        "is too different\n", obs);

                    cameras[num_cameras].f = 0.5 * (Kinit[0] + Kinit[4]);   
                }
            } else {
                cameras[num_cameras].f = 0.5 * (Kinit[0] + Kinit[4]);
            }
        } else if (m_image_data[image_idx].m_has_init_focal) {
            cameras[num_cameras].f = m_image_data[image_idx].m_init_focal;
        } else {
            if (parent_idx != -1) {
                cameras[num_cameras].f = cameras[parent_idx].f;
            } else {
                printf("[BundleInitializeImageFullBundle] "
                    "Error: parent_idx = -1\n");                    
            }
        }
    }

    cameras[num_cameras].k[0] = cameras[num_cameras].k[1] = 0.0;

    delete [] points_solve;
    delete [] projs_solve;
    delete [] idxs_solve;
    delete [] keys_solve;

    fflush(stdout);

    if (!m_skip_full_bundle) {
        /* Run sfm again to update parameters */
        RunSFM(num_points, num_cameras + 1, 0, false,
            cameras, points, added_order, colors, pt_views);
    }
}

/* Initialize a single image */
void BundlerApp::BundleImage(char *filename, int parent_img)
{
    ImageData data;

    data.m_name = strdup(filename);
    data.m_img = NULL;
    data.m_thumb = NULL;
    data.m_thumb8 = NULL;
    // data.m_wximage = NULL;
    data.m_image_loaded = false;
    data.m_keys_loaded = false;    
    data.m_licensed = true;

    /* Check if we added the camera already */
    int img_idx = (int) m_image_data.size();
    if (data.ReadCamera() && data.ReadTracks(img_idx, m_point_data)) {
        int num_tracks = (int) data.m_visible_points.size();

        for (int i = 0; i < num_tracks; i++) {
            int pt_idx = data.m_visible_points[i];
            double *pos = m_point_data[pt_idx].m_pos;
            double proj[2];

            data.m_camera.Project(pos, proj);

            printf("proj[%d] = %0.3f, %0.3f\n", i, proj[0], proj[1]);
        }

        m_image_data.push_back(data);
        return;
    }

    if (parent_img != -1) {
        double pos[3];
        m_image_data[parent_img].m_camera.GetPosition(pos);
        memcpy(data.m_drop_pt, pos, 3 * sizeof(double));

        BundleRegisterImage(data, true);
    } else {
        // CoalesceFeatureDescriptors();
        // CoalesceFeatureDescriptorsMedian();
        BundleRegisterImage(data, false);	
    }

    data.WriteCamera();
    data.WriteTracks();
}

/* Initialize images read from a file */
void BundlerApp::BundleImagesFromFile(FILE *f)
{
    char buf[256];

    if (!m_add_images_fast)
        CoalesceFeatureDescriptors();
    // CoalesceFeatureDescriptorsMedian();

    while (fgets(buf, 256, f)) {
        ImageData data;

        data.InitFromString(buf, m_image_directory, false);
        data.m_licensed = true;

        /* Read or produce the data */
        int img_idx = (int) m_image_data.size();

        if (data.ReadCamera() && data.ReadTracks(img_idx, m_point_data)) {
            m_image_data.push_back(data);
        } else if (!m_add_images_fast && BundleRegisterImage(data, false)) {
            data.WriteCamera();
            data.WriteTracks();

            int num_keys = data.m_keys.size();
            int img_idx = (int) m_image_data.size();

            /* Add views */
            for (int i = 0; i < num_keys; i++) {
                if (data.m_keys[i].m_extra != -1) {
                    int pt_idx = data.m_keys[i].m_extra;
                    m_point_data[pt_idx].
                        m_views.push_back(ImageKey(img_idx, i));
                }
            }

            data.UnloadKeys();
            m_image_data.push_back(data);
        }
    }
}

std::vector<KeypointMatch> 
RemoveDuplicateMatches(const std::vector<KeypointMatch> &matches)
{
    int num_matches = (int) matches.size();
    std::vector<KeypointMatch> matches_new;

    int num_matches_new = 0;
    for (int i = 0; i < num_matches; i++) {
        int target = matches[i].m_idx2;

        bool duplicate = false;
        for (int j = 0; j < num_matches_new; j++) {
            if (matches_new[j].m_idx2 == target) {
                duplicate = true;
                break;
            }
        }

        if (!duplicate)
            matches_new.push_back(matches[i]);

        num_matches_new = (int) matches_new.size();
    }

    return matches_new;
}

/* Register a new image with the existing model */
bool BundlerApp::BundleRegisterImage(ImageData &data, bool init_location)
{
    printf("[BundleRegisterImage] Registering [%dx%d] image\n",
        data.GetWidth(), data.GetHeight());

    /* Read the keys for this image */
    data.LoadOrExtractKeys(m_sift_binary, !m_optimize_for_fisheye);

    if ((int) data.m_keys_desc.size() == 0) {
        printf("[BundleRegisterImage] "
            "Error: image has no keypoints\n");
        return false;
    }

    clock_t start = clock();

    /* **** Connect the new camera to any existing points **** */
    int num_pts_solve = 0;
    int num_keys = (int) data.m_keys_desc.size();
    v3_t *points_solve = new v3_t[num_keys];
    v2_t *projs_solve = new v2_t[num_keys];
    int *idxs_solve = new int[num_keys];
    int *keys_solve = new int[num_keys];

    for (int i = 0; i < num_keys; i++)
        data.m_keys_desc[i].m_extra = -1;

    if (init_location) {
        /* Find the closest registered cameras to this image */
#ifdef __USE_ANN__
        ANNkd_tree *tree = CreateCameraSearchTree();
#else
        BruteForceSearch *search = CreateCameraSearchTree();
#endif    

#define NUM_NNS 20
        int nn_idxs[NUM_NNS];

        v3_t q = 
            v3_new(data.m_drop_pt[0], data.m_drop_pt[1], data.m_drop_pt[2]);
#ifdef __USE_ANN__
        float dists[NUM_NNS];
        // query_ann_tree_3D_idx(tree, q, 0.0, NUM_NNS, nn_idxs, dists);
        float query[3] = { Vx(q), Vy(q), Vz(q) };
        tree->annkPriSearch(query, NUM_NNS, nn_idxs, dists, 0.0);
#else
        double dists[NUM_NNS];
        search->GetClosestPoints(q, NUM_NNS, nn_idxs, dists);
#endif

        /* Initialize point indices to -1 */
        printf("[BundleRegisterImage] Initializing indices...\n");
        for (int i = 0; i < NUM_NNS; i++) {
            int cam_idx = GetRegisteredCameraIndex(nn_idxs[i]);

            printf("NN[%d] = %d\n", i, cam_idx);

            m_image_data[cam_idx].LoadKeys(true, !m_optimize_for_fisheye);

            int num_keys2 = (int) m_image_data[cam_idx].m_keys.size();
            for (int j = 0; j < num_keys2; j++) {
                m_image_data[cam_idx].m_keys_desc[j].m_extra = -1;
            }
        }

#if 1
        /* Set the point indices */
        printf("[BundleRegisterImage] Set point indices...\n");

        for (int i = 0; i < NUM_NNS; i++) {
            int cam_idx = GetRegisteredCameraIndex(nn_idxs[i]);
            int num_visible_points = 
                m_image_data[cam_idx].m_visible_points.size();

            for (int j = 0; j < num_visible_points; j++) {
                int pt_idx = m_image_data[cam_idx].m_visible_points[j];

                int num_views = m_point_data[pt_idx].m_views.size();
                for (int k = 0; k < num_views; k++) {
                    int img = m_point_data[pt_idx].m_views[k].first;

                    if (img == cam_idx) {
                        int key_idx = m_point_data[pt_idx].m_views[k].second;
                        GetKey(cam_idx, key_idx).m_extra = pt_idx;
                    }
                }
            }
        }
#else
        for (int i = 0; i < (int) m_point_data.size(); i++) {
            for (int j = 0; j < (int) m_point_data[i].m_views.size(); j++) {
                ImageKey p = m_point_data[i].m_views[j];

                for (int k = 0; k < NUM_NNS; k++) {
                    int cam_idx = GetRegisteredCameraIndex(nn_idxs[k]);

                    if (cam_idx == p.first)
                        GetKey(p.first,p.second).m_extra = i;
                }
            }
        }
#endif

        /* Match the new keys to each of the nearby images */
        printf("[BundleRegisterImage] Matching images...\n");

        typedef std::vector<KeypointMatch> KeypointMatchList;
        std::vector<KeypointMatchList> match_lists;

        for (int i = 0; i < NUM_NNS; i++) {
            int cam_idx = GetRegisteredCameraIndex(nn_idxs[i]);

            printf("[BundleRegisterImage] "
                "Comparing to image %d...\n",
                cam_idx);

            // m_image_data[cam_idx].LoadKeys();

            KeypointMatchList matches;
            matches = 
                MatchKeys(data.m_keys_desc, 
                m_image_data[cam_idx].m_keys_desc, 
                true, 0.75);

            KeypointMatchList matches_sym;
            matches_sym = 
                MatchKeys(m_image_data[cam_idx].m_keys_desc, 
                data.m_keys_desc,
                false, 1.0);

            if (matches_sym.size() != m_image_data[cam_idx].m_keys_desc.size())
                printf("Error: not enough matches\n");

            /* Prune asymmetric matches */
            KeypointMatchList matches_new;
            for (int j = 0; j < (int) matches.size(); j++) {
                int idx1 = matches[j].m_idx1;
                int idx2 = matches[j].m_idx2;

                if (matches_sym[idx2].m_idx2 != idx1)
                    continue;

                matches_new.push_back(matches[j]);
            }

            matches = matches_new;

            /* Remove duplicates */
            printf("Found %d matches [before pruning]\n", 
                (int) matches.size());

            matches = RemoveDuplicateMatches(matches);
            printf("Found %d matches [after pruning]\n", 
                (int) matches.size());

            /* Estimate a fundamental matrix */
            double F[9];
            std::vector<int> inliers = 
                EstimateFMatrix(data.m_keys_desc, 
                m_image_data[cam_idx].m_keys_desc, 
                matches, 
                m_fmatrix_rounds, m_fmatrix_threshold, F);

            int num_inliers = (int) inliers.size();
            printf("Inliers[%d] = %d out of %d\n", 
                cam_idx, num_inliers, (int) matches.size());

            /* Refine the matches */
            KeypointMatchList new_matches;

            for (int i = 0; i < num_inliers; i++) {
                new_matches.push_back(matches[inliers[i]]);
            }

            match_lists.push_back(new_matches);

            fflush(stdout);
        }

        /* Now, set up the data structures for bundle adjustment */
        int curr_num_pts = (int) m_point_data.size();

        int *saw = new int[curr_num_pts];
        for (int i = 0; i < curr_num_pts; i++) {
            saw[i] = 0;
        }

        printf("[BundleRegisterImage] "
            "Connecting existing matches...\n");

        for (int i = 0; i < NUM_NNS; i++) {
            int other = GetRegisteredCameraIndex(nn_idxs[i]);

            for (int j = 0; j < (int) match_lists[i].size(); j++) {
                int idx1 = match_lists[i][j].m_idx1;
                int idx2 = match_lists[i][j].m_idx2;

                if (GetKey(other,idx2).m_extra < 0) {
                    continue;
                } else {
                    int pt_idx = GetKey(other,idx2).m_extra;

                    if (saw[pt_idx] == 1)
                        continue;

                    saw[pt_idx] = 1;

                    printf("  Connecting existing point "
                        "%d => %d [%d] (cam: %d)\n", 
                        idx1, idx2, pt_idx, other);

                    /* This is an old point */
                    data.m_keys_desc[idx1].m_extra = pt_idx;

                    /* Add the point to the set we'll use to solve for
                    * the camera position */
                    double *pt = m_point_data[pt_idx].m_pos;

                    points_solve[num_pts_solve] = 
                        v3_new(pt[0], pt[1], pt[2]);

                    if (m_optimize_for_fisheye) {
                        double x = data.m_keys_desc[idx1].m_x;
                        double y = data.m_keys_desc[idx1].m_y;
                        double x_u, y_u;
                        m_image_data[idx1].UndistortPoint(x, y, x_u, y_u);
                        projs_solve[num_pts_solve] = v2_new(x_u, y_u);
                    } else {
                        projs_solve[num_pts_solve] = 
                            v2_new(data.m_keys_desc[idx1].m_x, 
                            data.m_keys_desc[idx1].m_y);
                    }

                    idxs_solve[num_pts_solve] = pt_idx;
                    keys_solve[num_pts_solve] = idx1;

                    num_pts_solve++;
                }
            }
        }

        fflush(stdout);

        delete [] saw;

#ifdef __USE_ANN__
        delete tree;
#else
        delete search;
#endif
    } else { /* !init_location */
        /* Match keys to all points */
#if 0
        std::vector<KeypointMatch> matches = 
            MatchKeysToPoints(data.m_keys, 0.6);
#else
        std::vector<KeypointMatch> matches = 
            MatchPointsToKeys(data.m_keys_desc, 0.7);
#endif

        matches = RemoveDuplicateMatches(matches);

        num_pts_solve = (int) matches.size();
        printf("[BundleRegisterImage] Found %d matches\n",
            num_pts_solve);

        delete [] points_solve;
        delete [] projs_solve;
        delete [] idxs_solve;
        delete [] keys_solve;

        points_solve = new v3_t[num_pts_solve];
        projs_solve = new v2_t[num_pts_solve];
        idxs_solve = new int[num_pts_solve];
        keys_solve = new int[num_pts_solve];

        for (int i = 0; i < num_pts_solve; i++) {
#if 0
            int key_idx = matches[i].m_idx1;
            int pt_idx = matches[i].m_idx2;
#else
            int key_idx = matches[i].m_idx2;
            int pt_idx = matches[i].m_idx1;
#endif

            const PointData &pt = m_point_data[pt_idx];

            points_solve[i] = v3_new(pt.m_pos[0], pt.m_pos[1], pt.m_pos[2]);

            if (m_optimize_for_fisheye) {
                double x = data.m_keys_desc[key_idx].m_x;
                double y = data.m_keys_desc[key_idx].m_y;
                double x_u, y_u;
                data.UndistortPoint(x, y, x_u, y_u);
                projs_solve[i] = v2_new(x_u, y_u);
            } else {
                projs_solve[i] = 
                    v2_new(data.m_keys_desc[key_idx].m_x, 
                    data.m_keys_desc[key_idx].m_y);
            }

            keys_solve[i] = key_idx;
            idxs_solve[i] = pt_idx;
        }
    }

    /* **** Solve for the camera position **** */
    printf("[BundleRegisterImage] Initializing camera...\n");

    camera_params_t *camera_new = new camera_params_t;
    InitializeCameraParams(data, *camera_new);
    ClearCameraConstraints(camera_new);

    double Kinit[9], Rinit[9], tinit[3];
    std::vector<int> inliers, inliers_weak, outliers;
    bool success = 
        FindAndVerifyCamera(num_pts_solve, points_solve, projs_solve,
        idxs_solve, Kinit, Rinit, tinit, 
        m_projection_estimation_threshold, 
        2.0 * m_projection_estimation_threshold,
        inliers, inliers_weak, outliers);


    ClearCameraConstraints(camera_new);

    if (success) {
        /* Set up the new camera */
        memcpy(camera_new->R, Rinit, 9 * sizeof(double));

        matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, camera_new->t);
        matrix_scale(3, 1, camera_new->t, -1.0, camera_new->t);

        // camera_new->f = 0.5 * (Kinit[0] + Kinit[4]);
        if (m_fixed_focal_length) {
            camera_new->f = m_init_focal_length;
        } else {
            if (m_use_focal_estimate) {
                if (data.m_has_init_focal) {
                    double ratio;
                    double init = data.m_init_focal;
                    double obs = 0.5 * (Kinit[0] + Kinit[4]);

                    printf("[BundleRegisterImage] "
                        "Camera has initial focal length of %0.3f\n", init);

                    if (init > obs) ratio = init / obs;
                    else            ratio = obs / init;

                    if (ratio < 1.4 || m_trust_focal_estimate) {
                        camera_new->f = init;
                        if (m_constrain_focal)
                            SetFocalConstraint(data, camera_new);
                    } else {
                        printf("[BundleRegisterImage] "
                            "Estimated focal length of %0.3f "
                            "is too different\n", obs);

                        camera_new->f = 0.5 * (Kinit[0] + Kinit[4]);   
                    }
                } else {
                    camera_new->f = 0.5 * (Kinit[0] + Kinit[4]);
                }
            } else if (data.m_has_init_focal) {
                camera_new->f = data.m_init_focal;
            } else {
                camera_new->f = m_init_focal_length; // cameras[parent_idx].f;
            }
        }
    } else {
        delete camera_new;

        delete [] points_solve;
        delete [] projs_solve;
        delete [] idxs_solve;
        delete [] keys_solve;

        return false;
    }

    /* **** Finally, start the bundle adjustment **** */
    printf("[BundleRegisterImage] Adjusting [%d,%d]...\n",
        (int) inliers.size(), (int) inliers_weak.size());

    int num_points_final = (int) inliers_weak.size();
    v3_t *points_final = new v3_t[num_points_final];
    v2_t *projs_final = new v2_t[num_points_final];

    for (int i = 0; i < num_points_final; i++) {
        points_final[i] = points_solve[inliers_weak[i]];
        projs_final[i] = projs_solve[inliers_weak[i]];
    }

    std::vector<int> inliers_final;
#if 1
    inliers_final = RefineCameraParameters(data, num_points_final, 
        points_final, projs_final, 
        NULL, camera_new, 
        NULL, !m_fixed_focal_length, true,
        m_optimize_for_fisheye,
        m_estimate_distortion,
        m_min_proj_error_threshold,
        m_max_proj_error_threshold);
#else
    inliers = 
        RefineCameraAndPoints(data, num_points_final,
        points_final, projs_final, idxs_final, 
        cameras, added_order, pt_views, &camera_new, 
        true);
#endif

    int num_inliers = (int) inliers_final.size();

    success = false;
#define MIN_INLIERS_ADD_IMAGE 16
    if (num_inliers < 0.5 * num_points_final) {
        printf("[BundleRegisterImage] "
            "Threw out too many outliers [%d / %d]\n", 
            num_inliers, num_points_final);
    } else if (num_inliers >= MIN_INLIERS_ADD_IMAGE) {
        printf("[BundleRegisterImage] Final camera parameters:\n");
        printf("f: %0.3f\n", camera_new->f);
        printf("R:\n");
        matrix_print(3, 3, camera_new->R);
        printf("t: ");
        matrix_print(1, 3, camera_new->t);

        fflush(stdout);

#if 0
        data.LoadImage();

        for (int i = 0; i < num_points_final; i++) {
            img_draw_pt(data.m_img, 
                iround(Vx(projs_final[i]) + 0.5 * data.GetWidth()),
                iround(Vy(projs_final[i]) + 0.5 * data.GetHeight()), 
                6, 0xff, 0x0, 0x0);
        }

        img_write_bmp_file(data.m_img, "projs.bmp");
#endif

        data.m_camera.m_adjusted = true;

        memcpy(data.m_camera.m_R, camera_new->R, 9 * sizeof(double));

        data.m_drop_pt[0] = camera_new->t[0];
        data.m_drop_pt[1] = camera_new->t[1];
        data.m_drop_pt[2] = camera_new->t[2];

        data.m_camera.SetPosition(data.m_drop_pt);
        data.m_camera.m_focal = camera_new->f;
        data.m_camera.Finalize();

        /* Connect the keys to the points */
        int num_inliers_final = (int) inliers_final.size();

        for (int i = 0; i < num_keys; i++) {
            data.m_keys_desc[i].m_extra = -1;
        }

        for (int i = 0; i < num_inliers_final; i++) {
            int key_idx = keys_solve[inliers_weak[inliers_final[i]]];
            int pt_idx = idxs_solve[inliers_weak[inliers_final[i]]];

            printf("(k,p)[%d] = %d, %d\n", i, key_idx, pt_idx);
            double proj[2];
            data.m_camera.Project(m_point_data[pt_idx].m_pos, proj);
            printf("pr[%d] = %0.3f, %0.3f == %0.3f, %0.3f\n", 
                i, proj[0], proj[1], 
                data.m_keys_desc[key_idx].m_x, 
                data.m_keys_desc[key_idx].m_y);

            data.m_keys_desc[key_idx].m_extra = pt_idx;
            data.m_visible_points.push_back(pt_idx);
            data.m_visible_keys.push_back(key_idx);
        }

        success = true;
    }

    delete camera_new;

    delete [] points_final;
    delete [] projs_final;

    delete [] points_solve;
    delete [] projs_solve;
    delete [] idxs_solve;
    delete [] keys_solve;

    clock_t end = clock();

    printf("[BundleRegisterImage] Registration took %0.3fs\n",
        (double) (end - start) / (double) CLOCKS_PER_SEC);

    return success;
}

int BundlerApp::RemoveBadPointsAndCameras(int num_points, int num_cameras, 
                                          int *added_order, 
                                          camera_params_t *cameras,
                                          v3_t *points, v3_t *colors,
                                          std::vector<ImageKeyVector> &pt_views)
{
    int num_pruned = 0;

    for (int i = 0; i < num_points; i++) {
        double *pos = points[i].p;
        int num_views = (int) pt_views[i].size();

        if (num_views == 0)
            continue;

        double max_angle = 0.0;
        for (int j = 0; j < num_views; j++) {
            int v1 = pt_views[i][j].first;

            double r1[3];
            matrix_diff(3, 1, 3, 1, pos, cameras[v1].t, r1);
            double norm = matrix_norm(3, 1, r1);
            matrix_scale(3, 1, r1, 1.0 / norm, r1);

            for (int k = j+1; k < num_views; k++) {
                int v2 = pt_views[i][k].first;

                double r2[3];
                matrix_diff(3, 1, 3, 1, pos, cameras[v2].t, r2);
                double norm = matrix_norm(3, 1, r2);
                matrix_scale(3, 1, r2, 1.0 / norm, r2);

                double dot;
                matrix_product(1, 3, 3, 1, r1, r2, &dot);

                double angle = 
                    acos(CLAMP(dot, -1.0 + 1.0e-8, 1.0 - 1.0e-8));

                if (angle > max_angle) {
                    max_angle = angle;
                }
            }
        }

        if (RAD2DEG(max_angle) < 0.5 * m_ray_angle_threshold) {
            printf("[RemoveBadPointsAndCamera] "
                "Removing point %d with angle %0.3f\n",
                i, RAD2DEG(max_angle));

            for (int j = 0; j < num_views; j++) {
                // Set extra flag back to 0
                int v = pt_views[i][j].first;
                int k = pt_views[i][j].second;
                GetKey(added_order[v], k).m_extra = -1;
            }

            pt_views[i].clear();

            if (colors != NULL) {
                Vx(colors[i]) = 0x0;
                Vy(colors[i]) = 0x0;
                Vz(colors[i]) = 0xff;
            }

            num_pruned++;
        }
    }

    printf("[RemoveBadPointsAndCameras] Pruned %d points\n", num_pruned);

    return num_pruned;
}
