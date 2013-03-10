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

/* BundleFast.cpp */
/* Routines for faster structure from motion */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "defines.h"
#include "matrix.h"
#include "triangulate.h"
#include "util.h"

#include "BundlerApp.h"
#include "Bundle.h"

#define MIN_INLIERS_EST_PROJECTION 15 /* 30 */ /* This constant needs
						* adjustment */
#define INIT_REPROJECTION_ERROR 16.0 /* 6.0 */ /* 8.0 */

/* Quickly compute pose of all cameras */
void BundlerApp::BundleAdjustFast() 
{
    clock_t start = clock();

    /* Compute initial image information */
    ComputeGeometricConstraints();

#if 0
    /* Read keypoints */
    printf("[SifterApp::BundleAdjust] Reading key colors...\n");
    ReadKeyColors();
#endif

    /* Set track pointers to -1 */
    for (int i = 0; i < (int) m_track_data.size(); i++) {
	m_track_data[i].m_extra = -1;
    }

    /* For now, assume all images form one connected component */
    int num_images = GetNumImages();
    int *added_order = new int[num_images];
    int *added_order_inv = new int[num_images];


    /* **** Run bundle adjustment! **** */

    camera_params_t *cameras = new camera_params_t[num_images];
    int max_pts = (int) m_track_data.size(); // 1243742; /* HACK! */
    v3_t *points = new v3_t[max_pts];
    v3_t *colors = new v3_t[max_pts];
    std::vector<ImageKeyVector> pt_views;

#if 0
    if (use_constraints) {
	/* Constrain only the z-coordinate */
	char constrained[7] = { 0, 1, 1, 0, 0, 0, 0 };
	double constraints[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for (int i = 0; i < num_images; i++) {
	    memcpy(cameras[i].constrained, constrained, 7);
	    memcpy(cameras[i].constraints, constraints, 7 * sizeof(double));
	}
    }
#endif

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

	printf("[SifterApp::BundleAdjust] Adjusting cameras "
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

	/* Run sfm for the first time */
	double error0;
	error0 = RunSFM(curr_num_pts, 2, 0, false,
			cameras, points, added_order, colors, pt_views);

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
				    -1.0 };
		
		double proj2[3] = { GetKey(added_order[1],k2).m_x,
				    GetKey(added_order[1],k2).m_y,
				    -1.0 };

		double proj1_norm[3], proj2_norm[3];
		
		matrix_product(3, 3, 3, 1, K1inv, proj1, proj1_norm);
		matrix_product(3, 3, 3, 1, K2inv, proj2, proj2_norm);

		v2_t p = v2_new(proj1_norm[0] / proj1_norm[2],
				proj1_norm[1] / proj1_norm[2]);
	    
		v2_t q = v2_new(proj2_norm[0] / proj2_norm[2],
				proj2_norm[1] / proj2_norm[2]);

		double proj_error;
		
                double t1[3];
                double t2[3];
			
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
	    printf("[SifterApp::BundleAdjust] Error: initial good pair "
		   "not provided!\n");
	    printf("[SifterApp::BundleAdjust] Please specify a pair of "
		   "cameras with medium baseline using\n"
		   "  --init_pair1 <img1> and --init_pair2 <img2>\n");
	    exit(1);
	}
    
	good_pair_1 = added_order_inv[m_initial_pair[0]];
	good_pair_2 = added_order_inv[m_initial_pair[1]];

	if (good_pair_1 == -1 || good_pair_2 == -1) {
	    printf("[SifterApp::BundleAdjust] Error: initial pair haven't "
		   "been adjusted!\n");
	    printf("[SifterApp::BundleAdjust] Please specify another pair!\n");
	    exit(0);
	}
#endif

	curr_num_cameras = num_init_cams;
	pt_count = curr_num_pts = (int) m_point_data.size();
    }
    
    int round = 0;
    while (curr_num_cameras < num_images) {
	int parent_idx;
	int max_cam = 
            FindCameraWithMostMatches(curr_num_cameras, curr_num_pts, 
                                      added_order, parent_idx, 
                                      max_matches, pt_views);

	printf("[SifterApp::BundleAdjust] max_matches = %d\n", max_matches);

	if (max_matches < m_min_max_matches)
	    break; /* No more connections */

	/* Find all images with 90% of the matches of the maximum */
	std::vector<ImagePair> image_set;

        if (false && max_matches < 48) {
            image_set.push_back(ImagePair(max_cam, parent_idx));
        } else {
            // int nMatches = MIN(100, iround(0.75 /* 0.9 */ * max_matches));
            int nMatches = iround(0.75 /* 0.9 */ * max_matches);
	    image_set = 
                FindCamerasWithNMatches(nMatches,
                                        curr_num_cameras, curr_num_pts, 
                                        added_order, pt_views);
        }
        
	int num_added_images = (int) image_set.size();

	printf("[SifterApp::BundleAdjustFast] Registering %d images\n",
	       num_added_images);

	for (int i = 0; i < num_added_images; i++)
	    printf("[SifterApp::BundleAdjustFast] Adjusting camera %d\n",
		   image_set[i].first);

	/* Now, throw the new cameras into the mix */
        int image_count = 0;
	for (int i = 0; i < num_added_images; i++) {
	    int next_idx = image_set[i].first;
	    int parent_idx = image_set[i].second;

	    added_order[curr_num_cameras + image_count] = next_idx;

	    printf("[SifterApp::BundleAdjust[%d]] Adjusting camera %d "
		   "(parent = %d)\n", 
		   round, next_idx, 
                   (parent_idx == -1 ? -1 : added_order[parent_idx]));

	    /* **** Set up the new camera **** */
            bool success = false;
            camera_params_t camera_new = 
		BundleInitializeImage(m_image_data[next_idx], 
				      next_idx, curr_num_cameras + image_count,
				      curr_num_cameras, curr_num_pts,
				      added_order, points, 
				      NULL /*cameras + parent_idx*/, cameras, 
				      pt_views, &success);

            if (success) {
                cameras[curr_num_cameras+image_count] = camera_new;
                image_count++;
            } else {
                printf("[BundleAdjust] Couldn't initialize image %d\n",
                       next_idx);
                m_image_data[next_idx].m_ignore_in_bundle = true;
            }
	}
	

	/* Compute the distance between the first pair of cameras */
#if 0
	double dist0 = GetCameraDistance(cameras + good_pair_1, 
					 cameras + good_pair_2, 
					 m_explicit_camera_centers);
#else
        double dist0 = 0.0;
#endif

	printf("[SifterApp::BundleAdjust] Adding new matches\n");

	pt_count = curr_num_pts;
#if 0
	for (int i = 0; i < num_added_images; i++) {
	    pt_count = 
		BundleAdjustAddNewPoints(curr_num_cameras + i, 
					 pt_count, curr_num_cameras + i,
					 added_order, cameras, points, colors,
					 dist0, pt_views);
	}
#endif
	
	curr_num_cameras += image_count;

        if (!m_skip_add_points) {
            pt_count = 
                BundleAdjustAddAllNewPoints(pt_count, curr_num_cameras,
                                            added_order, cameras, 
                                            points, colors,
                                            dist0, pt_views);
        }
        
	curr_num_pts = pt_count;

	printf("[SifterApp::BundleAdjust] Number of points = %d\n", pt_count);
	fflush(stdout);

        if (!m_skip_full_bundle) {
            /* Run sfm again to update parameters */
            RunSFM(curr_num_pts, curr_num_cameras, 0, false,
                   cameras, points, added_order, colors, pt_views);

            /* Remove bad points and cameras */
            RemoveBadPointsAndCameras(curr_num_pts, curr_num_cameras + 1, 
                                      added_order, cameras, points, colors, 
                                      pt_views);

            printf("  focal lengths:\n");
	
            for (int i = 0; i < curr_num_cameras; i++) {
                if(m_image_data[added_order[i]].m_has_init_focal) {
                    printf("   [%03d] %0.3f (%0.3f) %s %d; %0.3e %0.3e\n", 
                           i, cameras[i].f, 
                           m_image_data[added_order[i]].m_init_focal,
                           m_image_data[added_order[i]].m_name,
                           added_order[i], cameras[i].k[0], cameras[i].k[1]);
                } else {
                    printf("   [%03d] %0.3f %s %d; %0.3e %0.3e\n", 
                           i, cameras[i].f, 
                           m_image_data[added_order[i]].m_name,
                           added_order[i], cameras[i].k[0], cameras[i].k[1]);
                }
                // printf("   [%03d] %0.3f\n", i, cameras[i].f);
            }
        }

#if 0
	printf("  extrinsics:\n");
	for (int i = 0; i < curr_num_cameras; i++) {
	    printf("   [%03d] %0.3f %0.3f %0.3f\n"
		   "         %0.3f %0.3f %0.3f\n"
		   "         %0.3f %0.3f %0.3f\n", i,
		   cameras[i].R[0], cameras[i].R[1], cameras[i].R[2],
		   cameras[i].R[3], cameras[i].R[4], cameras[i].R[5],
		   cameras[i].R[6], cameras[i].R[7], cameras[i].R[8]);
	    printf("         %0.3f %0.3f %0.3f\n", 
		   cameras[i].t[0], cameras[i].t[1], cameras[i].t[2]);
	}
#endif
	
	/* Dump output for this round */
	char buf[256];
	sprintf(buf, "points%03d.ply", curr_num_cameras);

	DumpPointsToPly(m_output_directory, buf, 
                        curr_num_pts, curr_num_cameras, 
                        points, colors, cameras);

	if (m_bundle_output_base != NULL) {
	    sprintf(buf, "%s%03d.out", m_bundle_output_base, 
                    curr_num_cameras);
	    DumpOutputFile(m_output_directory, buf, 
                           num_images, curr_num_cameras, curr_num_pts,
			   added_order, cameras, points, colors, pt_views);
	}

	round++;
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

    bool *image_mask = new bool[num_images];

    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted)
	    image_mask[i] = true;
	else 
	    image_mask[i] = false;
    }
}
