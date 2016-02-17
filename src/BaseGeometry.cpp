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

/* BaseGeometry.cpp */
/* Various geometric functions */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "BaseApp.h"
#include "BundleUtil.h"

#include "defines.h"
#include "fit.h"
#include "matrix.h"
#include "qsort.h"

bool BaseApp::ReadTrackPairs(const char *filename)
{
    FILE *f = fopen(filename, "rb");
    if (!f) {
        printf("[ReadTrackPairs] Error opening file %s for reading\n", 
               filename);
        return false;
    }
    
    int num_images;
    fread(&num_images, sizeof(int), 1, f);
    if (num_images != GetNumImages()) {
        printf("[ReadTrackPairs] Mismatch in number of images!\n");
        return false;
    }

    unsigned long num_track_pairs = 0;
    for (int i = 0; i < num_images; i++) {
        int num_nbrs = 0;
        fread(&num_nbrs, sizeof(int), 1, f);
        num_track_pairs += num_nbrs;

        for (int j = 0; j < num_nbrs; j++) {
            int nbr = 0;
            fread(&nbr, sizeof(int), 1, f);

            SetMatch(i, nbr);
            SetMatch(nbr, i);
        }
    }

    fclose(f);

    printf("[ReadTrackPairs] Read %lu track pairs\n", num_track_pairs);
    fflush(stdout);

    return true;
}

void BaseApp::WriteTrackPairs(const char *filename)
{
    FILE *f = fopen(filename, "wb");
    if (!f) {
        printf("[WriteTrackPairs] Error opening file %s for writing\n", 
               filename);
        return;
    }

    int num_images = GetNumImages();
    
    fwrite(&num_images, sizeof(int), 1, f);
    for (int i = 0; i < num_images; i++) {
        // std::list<unsigned int> &nbrs = m_matches.GetNeighbors(i);
        MatchAdjList &nbrs = m_matches.GetNeighbors(i);
        int num_nbrs = (int) nbrs.size();
        
        fwrite(&num_nbrs, sizeof(int), 1, f);
        // std::list<unsigned int>::iterator iter;
        MatchAdjList::iterator iter;
        for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
            // int nbr = *iter;
            int nbr = iter->m_index;
            fwrite(&nbr, sizeof(int), 1, f);
        }
    }

    fclose(f);
}

void BaseApp::ReadGeometricConstraints(const char *filename) {
    FILE *f = fopen(filename, "r");

    int num_images; // = GetNumImages();

    fscanf(f, "%d\n", &num_images);
    if (num_images != GetNumImages()) {
	printf("[ReadGeometricConstraints] Error: number of "
	       "images don't match!\n");
	return;
    }

#if 0
    if (m_transforms == NULL) {
	m_transforms = new TransformInfo * [num_images * num_images];

	for (int i = 0; i < num_images * num_images; i++) {
	    m_transforms[i] = NULL;
	}
    }
#endif

    m_transforms.clear();

    // m_match_lists.clear();
    // m_matches.Clear();
    RemoveAllMatches();

    // for (int i = 0; i < num_images; i++) {
    //    for (int j = 0; j < num_images; j++) {   
    //         int match;
    //         fscanf(f, "%d", &match);

    /* Read the number of transforms */
    unsigned long long num_transforms;
    fscanf(f, "%llu\n", &num_transforms);

    for (unsigned long int count = 0; count < num_transforms; count++) {
        int i, j;
        fscanf(f, "%d %d\n", &i, &j);

        MatchIndex idx = GetMatchIndex(i, j);

        SetMatch(i, j);
        // m_matches[idx] = true;
        m_transforms[idx] = TransformInfo();

        /* Read the transform information */
        m_transforms[idx].ReadFromFile(f);

        /* Read the matches */
        // m_match_lists[idx].clear();

        int num_matches;
        fscanf(f, "%d\n", &num_matches);

        for (int k = 0; k < num_matches; k++) {
            int idx1, idx2;

            fscanf(f, "%d %d\n", &idx1, &idx2);
            
#if 0
            KeypointMatch match;
            
            match.m_idx1 = idx1;
            match.m_idx2 = idx2;
		    
            m_match_lists[idx].push_back(match);
#endif
        }
    }

    /* Read the tracks */
    int num_tracks = 0;
    fscanf(f, "%d", &num_tracks);

    printf("[ReadGeometricConstraints] Reading %d tracks\n", 
	   num_tracks);

    m_track_data.clear();

    clock_t start = clock();
    int count = 0;
    for (int i = 0; i < num_tracks; i++) {
	TrackData track;
	track.Read(f);

        int num_views = (int) track.m_views.size();

        if (num_views < m_min_track_views)
            continue;

        if (num_views > m_max_track_views)
            continue;

#if 0
	for (int j = 0; j < num_views; j++) {
	    int img = track.m_views[j].first;
	    int key = track.m_views[j].second;

	    GetKey(img, key).m_track = i;
	}
#endif
	for (int j = 0; j < num_views; j++) {
	    int img = track.m_views[j].first;
            int key = track.m_views[j].second;

            // assert(key >= 0 && key < m_image_data[img].GetNumKeys());

            m_image_data[img].m_visible_points.push_back(count);
            m_image_data[img].m_visible_keys.push_back(key);
	}

	m_track_data.push_back(track);
        count++;
    }
    clock_t end = clock();
    printf("[ReadGeometricConstraints] Reading tracks took %0.3fs\n",
           (double) (end - start) / CLOCKS_PER_SEC);

#if 1 
    if (!ReadTrackPairs("track-pairs.txt")) {
        start = clock();
        num_tracks = (int) m_track_data.size();
        for (int i = 0; i < num_tracks; i++) {
            if ((i % 10000) == 0) {
                printf("[ReadGeometricConstraints] "
                       "Processing track %d...\n", i);
                fflush(stdout);
            }

            /* Set match flags */
            // for (int j = 0; j < num_views; j++) {
            TrackData &track = m_track_data[i];
            
            ImageKeyVector::iterator iter1, iter2;
            for (iter1 = track.m_views.begin(); 
                 iter1 != track.m_views.end(); iter1++) {
                int img1 = iter1->first; // track.m_views[j].first;
                
                // assert(img1 >= 0 && img1 < num_images);
                
                // for (int k = j+1; k < num_views; k++) {
                for (iter2 = iter1+1; iter2 != track.m_views.end(); iter2++) {
                    int img2 = iter2->first; // track.m_views[k].first;
                    
                    // assert(img2 >= 0 && img2 < num_images);
                
                    if (img1 < img2) {
                        SetMatch(img1, img2);
                        SetMatch(img2, img1);
                    }
                }
            }
        }

        WriteTrackPairs("track-pairs.txt");
        end = clock();
        
        printf("[ReadGeometricConstraints] Computing track pairs "
               "took %0.3fs\n", (double) (end - start) / CLOCKS_PER_SEC);
    }
#endif

    fclose(f);

    // SetMatchesFromTracks();

    // WriteGeometricConstraints("constraints_test.txt");
}

void BaseApp::WriteGeometricConstraints(const char *filename) {
    FILE *f = fopen(filename, "w");

    if (f == NULL) {
	printf("Error opening file %s for writing\n", filename);
	return;
    }

    unsigned int num_images = GetNumImages();

    fprintf(f, "%d\n", num_images);

    /* Count the number of transforms to write */
    unsigned long long num_transforms = 0;
    for (unsigned int i = 0; i < num_images; i++) {
	// for (int j = 0; j < num_images; j++) {
        // std::list<unsigned int> &nbrs = m_matches.GetNeighbors(i);
        MatchAdjList &nbrs = m_matches.GetNeighbors(i);
        int num_nbrs = (int) nbrs.size();
        printf("num_nbrs[%d] = %d\n", i, num_nbrs);

        // std::list<unsigned int>::iterator iter;
        MatchAdjList::iterator iter;
        for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
            unsigned int j = iter->m_index; // *iter; // nbrs[nbr];

            MatchIndex idx = GetMatchIndex(i, j);

            if (m_transforms.find(idx) != m_transforms.end()) {
                num_transforms++;
            }
        }
    }

    printf("[WriteGeometricConstraints] Writing %llu transforms\n", 
           num_transforms);

    fprintf(f, "%llu\n", num_transforms);

    for (unsigned int i = 0; i < num_images; i++) {
        // std::list<unsigned int> &nbrs = m_matches.GetNeighbors(i);
        MatchAdjList &nbrs = m_matches.GetNeighbors(i);

        // std::list<unsigned int>::iterator iter;
        MatchAdjList::iterator iter;
        for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
            unsigned int j = iter->m_index; // *iter; // nbrs[nbr];

            MatchIndex idx = GetMatchIndex(i, j);

            if (m_transforms.find(idx) != m_transforms.end()) {
		fprintf(f, "%d %d\n", i, j);

#if 0
		/* Write the transform information */
		if (m_transforms.find(idx) == m_transforms.end()) {
		    m_transforms[idx] = TransformInfo();
		}
#endif

		m_transforms[idx].WriteToFile(f);

#if 0
		/* Write the matches */
		int num_matches = (int) m_match_lists[idx].size();
		fprintf(f, "%d\n", num_matches);

		for (int k = 0; k < num_matches; k++) {
		    fprintf(f, "%d %d\n", 
			    m_match_lists[idx][k].m_idx1,
			    m_match_lists[idx][k].m_idx2);
		}
#else
                fprintf(f, "0\n");
#endif
                // } else {
		//      fprintf(f, "0\n");
	    }
	}
    }

    /* Write the tracks */
    int num_tracks = (int) m_track_data.size();
    fprintf(f, "%d\n", num_tracks);
    for (int i = 0; i < num_tracks; i++) {
	m_track_data[i].Write(f);
    }
    
    fclose(f);
}

void BaseApp::WriteTracks(char *filename) 
{
    FILE *f = fopen(filename, "w");

    if (f == NULL) {
	printf("Error opening file %s for writing\n", filename);
	return;
    }

    int num_images = GetNumImages();
    int num_tracks = (int) m_track_data.size();
    
    fprintf(f, "%d %d\n", num_images, num_tracks);

    /* Write the tracks */
    for (int i = 0; i < num_tracks; i++) {
        int num_views = (int) m_track_data[i].m_views.size();
        fprintf(f, "%d ", num_views);

        for (int j = 0; j < num_views; j++) {
            int img = m_track_data[i].m_views[j].first;
            int key = m_track_data[i].m_views[j].second;

            fprintf(f, "%d %d ", img, key);
        }
        fprintf(f, "\n");
    }
    
    fclose(f);
}

void BaseApp::WriteTracks2(char *filename) 
{
    FILE *f = fopen(filename, "w");

    if (f == NULL) {
	printf("Error opening file %s for writing\n", filename);
	return;
    }

    int num_images = GetNumImages();
    int num_tracks = (int) m_track_data.size();

    /* Load keys for all images */
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted)
            m_image_data[i].LoadKeys(false);
    }

    fprintf(f, "%d %d\n", num_images, num_tracks);

    /* Write the tracks */
    for (int i = 0; i < num_tracks; i++) {
        int num_views = (int) m_track_data[i].m_views.size();
        int num_good_views = 0;

        for (int j = 0; j < num_views; j++) {
            int img = m_track_data[i].m_views[j].first;

            if (!m_image_data[img].m_camera.m_adjusted)
                continue;

            num_good_views++;
        }

        if (num_good_views < 2)
            continue;

        fprintf(f, "%d ", num_good_views);

        for (int j = 0; j < num_views; j++) {
            int img = m_track_data[i].m_views[j].first;
            int key = m_track_data[i].m_views[j].second;

            if (!m_image_data[img].m_camera.m_adjusted)
                continue;

            double x = m_image_data[img].m_keys[key].m_x;
            double y = m_image_data[img].m_keys[key].m_y;

            fprintf(f, "%d %0.6e %0.6e ", img, x, y);
        }
        fprintf(f, "\n");
    }
    
    fclose(f);
}


/* Setup data structures storing the points and lines visible to
 * each image */
void BaseApp::SetupImagePoints(int min_views) 
{
    /* Clear the views */
    int num_images = GetNumImages();
    
    for (int i = 0; i < num_images; i++) {
	m_image_data[i].m_visible_points.clear();
    }

    int num_points = (int) m_point_data.size();

    for (int i = 0; i < num_points; i++) {
	PointData &p = m_point_data[i];
	int num_views = (int) p.m_views.size();

	if (num_views < min_views)
	    continue;

	for (int j = 0; j < num_views; j++) {
	    int img = p.m_views[j].first;
            int key = p.m_views[j].second;
	    m_image_data[img].m_visible_points.push_back(i);
            m_image_data[img].m_visible_keys.push_back(key);
	}
    }
}


/* Fix cameras (reflection bug) */
void BaseApp::FixReflectionBug() 
{
    printf("[FixReflectionBug] Reflecting scene...\n");

    int num_images = GetNumImages();

    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted)
	    m_image_data[i].m_camera.Reflect();
    }

    int num_points = (int) m_point_data.size();
    for (int i = 0; i < num_points; i++) {
	m_point_data[i].m_pos[2] = -m_point_data[i].m_pos[2];
    }
}

/* Compute orientations for each image */
void BaseApp::ComputeImageRotations() 
{
    double center[3], up[3], x_axis[3], z_axis[3], scale;
    SetupSceneGroundPlane(center, up, x_axis, z_axis, scale);

    int num_images = GetNumImages();
    
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted) continue;
	
	/* Project the up vector into the image */
	double R[9];
	m_image_data[i].m_camera.GetPose(R);

	double x[2] = { 1.0, 0.0 };
	double y[2] = { 0.0, 1.0 };
	double up_camera[3];

	matrix_transpose_product(3, 3, 3, 1, R, up, up_camera);

	double x_dot, y_dot;
	matrix_product(1, 2, 2, 1, up_camera, x, &x_dot);
	matrix_product(1, 2, 2, 1, up_camera, y, &y_dot);

	int rot;
	if (fabs(x_dot) > fabs(y_dot)) {
            printf("[ComputeImageRotations] "
                   "Rotating image %d [90]\n", i);

	    if (x_dot > 0.0) {
		rot = 3;
	    } else {
		rot = 1;
	    }
	} else {
	    if (y_dot > 0.0) {
		rot = 0;
	    } else {
                printf("[ComputeImageRotations] "
                       "Rotating image %d [180]\n", i);

		rot = 2;
	    }
	}

	m_image_data[i].m_rotation = rot;
    }
}


/* Estimate the world axes based on the image orientations */
void BaseApp::EstimateAxes(double *xaxis, double *yaxis, double *zaxis)
{
    printf("[EstimateAxes] Estimating axes\n");

    double RTR[9] = { 0.0, 0.0, 0.0, 
		      0.0, 0.0, 0.0, 
		      0.0, 0.0, 0.0 };

    int num_images = GetNumImages();

    bool found = false;
    double ref_axis[3];

    double min_deg = 80.0;
    double dot_threshold = cos(DEG2RAD(min_deg));

    if (m_up_image != -1) {
	found = true;
	memcpy(ref_axis, m_image_data[m_up_image].m_camera.m_R + 3, 
	       sizeof(double) * 3);
    } else {
        /* Find a good up image */
        int max_inliers = 0;
        int best_image = -1;
        for (int i = 0; i < num_images; i++) {
            if (!m_image_data[i].m_camera.m_adjusted) continue;
     
            double Ri[9];
            m_image_data[i].m_camera.
                GetUprightRotation(m_image_data[i].m_rotation, Ri);

            memcpy(ref_axis, Ri + 3, sizeof(double) * 3);

            int inliers = 0;
            for (int j = 0; j < num_images; j++) {
                if (!m_image_data[j].m_camera.m_adjusted) continue;
                if (i == j) continue;
            
                double Rj[9];
                
                m_image_data[j].m_camera.
                    GetUprightRotation(m_image_data[j].m_rotation, Rj);

                /* How orthogonal is Rj's xaxis to Ri's yaxis? */
                double dot;
                matrix_product(1, 3, 3, 1, ref_axis, Rj + 0, &dot);

                /* Within 10 degrees of 90? */
                if (fabs(dot) > dot_threshold)
                    continue;
                
                inliers++;
            }

            if (inliers > max_inliers) {
                max_inliers = inliers;
                best_image = i;
            }
        }

        printf("  best image: %d with %d inliers\n", best_image, max_inliers);

        double Rbest[9];
        m_image_data[best_image].m_camera.
            GetUprightRotation(m_image_data[best_image].m_rotation, Rbest);
        
        memcpy(ref_axis, Rbest + 3, sizeof(double) * 3);
    }
    
    /* Compute the x-axis moment matrix */
    std::vector<int> agree;
    
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted) continue;
		
	double R[9];
        m_image_data[i].m_camera.
            GetUprightRotation(m_image_data[i].m_rotation, R);

        double dot;
        matrix_product(1, 3, 3, 1, R + 0, ref_axis, &dot);
		    
        if (fabs(dot) > dot_threshold) {
            printf("  outlier: %d\n", i);
            continue;
        }
        
	agree.push_back(i);

        double T[9];
	matrix_transpose_product(1, 3, 1, 3, R, R, T);
	matrix_sum(3, 3, 3, 3, RTR, T, RTR);
    }

    /* Find the smallest eigenvalue/vector of RTR.  This is our y-axis */
    matrix_minimum_unit_norm_solution(3, 3, RTR, yaxis);

    /* Check if we should flip the axis */
    int num_agree = (int) agree.size();
    int num_pos = 0, num_neg = 0;

    for (int i = 0; i < num_agree; i++) {
	double dot;
	matrix_product(1, 3, 3, 1, m_image_data[agree[i]].m_camera.m_R + 3,
		       yaxis, &dot);	

	if (dot < -0.707106781186548)
	    num_neg++;
	else if (dot > 0.707106781186548)
	    num_pos++;

        double R[9];
        m_image_data[agree[i]].m_camera.
            GetUprightRotation(m_image_data[agree[i]].m_rotation, R);

	matrix_product(1, 3, 3, 1, R, yaxis, &dot);

        // double angle = acos(dot);
        // printf("  Twist angle[%d] = %0.3f\n", agree[i], RAD2DEG(angle));
        // fflush(stdout);
    }

    if (num_neg > num_pos) {
	printf("[EstimateAxes] Flipping y-axis (%d, %d)\n",
	       num_pos, num_neg);
	matrix_scale(3, 1, yaxis, -1.0, yaxis);
    } else {
	printf("[EstimateAxes] Not flipping\n");
    }

#if 0
    zaxis[0] = yaxis[1];
    zaxis[1] = -yaxis[0];
    zaxis[2] = 0.0;
#endif

    /* Compute the average z-axis */
    zaxis[0] = zaxis[1] = zaxis[2] = 0.0;

    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted) continue;

        double *zi = m_image_data[i].m_camera.m_R + 6;
        matrix_sum(1, 3, 1, 3, zaxis, zi, zaxis);
    }

    matrix_print(1, 3, zaxis);

    double dot;
    matrix_product(1, 3, 3, 1, yaxis, zaxis, &dot);
    printf("dot = %0.3f\n", dot / matrix_norm(3, 1, zaxis));

    matrix_cross(yaxis, zaxis, xaxis);

    double norm = matrix_norm(3, 1, xaxis);
    matrix_scale(3, 1, xaxis, 1.0 / norm, xaxis);

    matrix_cross(xaxis, yaxis, zaxis);
}


/* Analyze the scene and return various attributes */
void BaseApp::SetupSceneGroundPlane(double *center, double *up, 
                                    double *x_axis, double *z_axis,
                                    double &scale) 
{
    int num_images = GetNumOriginalImages();

    /* Estimate a plane through all the camera centers */
    int num_cameras = 0;
    
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted)
	    num_cameras++;
    }

    v3_t *cc = new v3_t[num_cameras];  /* Camera centers */

    /* Marshall the camera centers */
    int count = 0;
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
	    double pos[3];
	    m_image_data[i].m_camera.GetPosition(pos);

	    if (count >= num_cameras)
		printf("error!\n");

	    cc[count] = v3_new(pos[0], pos[1], pos[2]);
	    count++;
	}
    }

    /* Compute mean and "scale" of the camera center locations */
    double mean[3] = { 0.0, 0.0, 0.0 };
    
    count = 0;
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
	    double pos[3];
	    m_image_data[i].m_camera.GetPosition(pos);

	    mean[0] += pos[0];
	    mean[1] += pos[1];
	    mean[2] += pos[2];

	    count++;
	}
    }

    mean[0] /= num_cameras;
    mean[1] /= num_cameras;
    mean[2] /= num_cameras;
    
    double distance = 0.0;
    count = 0;
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
	    double pos[3];
	    m_image_data[i].m_camera.GetPosition(pos);

	    double dx = mean[0] - pos[0];
	    double dy = mean[1] - pos[1];
	    double dz = mean[2] - pos[2];

	    distance += dx * dx + dy * dy + dz * dz;

	    if (count >= num_cameras)
		printf("error!\n");

	    Vx(cc[count]) -= mean[0];
	    Vy(cc[count]) -= mean[1];
	    Vz(cc[count]) -= mean[2];

	    count++;
	}
    }

    double rms_distance = sqrt(distance / num_cameras);

    /* Fit the plane */
    double plane[4];
    double U[9], S[3], VT[9];
    double norm;

    if (m_up_image == -1) {	
	int inliers_out;

	// fit_3D_plane_through_origin(num_cameras, cc, plane);
#ifndef WIN32
	double scale_factor = 0.05;  /* For Trevi */
	// double scale_factor = 0.005;  /* For Notre Dame */
#else
	double scale_factor = 0.05;
#endif

	fit_3D_plane_ortreg_ransac(num_cameras, cc, 
				   1024, scale_factor * rms_distance, 
				   &inliers_out, plane);

	printf("inliers: %d / %d\n", inliers_out, num_cameras);

	v3_svd(num_cameras, cc, U, S, VT);

	norm = matrix_norm(3, 1, plane);
	
	plane[0] /= norm;
	plane[1] /= norm;
	plane[2] /= norm;
    } else {
	if (!m_image_data[m_up_image].m_camera.m_adjusted) {
	    printf("[SetupScene] Error: user refered to an "
		   "unadjusted camera\n");
	    exit(1);
	}

	/* Multiply the y-axis by the camera matrix */
	double yaxis[3] = { 0.0, 1.0, 0.0 };
	matrix_transpose_product(3, 3, 3, 1, m_image_data[m_up_image].m_camera.m_R, 
		       yaxis, plane);

	norm = matrix_norm(3, 1, plane);
	plane[0] /= norm;
	plane[1] /= norm;
	plane[2] /= norm;

	/* Project all the camera centers onto the plane */
	for (int i = 0; i < num_cameras; i++) {
	    double pt[3] = { Vx(cc[i]), Vy(cc[i]), Vz(cc[i]) };
	    double pt_par[3], pt_perp[3];
	    double dot;
	    
	    matrix_product(1, 3, 3, 1, pt, plane, &dot);
	    matrix_scale(3, 1, plane, dot, pt_par);
	    matrix_diff(3, 1, 3, 1, pt, pt_par, pt_perp);

	    cc[i] = v3_new(pt_perp[0], pt_perp[1], pt_perp[2]);
	}

	v3_svd(num_cameras, cc, U, S, VT);
    }

    center[0] = mean[0];
    center[1] = mean[1];
    center[2] = mean[2];


    /* Find out which orientation of the plane agrees with the most up
     * vectors */

    int num_neg = 0;
    int num_pos = 0;
    double y_axis[3] = { 0.0, 1.0, 0.0 };
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;

	double R[9];
	m_image_data[i].m_camera.GetPose(R);
	double up_cam[3];
	matrix_product(3, 3, 3, 1, R, y_axis, up_cam);
	
	double dot;
	matrix_product(1, 3, 3, 1, up_cam, plane, &dot);
	
	if (fabs(dot) < 0.8) continue;

	if (dot < 0.0) 
	    num_neg++;
	else
	    num_pos++;
    }

    printf("num_pos = %d\n", num_pos);
    printf("num_neg = %d\n", num_neg);

    if (num_pos >= num_neg) {
	up[0] = plane[0];
	up[1] = plane[1];
	up[2] = plane[2];
    } else {
	up[0] = -plane[0];
	up[1] = -plane[1];
	up[2] = -plane[2];	
    }

#if 0
    x_axis[0] = Vx(cc[0]) + mean[0];
    x_axis[1] = Vy(cc[0]) + mean[1];
    x_axis[2] = Vz(cc[0]) + mean[2];
#else
    int perm[3];
    qsort_perm(3, S, perm);
    int mid_idx = perm[1];

    x_axis[0] = -VT[mid_idx * 3 + 0];
    x_axis[1] = -VT[mid_idx * 3 + 1];
    x_axis[2] = -VT[mid_idx * 3 + 2];
#endif

    norm = matrix_norm(3, 1, x_axis);
    matrix_scale(3, 1, x_axis, 1.0 / norm, x_axis);

    /* Project m_x so that it is perpendicular to m_up */
    double x_par[3], x_perp[3];
    double dot;

    matrix_product(1, 3, 3, 1, up, x_axis, &dot);
    // printf("x_dot is %0.3e\n", dot);

    matrix_scale(3, 1, up, dot, x_par);
    matrix_diff(3, 1, 3, 1, x_axis, x_par, x_perp);
    norm = matrix_norm(3, 1, x_perp);
    matrix_scale(3, 1, x_perp, 1.0 / norm, x_axis);

    matrix_cross(x_axis, up, z_axis);

    if (!m_metric)
	scale = m_scale * rms_distance;
    else
	scale = 1000.0; /* 1000 meters */

    delete [] cc;
}


/* Analyze the scene and return various attributes */
void BaseApp::SetupScene(double *center, double *up, 
                         double *x_axis, double *z_axis,
                         double &scale) 
{
    SetupSceneGroundPlane(center, up, x_axis, z_axis, scale);

    if (m_estimate_up_vector_szeliski) {
	/* Recompute the axes with Rick's method */
	EstimateAxes(x_axis, up, z_axis);
    }
}


void BaseApp::TransformWorldReal()
{
#ifndef __DEMO__
    /* Transform the points */
    int num_points = (int) m_point_data.size();

    for (int i = 0; i < num_points; i++) {
	double *pos = m_point_data[i].m_pos;
	double p[4] = { pos[0], pos[1], pos[2], 1.0 }, Tp[4];

	matrix_product(4, 4, 4, 1, m_xform, p, Tp);
	memcpy(pos, Tp, 3 * sizeof(double));
    }

    /* Transform the cameras */
    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;

	double pose[9], pos[4];
	
	m_image_data[i].m_camera.GetPosition(pos);
	m_image_data[i].m_camera.GetPose(pose);
	
	pos[3] = 1.0;
	
	double M3x3[9];
	memcpy(M3x3 + 0, m_xform + 0, 3 * sizeof(double));
	memcpy(M3x3 + 3, m_xform + 4, 3 * sizeof(double));
	memcpy(M3x3 + 6, m_xform + 8, 3 * sizeof(double));
        double scale;
        matrix_product(1, 3, 3, 1, M3x3, M3x3, &scale);
        matrix_scale(3, 3, M3x3, 1.0 / sqrt(scale), M3x3);

	double pose_new[9], pos_new[4];
	matrix_product(3, 3, 3, 3, M3x3, pose, pose_new);
	matrix_product(4, 4, 4, 1, m_xform, pos, pos_new);

	/* Factor out the scaling */
	double mag = matrix_norm(3, 1, pose_new);
	matrix_scale(3, 3, pose_new, 1.0 / mag, pose_new);

	m_image_data[i].m_camera.SetPose(pose_new);
	m_image_data[i].m_camera.SetPosition(pos_new);
	m_image_data[i].m_camera.Finalize();
    }
#endif /* __DEMO__ */
}

void BaseApp::WriteCameras(const char *filename)
{
    FILE *f=fopen(filename,"w");
    if (f == NULL) {
        printf("[WriteCameras] Error opening file %s for writing\n", filename);
        return;
    }
    
    for(int i=0;i<GetNumImages();i++) m_image_data[i].m_camera.WriteFile(f);
    fclose(f);
}

void BaseApp::WritePoints(const char *filename)
{
    FILE *f=fopen(filename,"w");
    for(int i=0;i<(int)m_point_data.size();i++) 
        m_point_data[i].WriteCoordinates(f);
    fclose(f);
}

void BaseApp::RepositionScene(double *center_out, double *R_out, 
                              double &scale_out)
{
    printf("[RepositionScene] Repositioning scene\n");

    double center[3], up[3], x_axis[3], z_axis[3];
    double scale;

    SetupScene(center, up, x_axis, z_axis, scale);

    double R[9] = { x_axis[0], x_axis[1], x_axis[2],
		    up[0], up[1], up[2],
		    z_axis[0], z_axis[1], z_axis[2] };

    memcpy(center_out, center, 3 * sizeof(double));
    memcpy(R_out, R, 9 * sizeof(double));
    scale_out = scale;

    printf("center = [%0.3f, %0.3f, %0.3f]\n", 
	   center[0], center[1], center[2]);

    printf("up = [%0.3f, %0.3f, %0.3f]\n", up[0], up[1], up[2]);
    printf("x_axis = [%0.3f, %0.3f, %0.3f]\n", 
	   x_axis[0], x_axis[1], x_axis[2]);
    printf("z_axis = [%0.3f, %0.3f, %0.3f]\n", 
	   z_axis[0], z_axis[1], z_axis[2]);

    printf("scale = %0.3f\n", scale);

    int num_points = (int) m_point_data.size();

    for (int i = 0; i < num_points; i++) {
	matrix_diff(3, 1, 3, 1, m_point_data[i].m_pos, center, 
		    m_point_data[i].m_pos);
	matrix_scale(3, 1, m_point_data[i].m_pos, m_scale / scale, 
                     m_point_data[i].m_pos);

	double tmp[3];
	matrix_product(3, 3, 3, 1, R, m_point_data[i].m_pos, tmp);

	m_point_data[i].m_pos[0] = tmp[0];
	m_point_data[i].m_pos[1] = tmp[1];
	m_point_data[i].m_pos[2] = tmp[2];
    }


    double R16[16], t16[16], M[16];
    matrix_ident(4, R16);
    matrix_ident(4, t16);
    
    memcpy(R16 + 0, R + 0, 3 * sizeof(double));
    memcpy(R16 + 4, R + 3, 3 * sizeof(double));
    memcpy(R16 + 8, R + 6, 3 * sizeof(double));

    t16[3] = -center[0];
    t16[7] = -center[1];
    t16[11] = -center[2];

    matrix_product(4, 4, 4, 4, R16, t16, M);

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
	    double pos[3];
	    m_image_data[i].m_camera.GetPosition(pos);

	    matrix_diff(3, 1, 3, 1, pos, center, pos);
	    matrix_scale(3, 1, pos, m_scale / scale, pos);

            double norm = matrix_norm(3, 1, pos);
            // printf("Camera %d has distance %0.3f\n", i, norm);

	    double tmp[3];
	    matrix_product(3, 3, 3, 1, R, pos, tmp);

	    double pose[9];
	    m_image_data[i].m_camera.GetPose(pose);

	    double pose_new[9];
	    matrix_product(3, 3, 3, 3, R, pose, pose_new);

	    m_image_data[i].m_camera.SetPose(pose_new);
	    m_image_data[i].m_camera.SetPosition(tmp);
	    m_image_data[i].m_camera.Finalize();

            /* Reposition the plane */
            m_image_data[i].m_fit_plane.Transform(M);


#if 0
            double Ri[9];
            m_image_data[i].m_camera.GetPose(Ri);

            // camera x_axis: u, y_axis: v, z_axis (viewing direction): w
            // Twist = u_x * w_z - u_z * w_x / sqrt(1 - w_y^2)
            // (twist is undefined if camera is pointed directly up or down)
            double c_twist = 
                (Ri[0] * Ri[8] - Ri[6] * Ri[2]) / sqrt(1 - Ri[5] * Ri[5]);
                // Ri[0] * Ri[6] - Ri[3] * Ri[8] / sqrt(1 - Ri[7] * Ri[7]);

            double twist = acos(c_twist);
            printf("  Real twist angle[%d] = %0.3f\n", i, RAD2DEG(twist));
            matrix_print(3, 3, Ri);

            fflush(stdout);
#endif
	}
    }

    memcpy(m_repos_R, R, sizeof(double) * 9);
    memcpy(m_repos_d, center, sizeof(double) * 3);
    m_repos_scale = m_scale / scale;

#if 0
    FILE *f = fopen("overhead.txt", "w");
    
    if (f == NULL)
        return;
    
    fprintf(f, "%d %d\n", num_images, num_points);
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
            double pos[3];
            m_image_data[i].m_camera.GetPosition(pos);
            fprintf(f, "%0.6e %0.6e\n", pos[0], pos[2]);
        } else {
            fprintf(f, "0.0 0.0\n");
        }
    }

    for (int i = 0; i < num_points; i++) {
        double *pos = m_point_data[i].m_pos;
        fprintf(f, "%0.6e %0.6e\n", pos[0], pos[2]);
    }

    fclose(f);
#endif
}

void BaseApp::TransformSceneCanonical(int i1, int i2)
{
    double R[9];
    memcpy(R, m_image_data[i1].m_camera.m_R, sizeof(double) * 9);
    
    double origin[3];
    m_image_data[i1].m_camera.GetPosition(origin);

    double eye2[3], diff[3];
    m_image_data[i2].m_camera.GetPosition(eye2);
    matrix_diff(3, 1, 3, 1, eye2, origin, diff);
    
    double scale = 1.0 / matrix_norm(3, 1, diff);

    double proj[2];
    m_image_data[1].m_camera.Project(m_point_data[0].m_pos, proj);
    printf("Proj[before] = %0.3f, %0.3f\n", proj[0], proj[1]);

    /* Transform points */
    int num_points = (int) m_point_data.size();
    for (int i = 0; i < num_points; i++) {
	matrix_diff(3, 1, 3, 1, m_point_data[i].m_pos, origin, 
		    m_point_data[i].m_pos);

	double tmp[3];
	matrix_product(3, 3, 3, 1, R, m_point_data[i].m_pos, tmp);

	m_point_data[i].m_pos[0] = tmp[0];
	m_point_data[i].m_pos[1] = tmp[1];
	m_point_data[i].m_pos[2] = tmp[2];

         matrix_scale(3, 1, 
                      m_point_data[i].m_pos, scale, m_point_data[i].m_pos);
    }

    double R16[16], t16[16], M[16];
    matrix_ident(4, R16);
    matrix_ident(4, t16);
    
    memcpy(R16 + 0, R + 0, 3 * sizeof(double));
    memcpy(R16 + 4, R + 3, 3 * sizeof(double));
    memcpy(R16 + 8, R + 6, 3 * sizeof(double));

    t16[3] = -origin[0];
    t16[7] = -origin[1];
    t16[11] = origin[2];

    matrix_product(4, 4, 4, 4, R16, t16, M);

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
	    double pos[3], pos_new[3];
	    m_image_data[i].m_camera.GetPosition(pos);

	    matrix_diff(3, 1, 3, 1, pos, origin, pos);
	    matrix_product(3, 3, 3, 1, R, pos, pos_new);
            matrix_scale(3, 1, pos_new, scale, pos_new);

#if 0
	    double pose[9], pose_new[9];
	    m_image_data[i].m_camera.GetPose(pose);
	    matrix_product(3, 3, 3, 3, R, pose, pose_new);

	    m_image_data[i].m_camera.SetPose(pose_new);
	    memcpy(m_image_data[i].m_camera.m_R, pose_new, 9 * sizeof(double));
#else
            double *R2 = m_image_data[i].m_camera.m_R;
            double R_new[9];
            matrix_transpose_product2(3, 3, 3, 3, R2, R, R_new);
	    memcpy(m_image_data[i].m_camera.m_R, R_new, 9 * sizeof(double));
#endif

	    m_image_data[i].m_camera.SetPosition(pos_new);
	    m_image_data[i].m_camera.Finalize();

            /* Reposition the plane */
            m_image_data[i].m_fit_plane.Transform(M);
	}
    }

    m_image_data[1].m_camera.Project(m_point_data[0].m_pos, proj);
    printf("Proj[after] = %0.3f, %0.3f\n", proj[0], proj[1]);

#if 0
    FILE *f = fopen("overhead.txt", "w");
    
    if (f == NULL)
        return;
    
    fprintf(f, "%d %d\n", num_images, num_points);
    for (int i = 0; i < num_images; i++) {
	if (m_image_data[i].m_camera.m_adjusted) {
            double pos[3];
            m_image_data[i].m_camera.GetPosition(pos);
            fprintf(f, "%0.6e %0.6e\n", pos[0], pos[2]);
        } else {
            fprintf(f, "0.0 0.0\n");
        }
    }

    for (int i = 0; i < num_points; i++) {
        double *pos = m_point_data[i].m_pos;
        fprintf(f, "%0.6e %0.6e\n", pos[0], pos[2]);
    }

    fclose(f);
#endif
}

void BaseApp::UnscaleCameras(int start_camera)
{

#if 0
    double center[3], up[3], x_axis[3], z_axis[3];
    double scale;

    SetupScene(center, up, x_axis, z_axis, scale);
#endif

    int num_images = GetNumImages();
    for (int i = start_camera; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted) {
            printf("[UnscaleCameras] Unscaling camera %d (scale: %0.3f...)\n", 
                i, m_repos_scale);

            double pos[3];
            m_image_data[i].m_camera.GetPosition(pos);

            matrix_diff(3, 1, 3, 1, pos, m_repos_d, pos);
            matrix_scale(3, 1, pos, m_repos_scale, pos);

            double norm = matrix_norm(3, 1, pos);
            // printf("Camera %d has distance %0.3f\n", i, norm);

            double tmp[3];
            matrix_product(3, 3, 3, 1, m_repos_R, pos, tmp);

            double pose[9];
            m_image_data[i].m_camera.GetPose(pose);

            double pose_new[9];
            matrix_product(3, 3, 3, 3, m_repos_R, pose, pose_new);

            m_image_data[i].m_camera.SetPose(pose_new);
            m_image_data[i].m_camera.SetPosition(tmp);
            m_image_data[i].m_camera.Finalize();
        }
    }    
}

#ifdef __USE_ANN__
ANNkd_tree *Create3DSearchTree(int n, v3_t *v)
{
    ANNpointArray pts = annAllocPts(n, 3);

    for (int i = 0; i < n; i++) {
        pts[i][0] = Vx(v[i]);
        pts[i][1] = Vy(v[i]);
        pts[i][2] = Vz(v[i]);
    }

    /* Create a search tree for k2 */
    ANNkd_tree *tree = new ANNkd_tree(pts, n, 3, 4);

    return tree;
}
#endif /* __USE_ANN__ */

/* Create a search tree for camera centers */
#ifdef __USE_ANN__
ANNkd_tree *BaseApp::CreateCameraSearchTree()
#else
BruteForceSearch *BaseApp::CreateCameraSearchTree()
#endif
{
    int num_images = GetNumImages();

    /* Find the cameras */
    v3_t *cameras = new v3_t[num_images];
    int *index_map = new int[num_images];

    int num_cameras = 0;
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;
	
	double pos[3];
	m_image_data[i].m_camera.GetPosition(pos);

	cameras[num_cameras] = v3_new(pos[0], pos[1], pos[2]);
	index_map[num_cameras] = i;

	num_cameras++;
    }

    /* Create a search tree */
#ifdef __USE_ANN__
    ANNkd_tree *tree = Create3DSearchTree(num_cameras, cameras);
#else
    BruteForceSearch *search = new BruteForceSearch(num_cameras, cameras);
#endif

    delete [] cameras;
    delete [] index_map;

#ifdef __USE_ANN__
    return tree;
#else
    return search;
#endif
}

static double ComputeConfidence(int num_rays, const v3_t *rays)
{
    if (num_rays <= 2)
        return 0.0;

    /* Compute the average ray */
    v3_t avg = v3_unit(v3_mean(num_rays, rays));
    
    /* Compute the ray furthest from the average */
    v3_t ex1 = v3_extremum(num_rays, rays, avg);
    
    /* Compute the ray furthest from the first */
    v3_t ex2 = v3_extremum(num_rays, rays, ex1);

    /* Compute the ray furthest from the first two */
    v3_t ex3 = v3_extremum2(num_rays, rays, ex1, ex2);

    /* Compute the minimum dot product between the rays */
    double dot12 = v3_dotp(ex1, ex2);
    double dot23 = v3_dotp(ex2, ex3);
    double dot13 = v3_dotp(ex1, ex3);
    
    double max_dot = MAX(MAX(dot12, dot23), dot13);
    double angle = acos(max_dot);
    
    const double max_angle = 20.0;

    return CLAMP(RAD2DEG(angle) / max_angle, 0.0, 1.0);
}

void BaseApp::EstimatePointNormalsConfidence()
{
    int num_points = (int) m_point_data.size();
    
    for (int i = 0; i < num_points; i++) {
        PointData &p = m_point_data[i];
        
        int num_views = (int) p.m_views.size();
        v3_t *rays = new v3_t[num_views];
        
        /* Estimate normal */
        double normal[3] = { 0.0, 0.0, 0.0 };

        for (int j = 0; j < num_views; j++) {
            int cam = p.m_views[j].first;
            
            double cam_pos[3];
            double *ray = rays[j].p;
            m_image_data[cam].m_camera.GetPosition(cam_pos);
            matrix_diff(3, 1, 3, 1, p.m_pos, cam_pos, ray);
            
            double norm = matrix_norm(3, 1, ray);
            matrix_scale(3, 1, ray, 1.0 / norm, ray);

            matrix_sum(3, 1, 3, 1, normal, ray, normal);
        }

        double norm = matrix_norm(3, 1, normal);
        matrix_scale(3, 1, normal, -1.0 / norm, p.m_norm);

        /* Compute the confidence in the point */
        p.m_conf = ComputeConfidence(num_views, rays);

        printf("conf[%d] = %0.3f\n", i, p.m_conf);

        delete [] rays;
    }
}


void BaseApp::EstimatePointNormals() 
{
#define NUM_NNS 32 // 100
    int num_points = (int) m_point_data.size();

    v3_t *pts = new v3_t[num_points];

    for (int i = 0; i < num_points; i++) {
        pts[i] = v3_new(m_point_data[i].m_pos[0],
            m_point_data[i].m_pos[1],
            m_point_data[i].m_pos[2]);

        m_point_data[i].m_norm[0] = 0.0;
        m_point_data[i].m_norm[1] = 0.0;
        m_point_data[i].m_norm[2] = 0.0;
    }

#ifdef __USE_ANN__
    ANNkd_tree *tree = Create3DSearchTree(num_points, pts);
#else
    BruteForceSearch *search = new BruteForceSearch(num_points, pts);
#endif

    int num_points_inc = 0;
    for (int i = 0; i < num_points; i++) {
        if (i % 500 == 0) {
            printf(".");
            fflush(stdout);
        }

        /* For each point, find the NUM_NNS nearest neighbors */
        PointData &p = m_point_data[i];
        v3_t q = v3_new(p.m_pos[0], p.m_pos[1], p.m_pos[2]);

        int nn_idxs[NUM_NNS];

#ifdef __USE_ANN__
        float dists[NUM_NNS];
        float query[3] = { Vx(q), Vy(q), Vz(q) };
        tree->annkPriSearch(query, NUM_NNS, nn_idxs, dists, 0.0);
#else
        double dists[NUM_NNS];
        search->GetClosestPoints(q, NUM_NNS, nn_idxs, dists);
#endif

        v3_t nns[NUM_NNS];

        for (int j = 0; j < NUM_NNS; j++) {
            PointData &nn = m_point_data[nn_idxs[j]];
            nns[j] = v3_new(nn.m_pos[0], nn.m_pos[1], nn.m_pos[2]);
        }

        /* Fit a plane to the nns */
        double params[4];

        int num_inliers;
        // double error = 
        fit_3D_plane_ortreg_ransac(NUM_NNS, nns, 64, 0.1, 
            &num_inliers, params);

#if 0
#define INLIER_THRESHOLD (0.7 * NUM_NNS)
        if (num_inliers < INLIER_THRESHOLD)
            continue;
#endif

        // printf("inliers = %d\n", num_inliers);

#if 0
#define ERROR_THRESHOLD 0.1
        if (error > ERROR_THRESHOLD)
            continue;
#endif

        /* Project the normal onto the x-z plane */
        double normal[3] = { params[0], params[1], params[2] };

        int ref_image = -1;
        double max_dot = 0.0;
        double max_dot_signed = 0.0;
        int num_views = (int) m_point_data[i].m_views.size();

        for (int j = 0; j < num_views; j++) {
            int idx = m_point_data[i].m_views[j].first;

            double dir[3];
            m_image_data[idx].m_camera.GetViewDirection(dir);

            double dot;
            matrix_product(1, 3, 3, 1, normal, dir, &dot);

            if (fabs(dot) > max_dot) {
                max_dot = fabs(dot);
                max_dot_signed = dot;
                ref_image = idx;
            }
        }

        if (ref_image == -1) {
            printf("ref image = -1\n");
        }

        if (max_dot_signed < 0.0) {
            matrix_scale(3, 1, normal, -1.0, normal);
        }

        memcpy(m_point_data[i].m_norm, normal, 3 * sizeof(double));
        m_point_data[i].m_ref_image = ref_image;

#if 0
        double dot;
        matrix_product(1, 3, 3, 1, normal, up, &dot);

        double parallel[3];
        matrix_scale(3, 1, up, dot, parallel);

        double perp[3];
        matrix_diff(3, 1, 3, 1, normal, parallel, perp);

#define NORM_THRESHOLD 0.95 // 0.9
        double norm = matrix_norm(3, 1, perp);

        if (norm < NORM_THRESHOLD)
            continue;

        matrix_scale(3, 1, perp, 1.0 / norm, perp);

        /* Compute the angle of the new normal */
        matrix_product(1, 3, 3, 1, perp, xaxis, &dot);
        double theta = acos(dot);

        double cross[3];
        matrix_cross(perp, xaxis, cross);
        matrix_product(1, 3, 3, 1, cross, up, &dot);	
#endif

        num_points_inc++;
    }

    printf("\n");

    delete [] pts;

#ifdef __USE_ANN__
    delete tree;
#else
    delete search;
#endif
}


void BaseApp::RemoveBadImages(int min_num_points)
{
    int num_images = GetNumImages();

    for (int i = 0; i < num_images; i++) {
        if (!m_image_data[i].m_camera.m_adjusted)
            continue;

        int num_visible = (int) m_image_data[i].m_visible_points.size();

        if (num_visible < min_num_points) {
            printf("[RemoveBadImages] Removing image %d (%d points visible)\n",
                   i, num_visible);
            m_image_data[i].m_camera.m_adjusted = false;

            for (int j = 0; j < num_visible; j++) {
                int idx = m_image_data[i].m_visible_points[j];
                
                PointData &pt = m_point_data[idx];
                int num_views = pt.m_views.size();
                
                for (int k = 0; k < num_views; k++) {
                    if (pt.m_views[k].first == i) {
                        printf("  Erasing from point %d\n", k);
                        pt.m_views.erase(pt.m_views.begin() + k);
                        break;
                    }
                }
            }
        }
    }
}

bool BaseApp::ImagesPartOfPanorama(int i1, int i2) 
{
    if (!m_image_data[i1].m_camera.m_adjusted || 
        !m_image_data[i2].m_camera.m_adjusted)
        return false;

    double pos1[3], pos2[3];
    m_image_data[i1].m_camera.GetPosition(pos1);
    m_image_data[i2].m_camera.GetPosition(pos2);

    /* The images have to have at least one point in common */
    std::vector<int> isect = 
        GetVectorIntersection(m_image_data[i1].m_visible_points,
                              m_image_data[i2].m_visible_points);

    if (isect.size() == 0)
        return false;

    /* Over each point seen by either camera, compute the average
     * angle between the rays */ 

    unsigned int num_vis1 = m_image_data[i1].m_visible_points.size();
    unsigned int num_vis2 = m_image_data[i2].m_visible_points.size();
    assert(num_vis1 + num_vis2 != 0);

    double angle_sum = 0.0, dist1_sum = 0.0, dist2_sum = 0.0;
    for (unsigned int i = 0; i < num_vis1; i++) {
        int p = m_image_data[i1].m_visible_points[i];
        double ray1[3], ray2[3];
        matrix_diff(3, 1, 3, 1, 
                    m_point_data[p].m_pos, pos1, ray1);
        matrix_diff(3, 1, 3, 1, 
                    m_point_data[p].m_pos, pos2, ray2);
        
        double dot;
        matrix_product(1, 3, 3, 1, ray1, ray2, &dot);

        double dist1 = matrix_norm(3, 1, ray1);
        double dist2 = matrix_norm(3, 1, ray2);

        double mag = dist1 * dist2;
        double angle = acos(CLAMP(dot / mag, -1.0 + 1.0e-8, 1.0 - 1.0e-8));
        angle_sum += angle;

        dist1_sum += dist1;
        dist2_sum += dist2;
    }

    for (unsigned int i = 0; i < num_vis2; i++) {
        int p = m_image_data[i2].m_visible_points[i];
        double ray1[3], ray2[3];
        matrix_diff(3, 1, 3, 1, 
                    m_point_data[p].m_pos, pos1, ray1);
        matrix_diff(3, 1, 3, 1, 
                    m_point_data[p].m_pos, pos2, ray2);
        
        double dot;
        matrix_product(1, 3, 3, 1, ray1, ray2, &dot);

        double dist1 = matrix_norm(3, 1, ray1);
        double dist2 = matrix_norm(3, 1, ray2);

        double mag = dist1 * dist2;
        double angle = acos(CLAMP(dot / mag, -1.0 + 1.0e-8, 1.0 - 1.0e-8));
        angle_sum += angle;

        dist1_sum += dist1;
        dist2_sum += dist2;
    }

    double angle_avg = angle_sum / (num_vis1 + num_vis2);
    double dist1_avg = dist1_sum / (num_vis1 + num_vis2);
    double dist2_avg = dist2_sum / (num_vis1 + num_vis2);
    
    double disp_cams[3];
    matrix_diff(3, 1, 3, 1, pos1, pos2, disp_cams);
    
    double dist_cams = matrix_norm(3, 1, disp_cams);

    // printf("  angle_avg[%d,%d] = %0.3f\n", i1, i2, RAD2DEG(angle_avg));
    
    if (RAD2DEG(angle_avg) > 3.0 /* 2.5 */ /* 1.5 */)
        return false;
    
    if (dist_cams > 0.1 * dist1_avg)
        return false;
    
    if (dist_cams > 0.1 * dist2_avg)
        return false;
    
    return true;
}
