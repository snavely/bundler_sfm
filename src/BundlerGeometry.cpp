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

/* BundlerGeometry.cpp */
/* Basic routines for computing image geometry */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "BundlerApp.h"
#include "BundleUtil.h"
#include "Epipolar.h"
#include "Register.h"

#include "defines.h"
#include "horn.h"
#include "matrix.h"
#include "qsort.h"
#include "util.h"



double global_scale = 1.0;

int BundlerApp::GetNumCameraParameters() {
    int cnp = 6;
     
    if (!m_fixed_focal_length)
        cnp++;
    
    if (m_estimate_distortion)
        cnp += 2;

    return cnp;
}

/* File IO routines */
void TransformInfo::ReadFromFile(FILE *f) {
    /* Homography */
    fscanf(f, "%lf %lf %lf "
	      "%lf %lf %lf "
	      "%lf %lf %lf",
	   m_H + 0, m_H + 1, m_H + 2,
	   m_H + 3, m_H + 4, m_H + 5,
	   m_H + 6, m_H + 7, m_H + 8);

    /* F-Matrix */
    fscanf(f, "%lf %lf %lf "
	      "%lf %lf %lf "
	      "%lf %lf %lf",
	   m_fmatrix + 0, m_fmatrix + 1, m_fmatrix + 2,
	   m_fmatrix + 3, m_fmatrix + 4, m_fmatrix + 5,
	   m_fmatrix + 6, m_fmatrix + 7, m_fmatrix + 8);

    /* Inlier info */
    fscanf(f, "%lf\n", &m_inlier_ratio);
    fscanf(f, "%d\n", &m_num_inliers);    
}

void TransformInfo::WriteToFile(FILE *f) {
    /* Homography */
    fprintf(f, "%0.6le %0.6le %0.6le "
	       "%0.6le %0.6le %0.6le "
	       "%0.6le %0.6le %0.6le\n",
	    m_H[0], m_H[1], m_H[2],
	    m_H[3], m_H[4], m_H[5],
	    m_H[6], m_H[7], m_H[8]);
    
    /* F-Matrix */
    fprintf(f, "%0.6le %0.6le %0.6le "
	       "%0.6le %0.6le %0.6le "
	       "%0.6le %0.6le %0.6le\n",
	    m_fmatrix[0], m_fmatrix[1], m_fmatrix[2],
	    m_fmatrix[3], m_fmatrix[4], m_fmatrix[5],
	    m_fmatrix[6], m_fmatrix[7], m_fmatrix[8]);

    /* Inlier info */
    fprintf(f, "%0.16le\n", m_inlier_ratio);
    fprintf(f, "%d\n", m_num_inliers);
}

void BundlerApp::ComputeGeometricConstraints(bool overwrite, 
                                             int new_image_start) 
{
    int num_images = GetNumImages();

    /* Read information from files if they exist */
    const char *filename = "constraints.txt";
    if (!overwrite && FileExists(filename)) {
	ReadGeometricConstraints(filename);
        return;
    } else {
        LoadMatches();

        if (num_images < 40000) 
            WriteMatchTable(".prune");

        if (!m_skip_fmatrix || !m_skip_homographies || 
            m_keypoint_border_width > 0 || m_keypoint_border_bottom > 0)
            LoadKeys(false);

        if (m_keypoint_border_width > 0) {
            for (int i = 0; i < num_images; i++) {
                for (int j = i+1; j < num_images; j++) {
                    if (!ImagesMatch(i, j))
                        continue;

                    RemoveMatchesNearBorder(i, j, m_keypoint_border_width);
                }
            }
        }

        if (m_keypoint_border_bottom > 0) {
            for (int i = 0; i < num_images; i++) {
                for (int j = i+1; j < num_images; j++) {
                    if (!ImagesMatch(i, j))
                        continue;

                    RemoveMatchesNearBottom(i, j, m_keypoint_border_bottom);
                }
            }
        }

        if (!m_skip_fmatrix) {
            ComputeEpipolarGeometry(true, new_image_start);
        }

        if (!m_skip_homographies) {
            ComputeTransforms(false, new_image_start);
        }

	MakeMatchListsSymmetric();

        if (num_images < 40000)
            WriteMatchTable(".ransac");

        // RemoveAllMatches();
	ComputeTracks(new_image_start);

        // ClearMatches();
        RemoveAllMatches();
        // SetMatchesFromTracks();

#if 1
        /* Set match flags */
        int num_tracks = (int) m_track_data.size();
        for (int i = 0; i < num_tracks; i++) {
            TrackData &track = m_track_data[i];
            int num_views = (int) track.m_views.size();

            for (int j = 0; j < num_views; j++) {
                int img1 = track.m_views[j].first;

                assert(img1 >= 0 && img1 < num_images);

                for (int k = j+1; k < num_views; k++) {
                    int img2 = track.m_views[k].first;

                    assert(img2 >= 0 && img2 < num_images);
                    
                    SetMatch(img1, img2);
                    SetMatch(img2, img1);
                }
            }
        }
#endif

        WriteGeometricConstraints(filename);

        if (num_images < 40000)
            WriteMatchTable(".corresp");

	// ComputeMatchPoints(new_image_start);

        // FilterTracks();
    }
}

/* Compute a transform between a given pair of images */
bool BundlerApp::ComputeTransform(int idx1, int idx2, bool removeBadMatches)
{
    assert(m_image_data[idx1].m_keys_loaded);
    assert(m_image_data[idx2].m_keys_loaded);

    MatchIndex offset = GetMatchIndex(idx1, idx2);

    double M[9];

    if (idx1 == idx2) {
	printf("[ComputeTransform] Error: computing tranform "
	       "for identical images\n");
	return false;
    }
    
    std::vector<KeypointMatch> &list = m_matches.GetMatchList(offset);

    std::vector<int> inliers = 
	EstimateTransform(m_image_data[idx1].m_keys, 
			  m_image_data[idx2].m_keys, 
			  list, MotionHomography,
			  m_homography_rounds, 
			  m_homography_threshold,
			  /* 15.0 */ /* 6.0 */ /* 4.0 */ M);
    
    int num_inliers = (int) inliers.size();

    printf("Inliers[%d,%d] = %d out of %d\n", idx1, idx2, num_inliers, 
           (int) list.size());

    if (removeBadMatches) {
	/* Refine the matches */
	std::vector<KeypointMatch> new_match_list;
    
	for (int i = 0; i < num_inliers; i++) {
	    new_match_list.push_back(list[inliers[i]]);
	}

	// m_match_lists[offset].clear();
	// m_match_lists[offset] = new_match_list;
        list.clear();
        list = new_match_list;
    }

#define MIN_INLIERS 10
    if (num_inliers >= MIN_INLIERS) {
	m_transforms[offset].m_num_inliers = num_inliers;
	m_transforms[offset].m_inlier_ratio = 
	    ((double) num_inliers) / ((double) list.size());

	memcpy(m_transforms[offset].m_H, M, 9 * sizeof(double));
	// m_transforms[offset]->m_scale = sqrt(M[0] * M[0] + M[1] * M[1]);
#if 1
	printf("Inliers[%d,%d] = %d out of %d\n", idx1, idx2, num_inliers, 
               (int) m_matches.GetNumMatches(offset));
	       // (int) m_match_lists[offset].size());
	printf("Ratio[%d,%d] = %0.3e\n", 
	       idx1, idx2, m_transforms[offset].m_inlier_ratio);

	matrix_print(3, 3, m_transforms[offset].m_H);
#endif

	return true;
    } else {
	return false;
    }
}


/* Compute rigid transforms between all matching images */
void BundlerApp::ComputeTransforms(bool removeBadMatches, int new_image_start) 
{
    unsigned int num_images = GetNumImages();

    m_transforms.clear();

    for (unsigned int i = 0; i < num_images; i++) {
        MatchAdjList::iterator iter;
        for (iter = m_matches.Begin(i); iter != m_matches.End(i); iter++) {
            unsigned int j = iter->m_index;

            assert(ImagesMatch(i, j));

            MatchIndex idx = GetMatchIndex(i, j);
            MatchIndex idx_rev = GetMatchIndex(j, i);

            m_transforms[idx] = TransformInfo();
            m_transforms[idx_rev] = TransformInfo();

            bool connect12 = ComputeTransform(i, j, removeBadMatches);

            if (!connect12) {
                if (removeBadMatches) {
                    // RemoveMatch(i, j);
                    // RemoveMatch(j, i);

                    // m_match_lists[idx].clear();
                    m_matches.RemoveMatch(idx);
                    m_matches.RemoveMatch(idx_rev);
                    // m_match_lists.erase(idx);

                    m_transforms.erase(idx);
                    m_transforms.erase(idx_rev);
                }
            } else {
                matrix_invert(3, m_transforms[idx].m_H, 
                              m_transforms[idx_rev].m_H);
            }
        }
    }

    /* Print the inlier ratios */
    FILE *f = fopen("pairwise_scores.txt", "w");
    
    for (unsigned int i = 0; i < num_images; i++) {
        MatchAdjList::iterator iter;

        for (iter = m_matches.Begin(i); iter != m_matches.End(i); iter++) {
            unsigned int j = iter->m_index; // first;
 
            assert(ImagesMatch(i, j));

            // MatchIndex idx = *iter;
            MatchIndex idx = GetMatchIndex(i, j);
            fprintf(f, "%d %d %0.5f\n", i, j, 
                    m_transforms[idx].m_inlier_ratio);
        }
    }

    fclose(f);
}

/* Compute epipolar geometry between a given pair of images */
bool BundlerApp::ComputeEpipolarGeometry(int idx1, int idx2, 
                                         bool removeBadMatches) 
{
    assert(m_image_data[idx1].m_keys_loaded);
    assert(m_image_data[idx2].m_keys_loaded);

    MatchIndex offset = GetMatchIndex(idx1, idx2);
    MatchIndex offset_rev = GetMatchIndex(idx2, idx1);
    std::vector<KeypointMatch> &list = m_matches.GetMatchList(offset);

    double F[9];
    
    std::vector<int> inliers = 
	EstimateFMatrix(m_image_data[idx1].m_keys, 
			m_image_data[idx2].m_keys, 
			list,
			m_fmatrix_rounds, 
			m_fmatrix_threshold /* 20.0 */ /* 9.0 */, F);
    
    int num_inliers = (int) inliers.size();

    printf("Inliers[%d,%d] = %d out of %d\n", idx1, idx2, num_inliers, 
	   (int) list.size());

    if (removeBadMatches) {
	/* Refine the matches */
	std::vector<KeypointMatch> new_match_list;
    
	for (int i = 0; i < num_inliers; i++) {
	    new_match_list.push_back(list[inliers[i]]);
	}

	// m_match_lists[offset].clear();
	// m_match_lists[offset] = new_match_list;
        list.clear();
        list = new_match_list;
    }
    
#define MIN_INLIERS_EPIPOLAR 16
    if (num_inliers >= m_min_num_feat_matches /*MIN_INLIERS_EPIPOLAR*/) {
	// if (m_transforms[offset] == NULL) 
        m_transforms[offset] = TransformInfo();
        m_transforms[offset_rev] = TransformInfo();

	memcpy(m_transforms[offset].m_fmatrix, F, 9 * sizeof(double));
	// m_transforms[offset]->m_scale = sqrt(M[0] * M[0] + M[1] * M[1]);
	printf("Inliers[%d,%d] = %d out of %lu\n", idx1, idx2, num_inliers, 
               list.size());

	return true;
    } else {
	return false;
    }
}

/* Compute epipolar geometry between all matching images */
void BundlerApp::ComputeEpipolarGeometry(bool removeBadMatches, 
                                         int new_image_start) 
{
    unsigned int num_images = GetNumImages();
    
    m_transforms.clear();

    std::vector<MatchIndex> remove;

    for (unsigned int i = 0; i < num_images; i++) {
        MatchAdjList::iterator iter;

        for (iter = m_matches.Begin(i); iter != m_matches.End(i); iter++) {
            unsigned int j = iter->m_index; // first;

            assert(ImagesMatch(i, j));

            MatchIndex idx = GetMatchIndex(i, j);
            MatchIndex idx_rev = GetMatchIndex(j, i);

            bool connect12 = 
                ComputeEpipolarGeometry(i, j, removeBadMatches);

            if (!connect12) {
                if (removeBadMatches) {
                    // RemoveMatch(i, j);
                    // RemoveMatch(j, i);
                    remove.push_back(idx);
                    remove.push_back(idx_rev);

                    // m_match_lists[idx].clear();
                    // m_match_lists.erase(idx);

                    m_transforms.erase(idx);
                    m_transforms.erase(idx_rev);
                }
            } else {
                matrix_transpose(3, 3, 
                                 m_transforms[idx].m_fmatrix, 
                                 m_transforms[idx_rev].m_fmatrix);
            }
        }
    }

    int num_removed = (int) remove.size();
    
    for (int i = 0; i < num_removed; i++) {
        int img1 = remove[i].first;
        int img2 = remove[i].second;

        // RemoveMatch(img1, img2);
        m_matches.RemoveMatch(GetMatchIndex(img1, img2));
    }
}


/* Coalesce feature descriptors for each feature point */
void BundlerApp::CoalesceFeatureDescriptorsMedian()
{
    int num_points = (int) m_point_data.size();

    for (int i = 0; i < num_points; i++) {
	m_point_data[i].m_desc = new float[128];
	for (int j = 0; j < 128; j++) 
	    m_point_data[i].m_desc[j] = 0.0f;
    }

    printf("[CoalesceFeatureDescriptorsMedian] Loading keys\n");

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;

	/* Load the keys for this image */
	m_image_data[i].LoadKeys(true);
    }
    
    printf("[CoalesceFeatureDescriptors] Coalescing features\n");

    for (int i = 0; i < num_points; i++) {
	int num_views = (int) m_point_data[i].m_views.size();

	if (num_views == 0)
	    continue;

	if (num_views <= 2) {
	    /* Average the keys */
	    for (int j = 0; j < num_views; j++) {
		int view = m_point_data[i].m_views[j].first;
		int key = m_point_data[i].m_views[j].second;
		
		for (int k = 0; k < 128; k++) {
		    m_point_data[i].m_desc[k] += 
                        GetKeyWithDesc(view, key).m_d[k];
		}
	    }

	    for (int k = 0; k < 128; k++) {
		m_point_data[i].m_desc[k] /= num_views;
	    }	    
	} else {
	    /* Find the 'median' key */
	    double min_dist = DBL_MAX;
	    int best_view = -1, best_key = -1;
	    for (int j = 0; j < num_views; j++) {
		double dist = 0.0;
		int v1 = m_point_data[i].m_views[j].first;
		int k1 = m_point_data[i].m_views[j].second;

		for (int k = 0; k < num_views; k++) {
		    if (j == k) continue;

		    int v2 = m_point_data[i].m_views[k].first;
		    int k2 = m_point_data[i].m_views[k].second;
		    
		    for (int l = 0; l < 128; l++) {
			dist += 
			    fabs((double) (GetKeyWithDesc(v1,k1).m_d[l] - 
                                           GetKeyWithDesc(v2,k2).m_d[l]));
		    }
		}

		if (dist < min_dist) {
		    min_dist = dist;
		    best_view = v1;
		    best_key = k1;
		}
	    }

	    /* Save the best view */
	    for (int k = 0; k < 128; k++) {
		m_point_data[i].m_desc[k] = 
		    GetKeyWithDesc(best_view, best_key).m_d[k];
	    }
	}	
    }
}


/* Coalesce feature descriptors for each feature point */
void BundlerApp::CoalesceFeatureDescriptors()
{
    if (m_features_coalesced)
        return;

    int num_points = (int) m_point_data.size();
    int *num_observations = new int[num_points];

    for (int i = 0; i < num_points; i++) {
	num_observations[i] = 0;
	m_point_data[i].m_desc = new float[128];
	for (int j = 0; j < 128; j++) 
	    m_point_data[i].m_desc[j] = 0.0f;
    }

    float *xx = new float[num_points * 128];

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;

	printf("[CoalesceFeatureDescriptors] Adding features from "
	       "image %d\n", i);
        fflush(stdout);

	m_image_data[i].LoadKeys(true);

	for (int j = 0; j < num_points; j++) {
	    int num_views = (int) m_point_data[j].m_views.size();
	    int key_seen = -1;

	    if (num_views == 0)
		continue;

	    for (int k = 0; k < num_views; k++) {
		if (m_point_data[j].m_views[k].first == i) {
		    key_seen = m_point_data[j].m_views[k].second;
		    break;
		}
	    }

	    if (key_seen >= 0) {
                KeypointWithDesc &key = GetKeyWithDesc(i, key_seen);
		for (int k = 0; k < 128; k++) {
                    float x = (float) key.m_d[k];
		    m_point_data[j].m_desc[k] += key.m_d[k];
                    xx[128 * j + k] += x * x;
		}

		num_observations[j]++;
	    }
	}

	m_image_data[i].UnloadKeys();
    }

    double *variance = new double[num_points];

    for (int i = 0; i < num_points; i++) {
	if (num_observations[i] != (int) m_point_data[i].m_views.size())
	    printf("[CoalesceFeatureDescriptors] "
		   "Mismatch in observation count\n");

	if (num_observations[i] == 0) continue;
	
	for (int j = 0; j < 128; j++) {
	    m_point_data[i].m_desc[j] /= num_observations[i];
            xx[128 * i + j] /= num_observations[i];
        }

        variance[i] = 0.0;
        
        for (int j = 0; j < 128; j++) {
            double x = m_point_data[i].m_desc[j];
            variance[i] += xx[128 * i + j] - x * x;
        }
    }

    double median = 
        kth_element(num_points, iround(0.5 * num_points), variance);

    printf("[CoalesceFeatureDescriptors] Median squared variance: "
           "%0.3f\n", median);
    fflush(stdout);

    delete [] num_observations;
    delete [] xx;
    delete [] variance;

    m_features_coalesced = true;
}

#ifndef __DEMO__
/* Compute likely matches between a set of keypoints and the
 * reconstructed points */
std::vector<KeypointMatch> 
   BundlerApp::MatchKeysToPoints(const std::vector<KeypointWithDesc> &k1, 
                                 double ratio) 
{
    ann_1_1_char::annMaxPtsVisit(20000);

    int num_points = (int) m_point_data.size();
    std::vector<KeypointMatch> matches;

    /* Create a new array of points */
    ann_1_1_char::ANNpointArray pts = ann_1_1_char::annAllocPts(num_points, 128);

    for (int i = 0; i < num_points; i++) {
	for (int j = 0; j < 128; j++) {
	    pts[i][j] = m_point_data[i].m_desc[j];
	}
    }
    
    clock_t start = clock();
    /* Create a search tree for k2 */
    ann_1_1_char::ANNkd_tree *tree = new ann_1_1_char::ANNkd_tree(pts, num_points, 128, 4);
    clock_t end = clock();
    
    printf("[MatchKeysToPoints] Building tree took %0.3fs\n", 
	   (end - start) / ((double) CLOCKS_PER_SEC));

    /* Now do the search */
    ann_1_1_char::ANNpoint query = ann_1_1_char::annAllocPt(128);
    start = clock();
    for (int i = 0; i < (int) k1.size(); i++) {
	int j;

	for (j = 0; j < 128; j++) {
	    query[j] = k1[i].m_d[j];
	}

	ann_1_1_char::ANNidx nn_idx[2];
	ann_1_1_char::ANNdist dist[2];

	tree->annkPriSearch(query, 2, nn_idx, dist, 0.0);

	if (sqrt(((double) dist[0]) / ((double) dist[1])) < ratio) {
	    matches.push_back(KeypointMatch(i, nn_idx[0]));
	}
    }
    end = clock();

    printf("[MatchKeysToPoints] Searching tree took %0.3fs\n",
	   (end - start) / ((double) CLOCKS_PER_SEC));

    int num_matches = (int) matches.size();

    printf("[MatchKeys] Found %d matches\n", num_matches);

    /* Cleanup */
    ann_1_1_char::annDeallocPts(pts);
    ann_1_1_char::annDeallocPt(query);

    delete tree;

    return matches;
}

/* Compute likely matches between the reconstructed points and a set
 * of keypoints and the reconstructed points */
std::vector<KeypointMatch> 
   BundlerApp::MatchPointsToKeys(const std::vector<KeypointWithDesc> &keys, 
                                 double ratio) 
{
    ann_1_1_char::annMaxPtsVisit(200);

    int num_points = (int) keys.size();
    std::vector<KeypointMatch> matches;

    /* Create a new array of points */
    ann_1_1_char::ANNpointArray pts = ann_1_1_char::annAllocPts(num_points, 128);

    for (int i = 0; i < num_points; i++)
	for (int j = 0; j < 128; j++)
	    pts[i][j] = keys[i].m_d[j];
    
    clock_t start = clock();
    /* Create a search tree for k2 */
    ann_1_1_char::ANNkd_tree *tree = new ann_1_1_char::ANNkd_tree(pts, num_points, 128, 4);
    clock_t end = clock();
    
    printf("[MatchPointsToKeys] Building tree took %0.3fs\n", 
	   (end - start) / ((double) CLOCKS_PER_SEC));

    /* Now do the search */
    ann_1_1_char::ANNpoint query = ann_1_1_char::annAllocPt(128);
    start = clock();
    for (int i = 0; i < (int) m_point_data.size(); i++) {
	int j;

	for (j = 0; j < 128; j++) {
	    query[j] = m_point_data[i].m_desc[j];
	}

	ann_1_1_char::ANNidx nn_idx[2];
	ann_1_1_char::ANNdist dist[2];

	tree->annkPriSearch(query, 2, nn_idx, dist, 0.0);

	if (sqrt(((double) dist[0]) / ((double) dist[1])) < ratio) {
	    matches.push_back(KeypointMatch(i, nn_idx[0]));
	}
    }
    end = clock();

    printf("[MatchPointsToKeys] Searching tree took %0.3fs\n",
	   (end - start) / ((double) CLOCKS_PER_SEC));

    int num_matches = (int) matches.size();

    printf("[MatchPointsToKeys] Found %d matches\n", num_matches);
    fflush(stdout);
    
    /* Cleanup */
    ann_1_1_char::annDeallocPts(pts);
    ann_1_1_char::annDeallocPt(query);

    delete tree;

    return matches;
}
#endif /* __DEMO__ */

/* Remove matches close to the edges of the two given images */
void BundlerApp::RemoveMatchesNearBorder(int i1, int i2, int border_width)
{
    MatchIndex idx = GetMatchIndex(i1, i2);

    // assert(m_match_lists.find(idx) != m_match_lists.end());
    m_matches.Contains(idx);
    assert(m_image_data[i1].m_keys_loaded);
    assert(m_image_data[i2].m_keys_loaded);
    
    std::vector<KeypointMatch> &list = m_matches.GetMatchList(idx); 
    // m_match_lists[idx];
    int num_matches = (int) list.size();

    int w1 = m_image_data[i1].GetWidth();
    int h1 = m_image_data[i1].GetHeight();
    int w2 = m_image_data[i2].GetWidth();
    int h2 = m_image_data[i2].GetHeight();
    
    double w1_min = -0.5 * w1 + border_width;
    double w1_max = 0.5 * w1 - border_width;
    double h1_min = -0.5 * h1 + border_width;
    double h1_max = 0.5 * h1 - border_width;

    double w2_min = -0.5 * w2 + border_width;
    double w2_max = 0.5 * w2 - border_width;
    double h2_min = -0.5 * h2 + border_width;
    double h2_max = 0.5 * h2 - border_width;

    int num_removed = 0;
    for (int k = 0; k < num_matches; k++) {
        KeypointMatch &m = list[k];
        
        const Keypoint &k1 = m_image_data[i1].m_keys[m.m_idx1];
        const Keypoint &k2 = m_image_data[i2].m_keys[m.m_idx2];

        if (k1.m_x < w1_min || k1.m_x > w1_max || 
            k1.m_y < h1_min || k1.m_y > h1_max ||
            k2.m_x < w2_min || k2.m_x > w2_max || 
            k2.m_y < h2_min || k2.m_y > h2_max) {
            
            /* Erase this match */
            list.erase(list.begin() + k);
            k--;
            num_matches--;

            num_removed++;
        }
    }

    printf("[RemoveMatchesNearBorder] Removed %d matches from pair (%d,%d)\n",
           num_removed, i1, i2);
}

/* Remove matches close to the bottom edge of the two given images */
void BundlerApp::RemoveMatchesNearBottom(int i1, int i2, int border_width)
{
    MatchIndex idx = GetMatchIndex(i1, i2);

    // assert(m_match_lists.find(idx) != m_match_lists.end());
    assert(m_matches.Contains(idx));
    assert(m_image_data[i1].m_keys_loaded);
    assert(m_image_data[i2].m_keys_loaded);
    
    std::vector<KeypointMatch> &list = m_matches.GetMatchList(idx); 
    // m_match_lists[idx];
    int num_matches = (int) list.size();

    int h1 = m_image_data[i1].GetHeight();
    int h2 = m_image_data[i2].GetHeight();
    
    double h1_min = -0.5 * h1 + border_width;
    double h2_min = -0.5 * h2 + border_width;

    int num_removed = 0;
    for (int k = 0; k < num_matches; k++) {
        KeypointMatch &m = list[k];
        
        const Keypoint &k1 = m_image_data[i1].m_keys[m.m_idx1];
        const Keypoint &k2 = m_image_data[i2].m_keys[m.m_idx2];

        if (k1.m_y < h1_min || k2.m_y < h2_min) {
            
            /* Erase this match */
            list.erase(list.begin() + k);
            k--;
            num_matches--;

            num_removed++;
        }
    }

    printf("[RemoveMatchesNearBottom] Removed %d matches from pair (%d,%d)\n",
           num_removed, i1, i2);
}
