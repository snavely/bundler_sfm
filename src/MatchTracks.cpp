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

/* MatchTracks.cpp */
/* Code for processing matches and tracks */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "BaseApp.h"
#include "BundleUtil.h"

int BaseApp::SetTracksFromPoints(int image) 
{
    if (!m_image_data[image].m_keys_loaded) {
        printf("[BaseApp::SetTrackFromPoints] Error: keypoints have not "
               "been loaded\n");

        return 0;
    }

    int num_keys = (int) m_image_data[image].m_keys.size();
    int num_points = (int) m_point_data.size();
    int num_tracks = 0;

    for (int i = 0; i < num_keys; i++) 
        m_image_data[image].m_keys[i].m_track = -1;

    
    for (int i = 0; i < num_points; i++) {
        int num_views = (int) m_point_data[i].m_views.size();
        
        for (int j = 0; j < num_views; j++) {
            ImageKey &k = m_point_data[i].m_views[j];

            if (k.first == image) {
                m_image_data[k.first].m_keys[k.second].m_track = i;
                num_tracks++;
            }
        }
    }

    return num_tracks;
}

void BaseApp::CreateTracksFromPoints()
{
    m_track_data.clear();
    
    int num_images = GetNumImages();
    
    for (int i = 0; i < num_images; i++) {
        m_image_data[i].m_visible_points.clear();
        m_image_data[i].m_visible_keys.clear();
    }

    int num_points = m_point_data.size();
    for (int i = 0; i < num_points; i++) {
        TrackData track;
        track.m_views = m_point_data[i].m_views;

        m_track_data.push_back(track);

        int num_views = (int) m_point_data[i].m_views.size();
        
        for (int j = 0; j < num_views; j++) {
            int v = m_point_data[i].m_views[j].first;
            int k = m_point_data[i].m_views[j].second;
            
            m_image_data[v].m_visible_points.push_back(i);
            m_image_data[v].m_visible_keys.push_back(k);
        }
    }
}

void BaseApp::SetTracksFromPoints()
{
    int num_images = GetNumImages();
    int num_points = (int) m_point_data.size();
    
    for (int i = 0; i < num_images; i++) {
        int num_keys = (int) m_image_data[i].m_keys.size();
        
        for (int j = 0; j < num_keys; j++) {
            m_image_data[i].m_keys[j].m_track = -1;
        }
    }

    for (int i = 0; i < num_points; i++) {
        int num_views = (int) m_point_data[i].m_views.size();
        
        for (int j = 0; j < num_views; j++) {
            ImageKey &k = m_point_data[i].m_views[j];
            
            m_image_data[k.first].m_keys[k.second].m_track = i;
        }
    }
}

void BaseApp::SetTracks(int image) 
{
    printf("[BaseApp::SetTracks] Setting tracks for image %d...\n", image);

    clock_t start = clock();

    // int num_tracks = (int) m_track_data.size();

    // std::vector<TrackData>::iterator t_iter;
    // int i = 0;

    ImageData &img_data = m_image_data[image];
    assert(img_data.m_keys_loaded);

    int num_tracks = (int) img_data.m_visible_points.size();
    
    for (int i = 0; i < num_tracks; i++) {
        int tr = img_data.m_visible_points[i];
        int key = img_data.m_visible_keys[i];

        assert(key < (int) img_data.m_keys.size());

        img_data.m_keys[key].m_track = tr;
    }

    clock_t end = clock();
    
    printf("[BaseApp::SetTracks] Finished in %0.3fs\n", 
	   (double) (end - start) / CLOCKS_PER_SEC);

    fflush(stdout);
}

int BaseApp::GetNumTrackMatches(int img1, int img2) 
{
    const std::vector<int> &tracks1 = m_image_data[img1].m_visible_points;
    const std::vector<int> &tracks2 = m_image_data[img2].m_visible_points;

    // std::vector<int> isect = GetVectorIntersection(tracks1, tracks2);
    // int num_isect = (int) isect.size();

    std::vector<int>::const_iterator iter;
    for (iter = tracks2.begin(); iter != tracks2.end(); iter++) {
        int track_idx = *iter;
        m_track_data[track_idx].m_extra = 0;
    }

    for (iter = tracks1.begin(); iter != tracks1.end(); iter++) {
        int track_idx = *iter;
        m_track_data[track_idx].m_extra = 1;
    }

    int num_isect = 0;
    for (iter = tracks2.begin(); iter != tracks2.end(); iter++) {
        int track_idx = *iter;
        num_isect += m_track_data[track_idx].m_extra;
    }

    return num_isect;
}

void BaseApp::SetMatchesFromTracks(int img1, int img2)
{
    std::vector<int> &tracks1 = m_image_data[img1].m_visible_points;
    std::vector<int> &tracks2 = m_image_data[img2].m_visible_points;

    std::vector<int> isect = GetVectorIntersection(tracks1, tracks2);
    
    int num_isect = (int) isect.size();

    if (num_isect == 0)
        return;
    
    MatchIndex idx = GetMatchIndex(img1, img2);

    std::vector<KeypointMatch> &matches = m_matches.GetMatchList(idx); 
    // m_match_lists[idx];

    matches.clear();
    matches.resize(num_isect);

    for (int i = 0; i < num_isect; i++) {
        int tr = isect[i];

#if 0
        int num_views = (int) m_track_data[tr].m_views.size();
        int k1 = -1, k2 = -1;

        for (int j = 0; j < num_views; j++) {
            if (m_track_data[tr].m_views[j].first == img1) {
                k1 = m_track_data[tr].m_views[j].second;
            } 

            if (m_track_data[tr].m_views[j].first == img2) {
                k2 = m_track_data[tr].m_views[j].second;
            } 
        }

        assert(k1 != -1 && k2 != -1);
#endif

        std::pair<std::vector<int>::const_iterator, 
            std::vector<int>::const_iterator> p;
        const std::vector<int> &pt1 = m_image_data[img1].m_visible_points;
        p = equal_range(pt1.begin(), pt1.end(), tr);
        assert(p.first != p.second);
        int offset = p.first - pt1.begin();
        int k1 = m_image_data[img1].m_visible_keys[offset];

        const std::vector<int> &pt2 = m_image_data[img2].m_visible_points;
        p = equal_range(pt2.begin(), pt2.end(), tr);
        assert(p.first != p.second);
        offset = p.first - pt2.begin();
        int k2 = m_image_data[img2].m_visible_keys[offset];

        matches[i] = KeypointMatch(k1, k2);
    }
}

void BaseApp::SetMatchesFromTracks() 
{
    /* Clear all matches */
    // ClearMatches();
    RemoveAllMatches();

    int num_tracks = (int) m_track_data.size();

    int num_tracks_used = 0;
    for (int i = 0; i < num_tracks; i++) {
	TrackData &t = m_track_data[i];
	
	int num_views = (int) t.m_views.size();
	
	if (num_views < m_min_track_views) 
	    continue; /* Not enough observations */

	if (num_views > m_max_track_views) 
	    continue; /* Too many observations */

	for (int j = 0; j < num_views; j++) {
	    for (int k = 0; k < num_views; k++) {
		if (j == k) continue;
		
		int v1 = t.m_views[j].first;
		int v2 = t.m_views[k].first;
		
		int k1 = t.m_views[j].second;
		int k2 = t.m_views[k].second;
		
		MatchIndex idx = GetMatchIndex(v1, v2);

		// m_matches[idx] = true;
                SetMatch(v1, v2);
		// m_match_lists[idx].push_back(KeypointMatch(k1,k2));

                m_matches.GetMatchList(idx).push_back(KeypointMatch(k1, k2));
	    }
	}

        num_tracks_used++;
    }

    printf("[BaseApp::SetMatchesFromTracks] Used %d tracks\n", 
           num_tracks_used);
}

/* Use the bundle-adjusted points to create a new set of matches */
void BaseApp::SetMatchesFromPoints(int threshold)
{
    printf("[BaseApp::SetMatchesFromPoints] Setting up matches...\n");

    /* Clear all matches */
    // ClearMatches();
    RemoveAllMatches();

    int num_points = (int) m_point_data.size();
    for (int i = 0; i < num_points; i++) {
	int num_views = (int) m_point_data[i].m_views.size();
        if (num_views < threshold)
            continue;

	for (int j = 0; j < num_views; j++) {
	    for (int k = 0; k < num_views; k++) {
		if (j == k) continue;
		
		ImageKey view1 = m_point_data[i].m_views[j];
		ImageKey view2 = m_point_data[i].m_views[k];

		int v1 = view1.first;
		int v2 = view2.first;

		int k1 = view1.second;
		int k2 = view2.second;

		KeypointMatch m;
		
		m.m_idx1 = k1;
		m.m_idx2 = k2;

		// m_matches[v1 * num_images + v2] = true;
                SetMatch(v1, v2);
                MatchIndex idx = GetMatchIndex(v1, v2);
		// m_match_lists[idx].push_back(m);
                m_matches.AddMatch(idx, m);
	    }
	}
    }

    printf("[BaseApp::SetMatchesFromPoints] Done!\n");
}

#if 0
void BaseApp::ClearMatches(MatchIndex idx)
{
    if (m_match_lists.find(idx) != m_match_lists.end()) {
        m_match_lists[idx].clear();
        m_match_lists.erase(m_match_lists.find(idx));
    }
}
#endif

/* Make match lists symmetric */
void BaseApp::MakeMatchListsSymmetric() 
{
    unsigned int num_images = GetNumImages();

    std::vector<MatchIndex> matches;

    for (unsigned int i = 0; i < num_images; i++) {
        MatchAdjList::const_iterator iter;
        for (iter = m_matches.Begin(i); iter != m_matches.End(i); iter++) {
            // unsigned int i = iter->first;
            unsigned int j = iter->m_index; // iter->second;

            if (j <= i)
                continue;

            assert(ImagesMatch(i, j));   

            // MatchIndex idx = *iter; 
            MatchIndex idx = GetMatchIndex(i, j);
            MatchIndex idx_rev = GetMatchIndex(j, i);
            // int num_matches = (int) m_match_lists[idx].size();

            const std::vector<KeypointMatch> &list = iter->m_match_list;
            unsigned int num_matches = list.size();

            // m_match_lists[idx_rev].clear();
            m_matches.SetMatch(idx_rev);
            m_matches.ClearMatch(idx_rev);

            for (unsigned int k = 0; k < num_matches; k++) {
                KeypointMatch m1, m2;
		
                m1 = list[k];
                
                m2.m_idx1 = m1.m_idx2;
                m2.m_idx2 = m1.m_idx1;

                // m_match_lists[idx_rev].push_back(m2);
                m_matches.AddMatch(idx_rev, m2);
            }

            matches.push_back(idx);
        }
    }
    
    unsigned int num_matches = matches.size();

    for (unsigned int i = 0; i < num_matches; i++) {
        unsigned int img1 = matches[i].first;
        unsigned int img2 = matches[i].second;
        SetMatch(img2, img1);
    }

    matches.clear();
}

/* Prune points that match to multiple targets */
void BaseApp::PruneDoubleMatches() 
{
    unsigned int num_images = GetNumImages();

    for (unsigned int i = 0; i < num_images; i++) {
        MatchAdjList::iterator iter;

        std::vector<unsigned int> remove;
        for (iter = m_matches.Begin(i); iter != m_matches.End(i); iter++) {
            HashSetInt seen;

            int num_pruned = 0;
            // MatchIndex idx = *iter; // GetMatchIndex(i, j);

            std::vector<KeypointMatch> &list = iter->m_match_list;

            /* Unmark keys */
            // int num_matches = (int) m_match_lists[idx].size();
            int num_matches = (int) list.size();

            for (int k = 0; k < num_matches; k++) {
                int idx2 = list[k].m_idx2;
		
                // if (GetKey(j,idx2).m_extra != -1) {
                if (seen.find(idx2) != seen.end()) {
                    /* This is a repeat */
                    // printf("[%d] Pruning repeat %d\n", i, idx2);
                    list.erase(list.begin() + k);
                    num_matches--;
                    k--;
                    
                    num_pruned++;
                } else {
                    /* Mark this key as matched */
                    // GetKey(j,idx2).m_extra = k;
                    seen.insert(idx2);
                }
            }

            // unsigned int i = iter->first;
            // unsigned int j = iter->second;
            unsigned int j = iter->m_index; // first;

            printf("[PruneDoubleMatches] Pruned[%d,%d] = %d / %d\n",
                   i, j, num_pruned, num_matches + num_pruned);

            if (num_matches < m_min_num_feat_matches) {
                /* Get rid of... */
                remove.push_back(iter->m_index); // first);
            }
        }

        for (unsigned int j = 0; j < remove.size(); j++) {
            int idx2 = remove[j];
            m_matches.RemoveMatch(GetMatchIndex(i, idx2));
            printf("[PruneDoubleMatches] Removing[%d,%d]\n", i, idx2);
        }
    }
}
