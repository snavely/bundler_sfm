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

/* BaseApp.h */
/* Base application */

#ifndef __baseapp_h__
#define __baseapp_h__


#ifndef __USE_ANN__
#include "BruteForceSearch.h"
#else
#include "anniface.h"
#endif /* __USE_ANN__ */

#include "ImageData.h"


#ifndef __DEMO__
#include "sfm.h"
#endif /* __DEMO__ */

#include "defines.h"

#include <assert.h>
#include <algorithm>
#include <list>

#ifndef WIN32
#include <unordered_map>
#include <unordered_set>
#else
#include <hash_map>
#include <hash_set>
#endif

// #ifndef __DEMO__
#ifdef __USE_BOOST__
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
using namespace boost;
#endif /* __USE_BOOST__ */
// #endif /* __DEMO__ */

#ifdef __USE_BOOST__
typedef adjacency_list<vecS, vecS, 
                       undirectedS, 
                       property<vertex_color_t, int>, 
                       property<edge_weight_t, int> > ImageGraph;
#endif

class TransformInfo {
public:

    /* File IO routines */
    void ReadFromFile(FILE *f);
    void WriteToFile(FILE *f);

    /* For object movies */
    double m_fmatrix[9];
    double m_ematrix[9];

    /* For homographies */
    double m_H[9];
    double m_inlier_ratio;
    int m_num_inliers;

    /* For color correction */
    double m_gain[3], m_bias[3];
};

typedef std::pair<unsigned long, unsigned long> MatchIndex;
// typedef unsigned long long MatchIndex;

#ifdef WIN32
namespace stdext {
    template<>
    class hash_compare<MatchIndex> {
	public:
		static const size_t bucket_size = 4;
        static const size_t min_buckets = 8;
        size_t
        operator()(const MatchIndex &__x) const
        { return __x.first * 1529 + __x.second; }

        bool operator()(const MatchIndex &__x1, const MatchIndex &__x2) const {
			return (__x1.first < __x2.first) || (__x1.first == __x2.first && __x1.second < __x2.second);
        }
    };
}
#else
namespace std {
    template<>
    struct hash<MatchIndex> {
        size_t
        operator()(MatchIndex __x) const
        { return __x.first * 1529 + __x.second; }
    };
}
#endif

/* Table containing information about which pairs of images match */
#if 0
#ifndef WIN32
class MatchTable : private  __gnu_cxx::hash_set<MatchIndex>
#else
class MatchTable : private  stdext::hash_set<MatchIndex>
#endif
{
public:
    typedef __gnu_cxx::hash_set<MatchIndex>::const_iterator const_iterator;

    void SetMatch(MatchIndex idx) { 
        insert(idx);
    }

    void RemoveMatch(MatchIndex idx) {
        erase(idx);
    }

    bool Contains(MatchIndex idx) {
        return (find(idx) != end());
    }

    void Clear() {
        clear();
    }

    iterator Begin() {
        return begin();
    }
    
    iterator End() {
        return end();
    }
};
#endif

#ifdef WIN32
typedef stdext::hash_map<unsigned int, std::vector<KeypointMatch> >
   MatchAdjTable;
#else
typedef unordered_map<unsigned int, std::vector<KeypointMatch> >
   MatchAdjTable;
#endif

class AdjListElem {
public:
    bool operator< (const AdjListElem &other) const {
        return m_index < other.m_index;
    }
    
    unsigned int m_index;
    std::vector<KeypointMatch> m_match_list;
};

typedef std::vector<AdjListElem> MatchAdjList;

class MatchTable
{
    // typedef __gnu_cxx::hash_set<MatchIndex>::const_iterator const_iterator;
public:

    MatchTable() { }

    MatchTable(int num_images) {
        m_match_lists.resize(num_images);
        // m_neighbors.resize(num_images);
    }

    void SetMatch(MatchIndex idx) { 
        if (Contains(idx))
            return;  // already set

        /* Create a new list */
        // m_match_lists[idx.first][idx.second] = std::vector<KeypointMatch> ();
        // m_match_lists[idx.first].insert(idx.second);
        // std::list<unsigned int> tmp;
        // tmp.push_back(idx.second);
        // m_neighbors[idx.first].merge(tmp);
#if 0
        MatchAdjList tmp;
        AdjListElem adjlist_elem;
        adjlist_elem.m_index = idx.second;
        tmp.push_back(adjlist_elem);
        m_match_lists[idx.first].merge(tmp);
#else
        /* Using vector */
        AdjListElem e;
        e.m_index = idx.second;
        MatchAdjList &l = m_match_lists[idx.first];
        MatchAdjList::iterator p = lower_bound(l.begin(), l.end(), e);
        l.insert(p, e);
#endif
    }

    void AddMatch(MatchIndex idx, KeypointMatch m) {
        assert(Contains(idx));
        // m_match_lists[idx.first][idx.second].push_back(m);
        GetMatchList(idx).push_back(m);
    }

    void ClearMatch(MatchIndex idx) { // but don't erase!
        if (Contains(idx)) {
            // m_match_lists[idx.first][idx.second].clear();
            GetMatchList(idx).clear();
        }
    }
    
    void RemoveMatch(MatchIndex idx) {
        if (Contains(idx)) {
            // m_match_lists[idx.first][idx.second].clear();
            // m_match_lists[idx.first].erase(idx.second);
            std::vector<KeypointMatch> &match_list = GetMatchList(idx);
            match_list.clear();

            // Remove the neighbor
#if 0
            std::list<unsigned int> &l = m_neighbors[idx.first];
            std::pair<std::list<unsigned int>::iterator,
                      std::list<unsigned int>::iterator> p = 
                equal_range(l.begin(), l.end(), idx.second);

            assert(p.first != l.end());
            
            l.erase(p.first, p.second);
#endif
            AdjListElem e;
            e.m_index = idx.second;
            MatchAdjList &l = m_match_lists[idx.first];
            std::pair<MatchAdjList::iterator, MatchAdjList::iterator> p = 
                equal_range(l.begin(), l.end(), e);

            assert(p.first != p.second); // l.end());
            
            l.erase(p.first, p.second);        
        }
    }

    unsigned int GetNumMatches(MatchIndex idx) {
        if (!Contains(idx))
            return 0;
        
        // return m_match_lists[idx.first][idx.second].size();
        return GetMatchList(idx).size();
    }

    std::vector<KeypointMatch> &GetMatchList(MatchIndex idx) {
        // assert(Contains(idx));
        // return m_match_lists[idx.first][idx.second];

        AdjListElem e;
        e.m_index = idx.second;
        MatchAdjList &l = m_match_lists[idx.first];
        std::pair<MatchAdjList::iterator, MatchAdjList::iterator> p = 
            equal_range(l.begin(), l.end(), e);
    
        assert(p.first != p.second); // l.end());
        
        return (p.first)->m_match_list;
    }
    
    bool Contains(MatchIndex idx) const {
        // return (m_match_lists[idx.first].find(idx.second) != 
        //         m_match_lists[idx.first].end());
        AdjListElem e;
        e.m_index = idx.second;
        const MatchAdjList &l = m_match_lists[idx.first];
        std::pair<MatchAdjList::const_iterator, 
            MatchAdjList::const_iterator> p = 
            equal_range(l.begin(), l.end(), e);
    
        return (p.first != p.second); // l.end());
    }

    void RemoveAll() {
        int num_lists = m_match_lists.size();

        for (int i = 0; i < num_lists; i++) {
            m_match_lists[i].clear();
            // m_neighbors[i].clear();
        }
    }

    unsigned int GetNumNeighbors(unsigned int i) {
        return m_match_lists[i].size();
    }

#if 0
    std::list<unsigned int> GetNeighbors(unsigned int i) {
        // return m_neighbors[i];
        std::list<unsigned int> nbrs;
        MatchAdjList::iterator p;
        for (p = Begin(i); p != End(i); p++) {
            nbrs.push_back(p->m_index);
        }
        return nbrs;
    }
#endif

    MatchAdjList &GetNeighbors(unsigned int i) {
        return m_match_lists[i];
    }

    MatchAdjList::iterator Begin(unsigned int i) {
        return m_match_lists[i].begin();
    }
    
    MatchAdjList::iterator End(unsigned int i) {
        return m_match_lists[i].end();
    }
    
private:
    // std::vector<MatchAdjTable> m_match_lists;
    // std::vector<KeypointMatchList> m_match_lists;
    // std::vector<std::list<unsigned int> > m_neighbors;
    std::vector<MatchAdjList> m_match_lists;
};

/* Return the match index of a pair of images */
MatchIndex GetMatchIndex(int i1, int i2);
MatchIndex GetMatchIndexUnordered(int i1, int i2);
// #include "GetMatchIndex.h"

class BaseApp 
{
public:
    virtual ~BaseApp() { }

    virtual bool OnInit() = 0;

    /* Process command line options */
    virtual void ProcessOptions(int argc, char **argv) = 0;

    /* Return the number of images */
    int GetNumImages();
    int GetNumOriginalImages();
    int GetNumMatches(int i1, int i2);

    int FindImageWithName(const char *name);

    /* Get match information */
    void SetMatch(int i1, int i2);
    void RemoveMatch(int i1, int i2);
    bool ImagesMatch(int i1, int i2);

    /* Get keys */
    Keypoint &GetKey(int img, int key);
    KeypointWithDesc &GetKeyWithDesc(int img, int key);
    int GetNumKeys(int img);
    /* Get the index of a registered camera */
    int GetRegisteredCameraIndex(int cam);

    /* Load a list of image names from a file */
    void LoadImageNamesFromFile(FILE *f);

    /* Load matches from files */
    void LoadMatches();
    void ReadMatchFile(int i, int j);
    void LoadMatchTable(const char *filename);
    void LoadMatchIndexes(const char *index_dir);
    /* Load keys from files */
    void LoadKeys(bool descriptor = true);
    void RemoveAllMatches();
    /* Prune points that match to multiple targets */
    void PruneDoubleMatches();

    /* Make match lists symmetric */
    void MakeMatchListsSymmetric();

    /* Setup data structures storing the points and lines visible to
     * each image */
    void SetupImagePoints(int min_views = 1);

    /* Output routines */

    /* IO routines */
    void ReadGeometricConstraints(const char *filename);
    void WriteGeometricConstraints(const char *filename);
    void WriteTracks(char *filename);
    void WriteTracks2(char *filename);
    void ReadCameraConstraints();
    void ReadPointConstraints();
    void ReadIntrinsicsFile();
    bool ReadTrackPairs(const char *filename);
    void WriteTrackPairs(const char *filename);

    /* Read the ignore file */
    void ReadIgnoreFile();

    /* Grab the color of each keypoint */
    void ReadKeyColors();

    /* Read in information about the world */
    void ReadBundleFile(const char *filename);
    void ReloadBundleFile (const char *filename);

    /* Read/write line segments */
    void ReadLines3D(const char *filename);
    void WriteLines3D(const char *filename);

    /* Clear the current model */
    void ClearModel();

    /* Read / write the match table */
    void ReadMatchTable(const char *append = "");
    void WriteMatchTable(const char *append = "");

    /* Initialize images read from a file without performing bundle
     * adjustment */
    void InitializeImagesFromFile(FILE *f);

    void SetTracks(int image);
    /* Set the track field of each keypoint given the 3D points */
    void CreateTracksFromPoints();
    void SetTracksFromPoints();
    int SetTracksFromPoints(int image);
    /* Use tracks to setup matches */
    void SetMatchesFromTracks();
    void SetMatchesFromTracks(int img1, int img2);
    // void ClearMatches(MatchIndex idx);
    int GetNumTrackMatches(int img1, int img2);

    /* Use the bundle-adjusted points to create a new set of matches */
    void SetMatchesFromPoints(int threshold = 0);

    void ReindexPoints();

    /* Create a search tree on cameras */
#ifdef __USE_ANN__
    ANNkd_tree *CreateCameraSearchTree();
#else
    BruteForceSearch *CreateCameraSearchTree();
#endif

#ifndef __DEMO__
    /* Write point files to a ply file */
    void DumpPointsToPly(const char *output_directory, const char *filename, 
                         int num_points, int num_cameras, 
			 v3_t *points, v3_t *colors, camera_params_t *cameras
                         /*bool reflect = true*/);

    /* Dump an output file containing information about the current
     * state of the world */
    void DumpOutputFile(const char *output_dir, const char *filename, 
			int num_images, int num_cameras, int num_points,
			int *added_order, 
			camera_params_t *cameras, v3_t *points, v3_t *colors,
			std::vector<ImageKeyVector> &pt_views
                        /*bool output_radial_distortion = false*/);

#endif /* __DEMO__ */

    /* XML output routines */
    void WritePointsXML(const char *filename);
    void WritePointsGeoXML(const char *filename);
    void WriteCamerasXML(const char *filename);
    void WriteCamerasGeoXML(const char *filename);

    /* Write POINTS in a FILE */
    void WritePoints(const char* filename);
    
    /* Write camera parameters in a FILE */
    void WriteCameras(const char* filename);

    /* Fix cameras (reflection bug) */
    void FixReflectionBug();

    void RemoveBadImages(int min_num_points);

    /* Estimate point normals and confidences */
    void EstimatePointNormalsConfidence();
    void EstimatePointNormals();

    /* Estimate the axes based on the image orientations */
    void EstimateAxes(double *xaxis, double *yaxis, double *zaxis);
    /* Analyze the scene and return various attributes */
    void SetupScene(double *center, double *up, 
		    double *x_axis, double *z_axis, double &scale);

    void SetupSceneGroundPlane(double *center, double *up, 
                               double *x_axis, double *z_axis, double &scale);

    void TransformWorldReal();
    virtual void RepositionScene(double *center_out, double *R_out, 
                                 double &scale_out);
    void TransformSceneCanonical(int c1, int c2);
    void UnscaleCameras(int start_camera);

    /* Compute orientations for each image */
    void ComputeImageRotations();

    /* Returns true if i1 and i2 are part of a panorama */
    bool ImagesPartOfPanorama(int i1, int i2);

    /* **** Data **** */

    double m_bundle_version;

    /* Geometry data */

    std::vector<ImageData> m_image_data;   /* Image data */
    int m_num_original_images;

    std::vector<PointData> m_point_data;   /* Information about 3D
					    * points in the scene */

    std::vector<int> m_num_views_orig;

    std::vector<TrackData> m_track_data;   /* Information about the
                                            * detected 3D tracks */

    ImageKeyVector m_outliers;             /* Outliers detected among
					    * the feature points */    



    PlaneData m_ground_plane;        /* Ground plane */

    double m_repos_R[9];
    double m_repos_d[3];
    double m_repos_scale;

    double m_xform[16];       /* Scene transform */

    /* Point constraint fields */
    bool m_use_point_constraints;
    v3_t *m_point_constraints;
    double m_point_constraint_weight;
    char *m_point_constraint_file;

    bool m_matches_loaded;    /* Have the matches been loaded? */
    bool m_matches_computed;  /* Have the matches been computed? */
    // bool *m_matches;          /* Match matrix */    

    MatchTable m_matches;

#if 0
#ifndef WIN32
    __gnu_cxx::hash_map<MatchIndex, std::vector<KeypointMatch> > m_match_lists;
#else
    stdext::hash_map<MatchIndex, std::vector<KeypointMatch> > m_match_lists;
#endif
#endif

#ifndef WIN32
    unordered_map<MatchIndex, TransformInfo> m_transforms;
#else
    stdext::hash_map<MatchIndex, TransformInfo> m_transforms;
#endif

    /* **** Options **** */

    double m_scale;              /* Scene scale */
    bool m_metric;               /* Is the scene in meters? */

    bool m_fisheye;              /* Use fisheye images? */
    char *m_fisheye_params;      /* Fisheye parameter file */

    char *m_ignore_file;         /* Contains the images to ignore during
                                  * bundle adjustment */
    bool m_use_intrinsics;
    char *m_intrinsics_file;     /* File containing intrinsics */     

    bool m_bundle_provided;      /* Was a bundle adjustment file given? */
    char *m_bundle_file;         /* Bundle file */

    const char *m_match_directory;     /* Which directory are the matches
                                        * stored in? */
    const char *m_match_index_dir;     /* Which directory are match indexes
                                        * stored in? */
    const char *m_match_table;         /* File where match table is stored */
    const char *m_key_directory;
    const char *m_image_directory;
    const char *m_sift_binary;         /* Where can we find the sift binary? */

    bool m_estimate_up_vector_szeliski;  /* Estimate the up vector
					  * using Rick's method? */

    int m_min_track_views;           /* Minimum number of observations
				      * for each track */
    int m_max_track_views;           /* Minimum number of observations
				      * for each track */

    int m_min_num_feat_matches;      /* Number of features matches for
                                      * an image pair to be considered
                                      * a match */

    int m_up_image;                  /* What image should we use to
				      * estimate the up vector (if
				      * any?) */

    int m_start_camera;              /* What image should we start
				      * viewing first? */

#ifdef __USE_BOOST__
    /* Graph fields */    
    ImageGraph m_image_graph;
    ImageGraph m_working_image_graph;
    std::vector<int> m_working_images;

    std::vector<int> m_graph_components;
    int m_max_graph_component;
#endif /* __DEMO__ */
};

#endif /* __baseapp_h__ */
