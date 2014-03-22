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

/* BundlerApp.h */
/* Bundler application */

#ifndef __bundlerapp_h__
#define __bundlerapp_h__

#include "BaseApp.h"
#include "LinkDirection.h"
#include "TwoFrameModel.h"

typedef std::pair<int,int> ImagePair;

class BundlerApp : public BaseApp
{
public:
    BundlerApp() {
        /* Set initial values */
        m_bundle_version = 0.3;

        m_fisheye = false;
        m_fixed_focal_length = true;
        m_estimate_distortion = false;
        m_construct_max_connectivity = false;
        m_bundle_provided = false;
        m_analyze_matches = false;

        m_estimate_ignored = false;

        m_use_constraints = false;
        m_constrain_focal = false;
        m_constrain_focal_weight = 100.0;
        m_distortion_weight = 1.0e2;
        
        m_use_point_constraints = false;
        m_point_constraint_weight = 0.0;
        m_point_constraints = NULL;
        m_point_constraint_file = NULL;

        m_only_bundle_init_focal = false;
        m_init_focal_length = 532.0;
        m_initial_pair[0] = -1;
        m_initial_pair[1] = -1;

        m_panorama_mode = false;
        m_homography_threshold = 6.0;
        m_homography_rounds = 256;
        m_fmatrix_threshold = 9.0;
        m_fmatrix_rounds = 2048;
        m_skip_fmatrix = false;
        m_skip_homographies = false;
        m_projection_estimation_threshold = 4.0; // 1.8;
        m_min_proj_error_threshold = 8.0;
        m_max_proj_error_threshold = 16.0;
        m_min_camera_distance_ratio = 0.0;
        m_baseline_threshold = -1.0;
        m_optimize_for_fisheye = false;
        m_use_focal_estimate = false;
        m_trust_focal_estimate = false;
        m_factor_essential = true;
        m_up_image = -1;
        // m_start_camera = -1;
        m_min_track_views = 2;
        m_max_track_views = 100000;
        m_min_num_feat_matches = 16;
        m_min_max_matches = 16;
        m_num_matches_add_camera = -1; /* No maximum by default */
        m_ray_angle_threshold = 2.0;

        m_keypoint_border_width = 0;
        m_keypoint_border_bottom = 0;

        m_fisheye_params = NULL;
        m_bundle_output_file = m_bundle_output_base = NULL;
        m_bundle_file = NULL;
        m_intrinsics_file = NULL;
        m_match_directory = ".";
        m_match_index_dir = NULL;
        m_match_table = NULL;
        m_key_directory = ".";
        m_image_directory = ".";
        m_output_directory = ".";
        m_use_intrinsics = false;
        
        m_matches_computed = false;
        m_ann_max_pts_visit = 400;
        
        m_matches_loaded = false;
        m_features_coalesced = false;

        m_assemble = false;
        m_run_bundle = false;
        m_rerun_bundle = false;
        m_fast_bundle = true;
#ifdef __USE_CERES__
        m_use_ceres = false;
#endif /* __USE_CERES__ */
        m_skip_full_bundle = false;
        m_skip_add_points = false;
        m_use_angular_score = false;

        m_compress_list = false;
        m_reposition_scene = false;
        m_prune_bad_points = false;
        m_predict_next_image = false;
        m_prediction_image = NULL;
        m_scale_focal = 1.0;
        m_scale_focal_file = NULL;
        m_rotate_cameras_file = NULL;
        m_output_relposes = false;
        m_output_relposes_file = NULL;

        m_compute_covariance = false;
        m_covariance_fix1 = -1;
        m_covariance_fix2 = -1;

        m_track_file = NULL;
        m_zero_distortion_params = false;
        m_enrich_points = false;
        m_fix_necker = false;

        

        m_ignore_file = NULL;
        m_add_image_file = NULL;
        m_add_images_fast = false;

        m_scale = 1.0;
        // matrix_ident(3, m_repos_R);
        m_repos_R[0] = 1.0; m_repos_R[1] = 0.0; m_repos_R[2] = 0.0;
        m_repos_R[3] = 0.0; m_repos_R[4] = 1.0; m_repos_R[5] = 0.0;
        m_repos_R[6] = 0.0; m_repos_R[7] = 0.0; m_repos_R[8] = 1.0;

        m_repos_d[0] = m_repos_d[1] = m_repos_d[2] = 0.0;
        m_repos_scale = 1.0;
    
        m_metric = false;

        m_estimate_up_vector_szeliski = false;
    
        // bool load_file = false;
    }

    virtual bool OnInit();

    /* Process command line options */
    virtual void ProcessOptions(int argc, char **argv);

    /* Create a search tree for all keypoints */


    /* Return the number of parameters used to model each camera */
    int GetNumCameraParameters();

    /* Enrich the set of correspondences */
    void EnrichCorrespondences(double alpha, double threshold);

    /* Prune image matches that are not well-supported */
    void PruneMatchesThreshold(int threshold);
    /* Remove matches close to the edges of the two given images */
    void RemoveMatchesNearBorder(int i1, int i2, int border_width);
    /* Remove matches close to the bottom edge of the two given images */
    void RemoveMatchesNearBottom(int i1, int i2, int border_width);

    /* Compute a transform between a given pair of images */
    bool ComputeTransform(int idx1, int idx2, bool removeBadMatches);

    /* Compute transforms between all matching images */
    void ComputeTransforms(bool removeBadMatches, int new_image_start = 0);

    /* Compute epipolar geometry between a given pair of images */
    bool ComputeEpipolarGeometry(int idx1, int idx2, bool removeBadMatches);

    /* Compute epipolar geometry between all matching images */
    void ComputeEpipolarGeometry(bool removeBadMatches, 
				 int new_image_start = 0);

    /* Compute a set of tracks that explain the matches */
    void ComputeTracks(int new_image_start = 0);

    /* Compute geometric information about image pairs */
    void ComputeGeometricConstraints(bool overwrite = false, 
				     int new_image_start = 0);

#ifndef __DEMO__
    /* Set constraints on cameras */
    void SetCameraConstraints(int cam_idx, camera_params_t *params);
    void SetFocalConstraint(const ImageData &data, camera_params_t *params);
    void ClearCameraConstraints(camera_params_t *params);
#endif /* __DEMO__ */

    void CheckPointKeyConsistency(const std::vector<ImageKeyVector> pt_views,
                                  int *added_order);

#ifndef __DEMO__
    /* Initialize the bundle adjustment procedure (loading an existing
     * model if one exists) */
    void InitializeBundleAdjust(int &num_init_cams,
				int *added_order,
				int *added_order_inv,
				camera_params_t *cameras,
				v3_t *points, v3_t *colors,
				std::vector<ImageKeyVector> &pt_views,
				bool use_constraints);

    /* Set up the matrix of projections and the visibility mask */
    void SetupProjections(int num_cameras, int num_points, 
			  int *added_order,
			  v2_t *projections, char *vmask);

    /* Find the camera with the most matches to existing points */
    int FindCameraWithMostMatches(int num_cameras, int num_points,
				  int *added_order,
				  int &parent_idx, int &max_matches,
				  const std::vector<ImageKeyVector> &pt_views);

    /* Find all cameras with at least N matches to existing points */
    std::vector<ImagePair> FindCamerasWithNMatches(int n, 
						   int num_cameras, 
						   int num_points,
						   int *added_order,
						   const std::vector<ImageKeyVector> &pt_views);

    /* Find the camera that would allow use to "grow" the scene as
     * much as possible */
    int FindCameraWithMostConnectivity(int num_cameras, int num_points,
				       int *added_order,
				       int &parent_idx, 
				       int &max_matches);

    /* Triangulate a subtrack */
    v3_t TriangulateNViews(const ImageKeyVector &views, 
			   int *added_order, camera_params_t *cameras,
			   double &error, bool explicit_camera_centers);

    v3_t GeneratePointAtInfinity(const ImageKeyVector &views, 
                                 int *added_order, 
                                 camera_params_t *cameras,
                                 double &error, 
                                 bool explicit_camera_centers);
    
    /* Add new points to the bundle adjustment */
    int BundleAdjustAddNewPoints(int camera_idx, 
				 int num_points, int num_cameras,
				 int *added_order,
				 camera_params_t *cameras,
				 v3_t *points, v3_t *colors,
				 double reference_baseline,
				 std::vector<ImageKeyVector> &pt_views);

    /* Add new points to the bundle adjustment */
    int BundleAdjustAddAllNewPoints(int num_points, int num_cameras,
				    int *added_order,
				    camera_params_t *cameras,
				    v3_t *points, v3_t *colors,
				    double reference_baseline,
				    std::vector<ImageKeyVector> &pt_views,
				    double max_reprojection_error = 16.0,
                                    int min_views = 2);
    
    /* Remove bad points and cameras from a reconstruction */
    int RemoveBadPointsAndCameras(int num_points, int num_cameras, 
                                  int *added_order, 
                                  camera_params_t *cameras,
                                  v3_t *points, v3_t *colors,
                                  std::vector<ImageKeyVector> &pt_views);

    /* Compute pose of all cameras */
    void BundleAdjust();

    /* Quickly compute pose of all cameras */
    void BundleAdjustFast();

    
    /* Estimate poses of all ignored cameras */
    void EstimateIgnoredCameras(int &curr_num_cameras,
                                camera_params_t *cameras,
                                int *added_order,
                                int &pt_count,
                                v3_t *points, 
                                v3_t *colors,
                                std::vector<ImageKeyVector> &pt_views);

    /* Pick a good initial pair of cameras to bootstrap the bundle
     * adjustment */
    void BundlePickInitialPair(int &i_best, int &j_best, 
                               bool use_init_focal_only);

    /* Setup the initial camera pair for bundle adjustment */
    int SetupInitialCameraPair(int i_best, int j_best,
			       double &init_focal_length_0,
			       double &init_focal_length_1,
			       camera_params_t *cameras,
			       v3_t *points, v3_t *colors,
			       std::vector<ImageKeyVector> &pt_views);

    /* Initialize a single image */
    void BundleImage(char *filename, int parent_img);
    /* Initialize images read from a file */
    void BundleImagesFromFile(FILE *f);


    /* Initialize an image for bundle adjustment */
    camera_params_t 
        BundleInitializeImage(ImageData &data, 
                              int image_idx, int camera_idx,
                              int num_cameras, int num_points,
                              int *added_order, v3_t *points,
                              camera_params_t *parent,
                              camera_params_t *cameras, 
                              std::vector<ImageKeyVector> &pt_views,
                              bool *success_out = NULL,
			      bool refine_cameras_and_points = false);

    /* Initialize an image for bundle adjustment (running a full
     * optimization) */
    void BundleInitializeImageFullBundle(int image_idx, int parent_idx,
					 int num_cameras,
					 int num_points, 
					 int *added_order,
					 camera_params_t *cameras,
					 v3_t *points, v3_t *colors,
					 std::vector<ImageKeyVector> 
					     &pt_views);


    /* Refine a set of 3D points */
    double RefinePoints(int num_points, v3_t *points, v2_t *projs,
			int *pt_idxs, camera_params_t *cameras,
			int *added_order,
			const std::vector<ImageKeyVector> &pt_views,
			camera_params_t *camera_out);

    /* Refine a given camera and the points it observes */
    std::vector<int> RefineCameraAndPoints(const ImageData &data, 
                                           int num_points,
					   v3_t *points, v2_t *projs,
					   int *pt_idxs, 
					   camera_params_t *cameras,
					   int *added_order,
					   const std::vector<ImageKeyVector> 
					      &pt_views,
					   camera_params_t *camera_out,
					     bool remove_outliers);

    
    void MatchCloseImagesAndAddTracks(ImageData &data, int this_cam_idx,
                                      int added_order_idx, 
                                      std::vector<ImageKeyVector> &pt_views);
    void RunSFMWithNewImages(int new_images, 
                             double *S = NULL, double *U = NULL, \
                             double *V = NULL, double *W = NULL);
    void ReRunSFM(double *S = NULL, double *U = NULL, double *V = NULL, 
                  double *W = NULL);

    /* Run bundle adjustment on a given reconstruction */
    double RunSFM(int num_pts, int num_cameras, int start_camera,
                  bool fix_points, camera_params_t *init_camera_params,
                  v3_t *init_pts, int *added_order, v3_t *colors,
                  std::vector<ImageKeyVector> &pt_views, 
                  int max_iter = 0, int max_iter2 = 0, 
                  int verbosity = 0, double eps2 = 1.0e-12,
                  double *S = NULL, double *U = NULL, double *V = NULL,
                  double *W = NULL, bool remove_outliers = true,
                  bool final_bundle = false, 
                  bool write_intermediate = false);

    double RunSFM_SBA(int num_pts, int num_cameras, int start_camera,
                      bool fix_points, camera_params_t *init_camera_params,
                      v3_t *init_pts, int *added_order, v3_t *colors,
                      std::vector<ImageKeyVector> &pt_views, 
                      double eps2 = 1.0e-12,
                      double *S = NULL, double *U = NULL, double *V = NULL,
                      double *W = NULL, bool remove_outliers = true);

#ifdef __USE_CERES__
    double RunSFM_Ceres(int num_pts, int num_cameras, int start_camera,
                        bool fix_points, camera_params_t *init_camera_params,
                        v3_t *init_pts, int *added_order, v3_t *colors,
                        std::vector<ImageKeyVector> &pt_views, 
                        int max_iter = 0, int max_iter2 = 0, 
                        int verbosity = 0, double eps2 = 1.0e-12,
                        double *S = NULL, double *U = NULL, double *V = NULL,
                        double *W = NULL, bool remove_outliers = true,
                        bool final_bundle = false, 
                        bool write_intermediate = false);
#endif /* __USE_CERES__ */

    double RunSFMNecker(int i1, int i2, 
                        camera_params_t *cameras, 
                        int num_points, v3_t *points, v3_t *colors,
                        std::vector<ImageKeyVector> &pt_views,
                        camera_params_t *cameras_new,
                        v3_t *points_new, 
                        double threshold);


#endif /* __DEMO__ */

    bool BundleTwoFrame(int i1, int i2, TwoFrameModel *model, 
                        double &angle_out, int &num_pts_out, 
                        bool bundle_from_tracks);
    bool EstimateRelativePose(int i1, int i2, 
                              camera_params_t &camera1, 
                              camera_params_t &camera2);

    bool EstimateRelativePose2(int i1, int i2, 
                               camera_params_t &camera1, 
                               camera_params_t &camera2);

    /* Register a new image with the existing model */
    bool BundleRegisterImage(ImageData &data, bool init_location);
    void RunBundleServer();


#ifdef __USE_BOOST__
    /* Graph operations */
    ImageGraph ComputeMSTWorkingGraph(std::vector<int> &interior);
    void PartitionGraph(ImageGraph &graph, std::vector<int> interior);
#endif /* __USE_BOOST__ */

    /* Output a compressed version of the bundle file */
    void OutputCompressed(const char *ext = "compressed");
    /* Other operations on bundle files */
    void ScaleFocalLengths(double focal);    
    void ScaleFocalLengths(char *focal_file);
    void RotateCameras(char *rotate_file);
    void PruneBadPoints();
    void ZeroDistortionParams();
    void OutputRelativePoses2D(const char *outfile);
    void OutputRelativePoses3D(const char *outfile);
    void ComputeCameraCovariance();

    /* Analyze point statistics */
    void AnalyzePoints();

    /* Find a ground plane in the scene */
    void FindGroundPlane();

    /* Find a sky plane in the scene */
    void FindSkyPlane();

    /* Coalesce feature descriptors for each feature point */
    void CoalesceFeatureDescriptors();
    void CoalesceFeatureDescriptorsMedian();

    /* Compute likely matches between a set of keypoints and the
     * reconstructed points */
    std::vector<KeypointMatch>
        MatchKeysToPoints(const std::vector<KeypointWithDesc> &k1, 
                          double ratio = 0.6);

    std::vector<KeypointMatch>
        MatchPointsToKeys(const std::vector<KeypointWithDesc> &keys, 
			  double ratio = 0.6);

    void ReadProjectivePoints();
    void ReadProjectiveCameras();

    /* Routines for predicting the next images that should be captured */
    bool ImageVerifiesRay(int img, const double *p0, const double *p1);
    std::vector<int> GetVerifiersForImage(int img);
    void GetPointCoverage(int img, double &left, double &right,
                          double &up, double &down);
    int PredictNextImage(LinkDirection &dir);
    void RenderPredictedImage(int idx, LinkDirection dir, 
                              const char *out_file);




    /* **** Bundler Options **** */

    bool m_panorama_mode;        /* Are we reconstructing a panorama? */
    bool m_add_images_fast;
    bool m_estimate_ignored;
    bool m_analyze_matches;      /* Analyze matches */

    int m_ann_max_pts_visit;     /* Maximum points to visit during
                                  * global matching */

    bool m_optimize_for_fisheye; /* Optimize for fisheye-distorted
                                  * points */

    int m_homography_rounds;     /* Homography RANSAC params */
    double m_homography_threshold;

    int m_fmatrix_rounds;        /* F-matrix RANSAC params */
    double m_fmatrix_threshold;
    bool m_skip_fmatrix;
    bool m_skip_homographies;
    bool m_use_angular_score;

    double m_projection_estimation_threshold;  /* RANSAC threshold
						* for estimating 
						* projection matrix */

    double m_min_proj_error_threshold;
    double m_max_proj_error_threshold;

    double m_min_camera_distance_ratio;  /* The minimum distance for a
					  * non-panorama */

    double m_baseline_threshold;     /* The smallest permissible
				      * distance between two
				      * camera centers */

    double m_ray_angle_threshold;    /* Ray angle threshold */

    bool m_use_focal_estimate;       /* Estimate focal length of
				      * cameras */

    bool m_trust_focal_estimate;

    int m_min_max_matches;           /* Minimum number of matches
                                      * needed to register an image */
    int m_num_matches_add_camera;    /* Number of matches needed to 
                                      * consider registering an image */

    const char *m_bundle_output_file;  /* Output file names for BA */
    const char *m_bundle_output_base;
    const char *m_output_directory;

    bool m_compute_covariance;   /* Compute the covariance of a
                                  * reconstruction */
    int m_covariance_fix1;       /* Image to fix when computing
                                  * covariance */
    int m_covariance_fix2;       /* Image to fix translation of when
                                  * computing covariance */

    int m_keypoint_border_width; /* Throw out keypoints too close to
                                  * the border of an image */
    int m_keypoint_border_bottom; /* Throw out keypoints too close to
                                   * the bottom of an image */
    double m_init_focal_length;  /* Initial focal length for BA */
    bool m_fixed_focal_length;   /* Is the focal length constant? */
    bool m_use_constraints;      /* Should we use camera constraints? */
    bool m_constrain_focal;      /* Should we constrain the focal
				  * length of calibrated cameras? */
    double m_constrain_focal_weight;  /* The weight for the focal
				       * length constraint */

    bool m_factor_essential;     /* Should the model be initialized by
                                  * factoring the essential matrix? */

    bool m_estimate_distortion;  /* Should we estimate distortion for
				  * each camera? */
    double m_distortion_weight;  /* Weight on distortion parameter
                                  * constraints */

    bool m_construct_max_connectivity;  /* Do bundle adjustment using
					 * the connectivity score? */

    bool m_only_bundle_init_focal;  /* Only bundle adjust camera with
				     * initialized focal lengths */

    bool m_fix_necker;           /* Fix Necker reversal during bundle
				  * adjustment? */ 

    int m_initial_pair[2];       /* Images to use as the initial pair
				  * during bundle adjustment */

    bool m_features_coalesced;   /* Have features been coalesced */

    bool m_assemble;             /* Assemble the scene from the bottom up */
    bool m_run_bundle;           /* Should we run bundle adjustment
				  * automatically? */
    bool m_rerun_bundle;         /* Should we rerun bundle adjustment
				  * automatically? */
    bool m_fast_bundle;          /* Should we run the fast version of
				  * bundle adjustment? */

#ifdef __USE_CERES__
    bool m_use_ceres;            /* Use Ceres solver for 
                                  * bundle adjustment? */
#endif /* __USE_CERES__ */

    bool m_skip_full_bundle;     /* Skip full optimization stages */
    bool m_skip_add_points;      /* Don't add new points to the
                                  * optimization */

    /* Operations on bundle files */
    bool m_compress_list;        /* Output a compressed list and
				  * bundle file */
    bool m_reposition_scene;     /* Should we reposition the scene? */
    bool m_prune_bad_points;     /* Should we prune bad points? */



    double m_scale_focal;        /* Amount by which to scale the focal
                                  * lengths */

    bool m_predict_next_image;
    char *m_prediction_image;

    char *m_add_image_file;      /* Additional images to add */
    char *m_scale_focal_file;
    char *m_rotate_cameras_file;
    char *m_track_file;

    bool m_output_relposes;
    char *m_output_relposes_file;

    bool m_enrich_points;        /* Enrich the point set? */
    bool m_zero_distortion_params; /* Set all distortion parameters to
                                    * zero */

    int argc;
    char **argv;
};

#endif /* __bundlerapp_h__ */
