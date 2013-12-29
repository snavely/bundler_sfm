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

/* ImageData.h */
/* Simple image storage and operations */

#ifndef __image_data_h__
#define __image_data_h__
#include <vector>
using namespace std;

#ifndef __BUNDLER__
#include "wx/wx.h"
#endif /* __BUNDLER__ */

#include "BoundingBox.h"
#include "Camera.h"
#include "Geometry.h"
#include "ParameterBound.h"
// #include "Polygon.h"
#include "keys.h"



#include "image.h"



enum ImageBoundary {
    BoundaryNorth,
    BoundaryEast,
    BoundarySouth,
    BoundaryWest,
};

class ImageNote {
public:
    void Read(FILE *f, double width, double height);

    char *m_text;
    BoundingBox m_bbox;
};

class ImageDate {
public:
    ImageDate() : m_known(false) { }

    void GetString(char *buf);
    void GetMonthString(char *buf);

    /* Return the date in the form of a double */
    double GetDateDouble();
    double GetDateOnlyDouble();
    /* Return the time in the form of a double */
    double GetTimeDouble();

    bool m_known;  /* Is the date/time known? */
    int m_month, m_day, m_year;      /* Date */
    int m_hour, m_minute, m_second;  /* Time */
};

int CompareImageDates(ImageDate *date1, ImageDate *date2);

class AnnotationData {
public:
    AnnotationData(int idx, BoundingBox bbox) : m_idx(idx), m_bbox(bbox) { }
    AnnotationData() { }

    int m_idx;
    BoundingBox m_bbox;
};

#ifndef __BUNDLER__
class TPSBasis {
public:
    TPSBasis() {
        m_LU = NULL;
        m_ipiv = NULL;
    }

    void Clear() {
        if (m_LU != NULL) {
            delete [] m_LU;
            m_LU = NULL;
        }

        if (m_ipiv != NULL) {
            delete [] m_ipiv;
            m_ipiv = NULL;
        }

        m_basis_points.clear();
        m_indices.clear();
    }

    double *m_LU;
    int *m_ipiv;
    std::vector<DPoint> m_basis_points;
    std::vector<int> m_indices;
};
#endif /* __BUNDLER__ */

class ImageData {
public:
    ImageData() { 
        m_img = NULL;
        m_thumb_fixed = NULL;
        m_thumb = NULL;
        m_thumb256 = NULL;


        m_marked = -1;
        m_canonical=false;
        m_canonical_pano=false;
        m_rahul_pano_index=-1;
        m_texture_index = m_back_texture_index = -1;
        m_thumb_texture_index = -1;
        m_thumb_fixed_texture_index = -1;
        m_texture_image_loaded = false;
        m_panorama_index = -1;
        m_keys_loaded = false; 
        m_keys_desc_loaded = false;
        m_keys_scale_rot_loaded = false;
        m_cached_dimensions = false; 
        m_cached_keys = false;
        //m_k0 = 1.0;
        //m_k1 = m_k2 = m_k3 = m_k4 = 0.0; 
        m_known_intrinsics = false;
        //m_rd_focal = 0.0;
        m_init_focal = 0.0;
        m_has_init_focal = false;
        m_rotation = 0;
        m_is_dragged = false;
        m_is_dropped = false;
        m_ignore_in_bundle = false;
        m_licensed = false;
        m_real_name = NULL;
        m_partition = -1;
        m_layout_x = m_layout_y = 0;
        m_lat = m_long = 0.0;
        m_day_photo = true;
        m_added = false;
        m_is_rep = false;
        m_geosupport = 0.0;
    }

    /* Initialize the image data given a string description */
    void InitFromString(char *buf, const char *path, 
                        bool fisheye_by_default);

    /* Create a pinhole view for the keys */
    void UndistortKeys();
    std::vector<Keypoint> UndistortKeysCopy();
    void DistortPoint(double x, double y, double *R, 
        double &x_out, double &y_out) const;
    void DistortPoint(double x, double y, double &x_out, double &y_out) const;
    void UndistortPoint(double x, double y, 
        double &x_out, double &y_out) const;

    /* Radial distortion routines */
    void DistortPointRD(double x, double y, 
        double &x_out, double &y_out) const;
    void UndistortPointRD(double x, double y, 
        double &x_out, double &y_out) const;

    img_t *UndistortImage(double theta = 0.0, double phi = 0.0, 
        int num_rot = 0);
    /* Undistort the image to a different size than the original */
    img_t *UndistortImageResize(int w_new, int h_new);

    void ToNDC(double x, double y, double &x_n, double &y_n);
    void FromNDC(double x_n, double y_n, double &x, double &y);

    /* Save the dimensions of the image */
    void CacheDimensions();
    void CacheNumKeys();

    /* Return the image dimensions */
    int GetWidth();
    int GetHeight();
    int GetArea();

    void GetBaseName(char *buf);    
    void GetName(char *buf);
    void GetCompressedTextureFilename(char *buf) const;
    bool CompressedTextureExists();

    /* Return the bounding box of the image */
    BoundingBox GetBoundingBox();

    /* Get the texture coordinates for a point */
    v2_t GetPointTexCoords(v2_t point);

    /* Get the texture coordinates of a set of points */
    std::vector<v2_t> GetPointTexCoords(const std::vector<v2_t> &pts);

    /* Returns true if the given pixel is in range */
    bool PixelInRange(double x, double y);

    void LoadImage();
    void UnloadImage();

    void LoadTexImage();
    void UnloadTexImage();

	void LoadBackTexImage();
	void UnloadBackTexImage();
	
	bool TexImageExists();

    /* Compute a Gaussian pyramid for the current image */

    bool LoadFeatureWeightMap();
    void SaveFeatureWeightMap();
    void ComputeFeatureWeightMap(int id, 
        const std::vector<PointData> &pt_data);
    void ComputeFeatureWeightMap(int id, 
        const std::vector<PointData> &pt_data,
        const std::vector<int> &used_points);

    void CheckLoadFloatingThumb();
    void LoadFloatingThumbnail();
    void UnloadFloatingThumbnail();

    void CheckLoadFixedThumb();
    void LoadFixedThumbnail(int w_max, int h_max, int rotation);
    void UnloadFixedThumbnail();

    void LoadThumb256();
    void UnloadThumb256();

    int GetNumKeys();
    void LoadOrExtractKeys(const char *sift_binary, bool undistort = true);
    void LoadKeys(bool descriptor = true, bool undistort = true);
    void LoadDescriptors(bool undistort);
    void LoadKeysWithScaleRot(bool descriptor = true, bool undistort = true);
    void UnloadKeys();
    void UnloadKeysWithScaleRot();
    void ExtractFeatures(const char *sift_binary, bool undistort);

    /* Find line segments in the image */
    void DetectLineSegments(double sigma, 
        double threshold1, double threshold2,
        double min_line_segment_size);

    /* Render line segments to an image */
    img_t *RenderLineSegments();

    /* Query whether a given line segment is supported by the image */
    bool LineSegmentSupportedByImage(LineSegment2D &line);

    /* Return true if the given line intersects the image */
    bool LineIntersectsImage(double *line, double *isect1, double *isect2,
        ImageBoundary &b1, ImageBoundary &b2);

    /* Find the gradient direction at the given position */
    void Gradient(double x, double y, double *grad);

    vector<bool> m_visinfo;
    vector<int> m_visindex;



    /* Read key colors */
    void ReadKeyColors();

    /* Read the camera */
    bool ReadCamera();
    /* Read the tracks */
    bool ReadTracks(int img_idx, std::vector<PointData> &pt_list);
    /* Write the camera */
    void WriteCamera();
    /* Write the camera in XML format */
    void WriteCameraXML(FILE *f);

    /* Write the tracks */
    void WriteTracks();

    void ReadMetadata();

    /* Get notes associated with this image */
    std::vector<ImageNote> GetNotes();

#ifndef __BUNDLER__
    /* Compute 3D points used for estimating the amount of distortion
     * in the image */
    void ComputeDistortionPoints(const std::vector<PointData> &pt_data,
                                 std::vector<DPoint3> &pts_world,
                                 std::vector<DPoint3> &pts_plane);

    void ComputeOrthoPlane(const std::vector<PointData> &pt_data);

    void ComputeTPSBasis(const std::vector<PointData> &pt_data);
    void FreeTPSBasis();
#endif /* __BUNDLER__ */

    void ComputeHistogram(int nbins, double *r_hist, 
        double *g_hist, double *b_hist);

    void ComputeHistogramLUV(int nbins, double *L_hist, 
        double *U_hist, double *V_hist);

    void ComputeHistogramLUVNorm(int nbins, double *L_hist, 
        double *U_hist, double *V_hist);

    int m_width, m_height;   /* Cached dimensions */
    int m_xmin,m_xmax,m_ymin,m_ymax; /*dimensions of bounding box of geometry projected onto the image*/
    bool m_cached_dimensions;

    int m_num_keys;          /* Cached number of keys */
    bool m_cached_keys;

    bool m_added;
    int m_marked; /* if true, then the camera is on some orbit*/
    bool m_canonical;
    bool m_canonical_pano;
    int m_rahul_pano_index;

    char *m_real_name;
    char *m_name;             /* Filename */
    char *m_key_name;         /* Key filename */
    char m_user_name[256];    /* User name */
    char m_flickr_index[256]; /* Flickr index */
    ImageDate m_date;         /* Date / time the photo was taken */
    bool m_day_photo;         /* Was this photo shot during the day? */

#ifndef __BUNDLER__
    wxImage *m_wximage;       /* WX image */
    wxBitmap *m_bitmap;       /* Bitmap */
#endif /* __BUNDLER__ */

    img_t *m_img;             /* Image */
    img_t *m_texture_img;     /* Texture image  for foreground -- rahul*/
	img_t *m_back_texture_img;/* Texture Image for background -- rahul*/


    img_t *m_thumb;           /* Thumbnail 1 */
    img_t *m_thumb8;          /* Thumbnail 2 */
    img_t *m_thumb_fixed;     /* Thumbnail 3 */
    img_t *m_thumb256;        /* Thumbnail 4 */
    img_t *m_feature_weight_map; /* Weighting for image comparison */



    bool m_image_loaded, m_keys_loaded, m_keys_desc_loaded, 
        m_keys_scale_rot_loaded;
    bool m_texture_image_loaded;
    bool m_licensed;         /* Can we actually use this image? */

    bool m_fisheye;            /* Is this a fisheye image? */
    double m_fCx, m_fCy;       /* Fisheye center */
    double m_fRad, m_fAngle;   /* Other fisheye parameters */
    double m_fFocal;           /* Desired focal length */

    double m_lat, m_long;      /* Latitude and longitude (if known) */
    double m_geocentric[3];    /* Geocentric position */
    double m_geoplanar[3];     /* Geoplanar */
    double m_geosupport;       /* What proportion of pairs vouch for me? */

    std::vector<AnnotationData> m_notes;

    /* Radial distortion parameters */
    bool m_known_intrinsics;
    double m_K[9];
    //double m_k0, m_k1, m_k2, m_k3, m_k4, m_rd_focal;
    double m_k[5];	

    bool m_ignore_in_bundle;  /* Ignore this image during bundle
                              * adjustment */

    /* Focal length parameters */
    bool m_has_init_focal;   /* Do we have an initial focal length? */
    double m_init_focal;     /* Initial focal length */

    CameraInfo m_camera;     /* Information on the camera used to 
                             * capture this image */

    std::vector<LineSegment2D> m_line_segments;

    /* Texture mapping data */
    int m_texture_index;               /* Texture index of this image */
	int m_back_texture_index;		   /* Texture index OF background image*/
    int m_thumb_texture_index;         /* Texture index of the thumbnail */
    int m_thumb_fixed_texture_index;   /* Texture index of the thumbnail */

    ParameterBound m_bounds;  /* Texture coord bounds */
    ParameterBound m_thumb_bounds;
    ParameterBound m_thumb_fixed_bounds;

    PlaneData m_fit_plane;    /* Plane fit to the observed points */
    // PlaneData m_ortho_plane;
    PlaneData m_current_plane;  /* Temporary plane valid for current
                                * view */

    std::vector<Keypoint> m_keys;              /* Keypoints in this image */
    std::vector<KeypointWithDesc> m_keys_desc; /* Keypoints with descriptors */
    std::vector<KeypointWithScaleRot> m_keys_scale_rot;
    std::vector<bool> m_key_flags;

    std::vector<int> m_visible_points;  /* Indices of points visible
                                         * in this image */
    std::vector<int> m_visible_keys;

    std::vector<int> m_visible_lines;   /* Indices of lines visible 
                                         * in this image */

#ifndef __BUNDLER__
    std::vector<DPoint3> m_dist_pts_world;
    std::vector<DPoint3> m_dist_pts_plane;
#endif /* __BUNDLER__ */

    std::vector<int> m_neighbors;       /* List of neighboring images */

    int m_rotation;                     /* Orientation of this image
                                         * (0,1,2,3) */

    int m_panorama_index;               /* Index of the panorama this
                                         * view belongs to */

    bool m_is_dragged;                  /* Is this image being dragged? */
    bool m_is_dropped;                  /* Has this image been dropped? */

    bool m_is_rep;                      /* Is this a representative image? */

    int m_drag_x, m_drag_y;             /* Screen coords of drag
                                         * position */

    double m_drop_pt[3];                /* World coords of drop
                                         * position */

    int m_partition;                    /* Partition this image belongs to */
    int m_layout_x, m_layout_y;         /* 2D layout positions */

#ifndef __BUNDLER__
    TPSBasis m_tps_basis;               /* Thin-plate spline basis */
#endif
};

#endif /* __image_data_h__ */
