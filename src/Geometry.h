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

/* Geometry.h */
/* Geometric primitives */

#ifndef __geometry_h__
#define __geometry_h__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "ParameterBound.h"



#include "vector.h"

typedef std::pair<int,int> ImageKey;
typedef std::vector<ImageKey> ImageKeyVector;

/* Data for tracks */
class TrackData {
public:
    TrackData() : m_extra(-1) {}
    TrackData(ImageKeyVector views) : m_views(views), m_extra(-1) { }

    /* Read/write routines */
    void Read(FILE *f);
    void Write(FILE *f);

    ImageKeyVector m_views;
    int m_extra;
};

class PlaneData;

/* Data for 3D points */
class PointData {
public:
    PointData() { m_fixed = false; }

    /* Write the point data in XML */
    void WriteXML(FILE *f);
    void WriteGeoXML(FILE *f);

    /* Write coordinates*/
    void WriteCoordinates(FILE *f);

    /* Create a planar patch for this point */
    void CreatePlanarPatch(double size, PlaneData &plane);

    double m_pos[3];  /* 3D position of the point */
    double m_norm[3]; /* Estimated normal for this point */
    float m_color[3]; /* Color of the point */
    double m_conf;    /* Confidence in this point */

    ImageKeyVector m_views;  /* View / keys corresponding to this point */
    bool m_fixed;      /* Should this point be fixed during bundle
			* adjustment? */

    float *m_desc;     /* Descriptor for this point */
    int m_num_vis;     /*number of images that see this point*/
    int m_ref_image;   /* Reference image */
};

class LineSegmentMatch {
public:
    int m_idx1, m_idx2;
    double m_t1, m_u1;
    double m_t2, m_u2;
    // int m_extra;
};

class LineSegment2D {
public:
    double m_p1[2], m_p2[2];
    int m_extra;  /* Extra 4 bytes of info */

    /* I/O */
    void Read(FILE *f);
    void Write(FILE *f);

    /* Return the length of the segment */
    double Length();

    /* Returns true (and stores the endpoints in t and u) if the
     * epipolar swath e1, e2 intersects this segment */
    bool IntersectsEpipolarSwath(double e1[3], double e2[3], 
				 double &t, double &u);

    /* Return the sample at the given parameter value 
     * (l(0) = p1, l(1) = p2) */
    void Sample(double t, double *p);
    
    /* Return the homogeneous form of the line */
    void Homogeneous(double *l);
};

class CameraInfo;

class LineSegment3D {
public:
    LineSegment3D() : m_ignore(false) {}

    void Read(FILE *f);
    void Write(FILE *f);

    void Render(const CameraInfo &camera, double max_width,
		int stroke_texture, ParameterBound stroke_bounds);

    double m_p1[3], m_p2[3];
    std::vector<int> m_views;  /* Views this line is visible in */
    bool m_ignore;
};

class Cube {
public:
    Cube() {}
    

    Cube(double *origin, double *x_axis, double *y_axis, double *z_axis,
	 double x_scale, double y_scale, double z_scale) 
    {
	memcpy(m_origin, origin, 3 * sizeof(double));
	memcpy(m_x_axis, x_axis, 3 * sizeof(double));
	memcpy(m_y_axis, y_axis, 3 * sizeof(double));
	memcpy(m_z_axis, z_axis, 3 * sizeof(double));

	m_x_scale = x_scale;
	m_y_scale = y_scale;
	m_z_scale = z_scale;

	Finalize();
    }
    
    void Finalize();
    void Render();

    double m_origin[3];
    double m_x_axis[3], m_y_axis[3], m_z_axis[3];
    double m_x_scale, m_y_scale, m_z_scale;

    double m_vertices[24];
};

class PlaneMask {
public:
    double m_u_extent, m_v_extent;
    int m_grid_w, m_grid_h;
    bool *m_mask;
};

class PlaneData {
public:
    PlaneData() : m_texture_index(-1) {}
    PlaneData(double normal[3],double dist); //constructor

    /* Project a point onto the plane */
    void Project(double *p, double *p_proj);
    /* Project a vector onto the u or v axis */
    double ProjectU(double *p, double *p_proj);
    double ProjectV(double *p, double *p_proj);

    void Transform(const double *M);

    /* Setup various planar aspects */
    void Setup(std::vector<PointData> &point_data, 
	       double *origin, double *up);

    /* Intersect a ray with the plane */
    double IntersectRay(double *pos, double *ray, double *pt) const;

    /* Check if a point is inside the plane*/
    bool CheckInside(double *point);

    bool SetCorners(double *R, double *t, double f, int w, int h, 
                    int ymin=0, int ymax=-1, int xmin=0, int xmax=-1);



    /* Read plane from file */
    void Read(FILE *f);
    void Write(FILE *f);

    /* Render the plane using OpenGL */
    void Render();

    /* Return the plane parameters */
    void GetParams(double *params);

    double m_normal[3];

    double m_origin[3];
    double m_uaxis[3];
    double m_vaxis[3];
    double m_corners[12];

    double m_color[4];

    double m_u_extent, m_v_extent;

    double m_dist;

    std::vector<int> m_views;  /* Views corresponding to this plane */
    std::vector<int> m_points; /* Points corresponding to this plane */

    double m_front[3];   /* Vector pointing out of the front of the
			  * plane */

    PlaneMask m_mask; /* Mask for this plane */

    int m_texture_index;
    ParameterBound m_bounds;
};

/* Return a unit vector pointing in the direction of the given line */
void LineToUnitVector(double *l, double *v);

v3_t FindRobustMean(const std::vector<v3_t> &points);

v3_t FindWeightedRobustMean(const std::vector<v3_t> &points,
			    const std::vector<double> &weights);

double FindRobustVariance(v3_t mean, const std::vector<v3_t> &points);

/* Fit a plane to the points at the given indices */
std::vector<int> FitPlaneToPoints(const std::vector<PointData> &points,
                                  const std::vector<int> &indices, 
                                  double *plane, 
                                  int ransac_rounds,
                                  double ransac_threshold,
                                  bool par_to_up, bool perp_to_up,
                                  double *up = NULL);

/* Return the projections of a set of points into a given camera */
std::vector<v2_t> GetPointProjections(const CameraInfo &cam, 
                                      const std::vector<PointData> &points,
                                      const std::vector<int> &indices,
                                      bool inside_only, 
                                      int &num_inside);

#endif /* __geometry_h__ */
