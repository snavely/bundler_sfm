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

/* Camera.h */
/* Camera classes */

#ifndef __camera_h__
#define __camera_h__

#ifndef __BUNDLER__
#include "Bezier.h"
#include "LinkDirection.h"
#endif

#include "BoundingBox.h"
#include "Geometry.h"

class CameraInfo {
public:
    CameraInfo() { 
        m_adjusted = false; 

        m_constrained[0] = m_constrained[1] = m_constrained[2] = false; 
        m_constrained[3] = m_constrained[4] = m_constrained[5] = false; 
        m_constrained[6] = false;

        m_constraints[0] = m_constraints[1] = m_constraints[2] = 0.0;
        m_constraints[3] = m_constraints[4] = m_constraints[5] = 0.0;
        m_constraints[6] = 0.0;

        m_constraint_weights[0] = 
            m_constraint_weights[1] = 
            m_constraint_weights[2] = 0.0;
        m_constraint_weights[3] = 
            m_constraint_weights[4] = 
            m_constraint_weights[5] = 0.0;
        m_constraint_weights[6] = 0.0;      

        m_k[0] = m_k[1] = 0.0;
        
#ifndef __BUNDLER__
        for (int i = 0; i < NUM_LINK_DIRECTIONS; i++) 
            m_links[i] = -1;
#endif
    }

    /* Finalize the camera */
    void Finalize();
    /* Get rigid transform for this camera */
    void GetRigid(double *T) const;
    /* Get a 4x4 rigid transform */
    void GetRigid4x4(double *T) const;
    /* Return the position of the camera */
    void inline GetPosition(double *pos) const {
        pos[0] = m_R[0] * m_t[0] + m_R[3] * m_t[1] + m_R[6] * m_t[2];
        pos[1] = m_R[1] * m_t[0] + m_R[4] * m_t[1] + m_R[7] * m_t[2];
        pos[2] = m_R[2] * m_t[0] + m_R[5] * m_t[1] + m_R[8] * m_t[2];

        pos[0] = -pos[0];
        pos[1] = -pos[1];
        pos[2] = -pos[2];
    }
    
    /* Set the position of the camera */
    void SetPosition(const double *pos);
    /* Return the pose of the camera as a rotation matrix */
    void GetPose(double *R) const;
    void GetPoseQuaternion(double *q) const;
    /* Set the pose of the camera */
    void SetPose(const double *R);
    /* Get upright rotation matrix */
    void GetUprightRotation(int rotate, double *R);
    /* Return the 3x3 intrinsic matrix */
    void GetIntrinsics(double *K) const;
    /* Get the field of view */
    double GetFOV() const;
    double GetFOVMax(int rotate) const;
    /* Set the field of view */
    void SetFOV(double fov);
    /* Project the point into the camera */
    bool Project(const double *p, double *proj) const;
    /* Compute the essential matrix between this camera and another one */
    void ComputeEssentialMatrix(const CameraInfo &cam, double *E, double *F);
    /* Flip the camera over the z-axis */
    void Reflect();
    /* Compute the distance to another camera */
    double CameraDistance(const CameraInfo &cam) const;
    /* Find the horizon line in the image */
    void ComputeHorizonLine(double *ground, double *up);
    /* Compute vanishing line */
    void ComputeVanishingLine(const PlaneData &plane, double *line);
    /* Returns true if the given point is in front of the camera */
    bool PointInFront(double *p);
    /* Returns true if the given 2D point is above the horizon */
    bool PointAboveHorizon(double *p);

#ifndef __BUNDLER__
    /* Return a bezier curve given another camera */
    Bezier ComputeBezier(const CameraInfo &cam) const;
#endif

    /* Convert a pixel position to a ray direction */
    void PixelToCameraRay(double x, double y, double *ray);
    /* Convert a pixel position to an (absolute) ray direction */
    void PixelToCameraRayAbsolute(double x, double y, double *ray);
    /* Point the camera in a different direction */
    void PointAt(double *ray);
    void PointAtAbsolute(double *ray);
    /* Get the bounding box for this camera */
    BoundingBox GetBoundingBox() const;
    /* Return the view direction */
    void GetViewDirection(double *view) const;
    /* Get the twist angle of the camera */
    double GetTwistAngleRadians() const;
    /* Return the halfspace in front of the camera */
    void GetFrontHalfspace(double *plane) const;
    bool PointInsideImage(const double *p) const;

    /* Get a camera whose up vector points in the right direction */
    CameraInfo GetUpCamera(double *up) const;

#ifndef __BUNDLER__
    /* Read/write the links to a file */
    void ReadLinks(FILE *f);
    void WriteLinks(FILE *f);
#endif

    /* Write in XML format */
    void WriteXML(FILE *f);



    /* Write params to file*/
    void WriteFile(FILE *f);

    bool m_adjusted;        /* Has this camera been adjusted? */
    double m_focal;         /* Focal length */
    double m_k[2];          /* Distortion parameters */
    double m_R[9], m_t[3];  /* Extrinsics */
    double m_Pmatrix[12];
    int m_width, m_height;  /* Image dimensions */

    /* Horizon line */
    double m_horizon[3];

    double m_RGB_transform[12];   /* Local affine transform for RGB
                                   * space */

    double m_color[3];

    /* Constraints on camera center location */
    bool m_constrained[7];
    double m_constraints[7];
    double m_constraint_weights[7];

#ifndef __BUNDLER__
    int m_links[NUM_LINK_DIRECTIONS];
#endif /* __BUNDLER__ */
};

/* Interpolate between two camera views */
CameraInfo InterpolateCameras(const CameraInfo &cam1, 
                              const CameraInfo &cam2, double t);

CameraInfo InterpolateCamerasThetaPhi(const CameraInfo &cam1, 
                                      const CameraInfo &cam2, double t,
                                      bool interp_fov = false);

#endif /* __camera_h__ */
