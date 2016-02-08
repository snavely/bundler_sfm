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

/* Camera.cpp */
/* Camera classes */

#include <string.h>
#include <math.h>

#include "Camera.h"
#include "BundleUtil.h"
#include "Geometry.h"

#include "defines.h"
#include "fit.h"
#include "matrix.h"
#include "vector.h"



/* Finalize the camera */
void CameraInfo::Finalize() {
    /* Compute projection matrix */
    double K[9];
    double Ptmp[12] = 
	{ m_R[0], m_R[1], m_R[2], m_t[0],
	  m_R[3], m_R[4], m_R[5], m_t[1],
	  m_R[6], m_R[7], m_R[8], m_t[2] };
	    
    GetIntrinsics(K);

    matrix_product(3, 3, 3, 4, K, Ptmp, m_Pmatrix);
}

/* Get rigid transform for this camera */
void CameraInfo::GetRigid(double *T) const
{
    T[0] = m_R[0];  T[1] = m_R[1];  T[2] = m_R[2];  T[3] = m_t[0];
    T[4] = m_R[3];  T[5] = m_R[4];  T[6] = m_R[5];  T[7] = m_t[1];
    T[8] = m_R[6];  T[9] = m_R[7];  T[10] = m_R[8]; T[11] = m_t[2];
}

/* Get a 4x4 rigid transform */
void CameraInfo::GetRigid4x4(double *T) const {
    GetRigid(T);
    T[12] = T[13] = T[14] = 0.0;
    T[15] = 1.0;
}

/* Return the position of the camera */
#if 0
void inline CameraInfo::GetPosition(double *pos) const {
    // double Rinv[9];

    // matrix_transpose(3, 3, (double *) m_R, Rinv);
    // matrix_transpose_product(3, 3, 3, 1, (double *) m_R, (double *) m_t, pos);

    pos[0] = m_R[0] * m_t[0] + m_R[3] * m_t[1] + m_R[6] * m_t[2];
    pos[1] = m_R[1] * m_t[0] + m_R[4] * m_t[1] + m_R[7] * m_t[2];
    pos[2] = m_R[2] * m_t[0] + m_R[5] * m_t[1] + m_R[8] * m_t[2];

    pos[0] = -pos[0];
    pos[1] = -pos[1];
    pos[2] = -pos[2];
}
#endif

/* Set the position of the camera */
void CameraInfo::SetPosition(const double *pos) {
    matrix_product(3, 3, 3, 1, m_R, (double *) pos, m_t);
    matrix_scale(3, 1, m_t, -1.0, m_t);
}

/* Return the pose of the camera as a rotation matrix */
void CameraInfo::GetPose(double *R) const {
    matrix_transpose(3, 3, (double *) m_R, R);
}

void CameraInfo::GetPoseQuaternion(double *q) const {
    double R[9];
    matrix_transpose(3, 3, (double *) m_R, R);
    matrix_to_quaternion(R, q);
}

/* Set the pose of the camera */
void CameraInfo::SetPose(const double *R) {
    matrix_transpose(3, 3, (double *) R, (double *) m_R);
}

/* Get upright rotation matrix */
void CameraInfo::GetUprightRotation(int rotate, double *R) {
    double R90[9] = 
        { 0.0, -1.0,  0.0,
          1.0,  0.0,  0.0,
          0.0,  0.0,  1.0 };
    
    double Rup[9];
    matrix_power(3, R90, rotate, Rup);

    matrix_product(3, 3, 3, 3, Rup, m_R, R);
}

/* Return the 3x3 intrinsic matrix */
void CameraInfo::GetIntrinsics(double *K) const {
    K[0] = m_focal; K[1] = 0.0;     K[2] = 0.0;
    K[3] = 0.0;     K[4] = m_focal; K[5] = 0.0;
    K[6] = 0.0;     K[7] = 0.0;     K[8] = 1.0;
}

/* Get the field of view */
double CameraInfo::GetFOV() const {
    return 2.0 * atan(m_width / (2.0 * m_focal));
}

double CameraInfo::GetFOVMax(int rotate) const {
    if (((rotate % 2) == 0 && m_width >= m_height) || 
        ((rotate % 2) == 1 && m_width < m_height)) {
        return 2.0 * atan(m_width / (2.0 * m_focal));
    } else {
        double vfov = 2.0 * atan(m_height / (2.0 * m_focal));
        double hfov = 2.0 * atan(tan(0.5 * vfov) * m_width / (double) m_height);
        printf("vfov = %0.3f, hfov = %0.3f\n", vfov, hfov);
        return hfov;
    }
}

/* Set the field of view */
void CameraInfo::SetFOV(double fov) {
    m_focal = 0.5 * m_width / tan(0.5 * DEG2RAD(fov));
}

/* Project the point into the camera */
bool CameraInfo::Project(const double *p, double *proj) const
{
    double p4[4] = { p[0], p[1], p[2], 1.0 };
    double proj3[3];

    matrix_product(3, 4, 4, 1, (double *) m_Pmatrix, p4, proj3);

    if (proj3[2] == 0.0)
        return false;

    proj[0] = proj3[0] / -proj3[2];
    proj[1] = proj3[1] / -proj3[2];

    if (m_k[0] == 0.0 && m_k[1] == 0.0)
        return (proj3[2] < 0.0);

    double rsq = (proj[0] * proj[0] + proj[1] * proj[1]) / (m_focal * m_focal);
    double factor = 1.0 + m_k[0] * rsq + m_k[1] * rsq * rsq;

    if (rsq > 8.0 || factor < 0.0) // bad extrapolation
        return (proj[2] < 0.0);

    proj[0] *= factor;
    proj[1] *= factor;
    
    return (proj3[2] < 0.0);
}

/* Compute the essential matrix between this camera and another one */
void CameraInfo::ComputeEssentialMatrix(const CameraInfo &cam, 
					double *E, double *F)
{
    /* Put the first camera at the canonical location */
    double P[16], Pinv[16];
    
    // memcpy(P, m_Pmatrix, sizeof(double) * 12);
    GetRigid(P);
    P[12] = P[13] = P[14] = 0.0;  P[15] = 1.0;
    
    matrix_invert(4, P, Pinv);

    double P2[16];
    cam.GetRigid(P2);
    P2[12] = P2[13] = P2[14] = 0.0;  P2[15] = 1.0;

    double P2new[16];
    matrix_product(4, 4, 4, 4, P2, Pinv, P2new);

    double R2new[9] = { P2new[0], P2new[1], P2new[2],
			P2new[4], P2new[5], P2new[6],
			P2new[8], P2new[9], P2new[10] };

    double t2new[3] = { P2new[3], P2new[7], P2new[11] };
    double t2new_cross[9] = { 0.0, -t2new[2], t2new[1],
			      t2new[2], 0.0, -t2new[0],
			      -t2new[1], t2new[0], 0.0 };

    matrix_product(3, 3, 3, 3, t2new_cross, R2new, E);

    /* Special black magic because we flipped the Z-axis */
    E[0] = -E[0];
    E[1] = -E[1];
    E[3] = -E[3];
    E[4] = -E[4];
    E[8] = -E[8];

    /* Compute F matrix */
    double K1[9], K2[9];
    GetIntrinsics(K1);
    cam.GetIntrinsics(K2);

    double K1inv[9], K2inv[9];
    matrix_invert(3, K1, K1inv);
    matrix_invert(3, K2, K2inv);

    double tmp[9];
    matrix_transpose_product(3, 3, 3, 3, K2inv, E, tmp);
    matrix_product(3, 3, 3, 3, tmp, K1inv, F);
}

/* Flip the camera over the z-axis */
void CameraInfo::Reflect() 
{
    m_R[2] = -m_R[2];
    m_R[5] = -m_R[5];
    m_R[6] = -m_R[6];
    m_R[7] = -m_R[7];

    m_t[2] = -m_t[2];

    Finalize();
}

/* Compute the distance to another camera */
double CameraInfo::CameraDistance(const CameraInfo &cam) const 
{
    double pos1[3], pos2[3];
    
    GetPosition(pos1);
    cam.GetPosition(pos2);

    double diff[3];
    matrix_diff(3, 1, 3, 1, pos1, pos2, diff);
    
    return matrix_norm(3, 1, diff);
}

/* Compute vanishing line */
void CameraInfo::ComputeVanishingLine(const PlaneData &plane, double *line)
{
    /* Find two points on the intersection of the plane with the plane
     * at infinity */

    double normal[3];
    memcpy(normal, plane.m_normal, 3 * sizeof(double));
    
    matrix_scale(3, 1, normal, 1.0 / matrix_norm(3, 1, normal), normal);

    double z1 = (-normal[0] - normal[1]) / normal[2];
    double z2 = (normal[0] - normal[1]) / normal[2];
    
    double p1[4] = {  1.0,  1.0, z1, 0.0 };
    double p2[4] = { -1.0,  1.0, z2, 0.0 };

    double p1_rot[3], p2_rot[3];

    /* Project the point into the image */
    matrix_product(3, 3, 3, 1, m_R, p1, p1_rot);
    matrix_product(3, 3, 3, 1, m_R, p2, p2_rot);

    double pr1[3] = { m_focal * p1_rot[0] / p1_rot[2],
		      m_focal * p1_rot[1] / p1_rot[2], -1.0 };

    double pr2[3] = { m_focal * p2_rot[0] / p2_rot[2],
		      m_focal * p2_rot[1] / p2_rot[2], -1.0 };

    matrix_cross(pr1, pr2, line);
}


/* Find the horizon line in the image */
void CameraInfo::ComputeHorizonLine(double *ground, double *up)
{
    double pos[3];
    GetPosition(pos);

    double R[9];
    GetPose(R);

    // double pos4[4] = { pos[0], pos[1], pos[2], 1.0 };

    /* Find the 3D positions of the corners of the image */
    double  upper_left[4] = { -0.5 * m_width,  0.5 * m_height, -m_focal, 1.0 };
    double  lower_left[4] = { -0.5 * m_width, -0.5 * m_height, -m_focal, 1.0 };
    double upper_right[4] = {  0.5 * m_width,  0.5 * m_height, -m_focal, 1.0 };
    double lower_right[4] = {  0.5 * m_width, -0.5 * m_height, -m_focal, 1.0 };
    
    matrix_scale(3, 1, upper_left, 1.0 / m_focal, upper_left);
    matrix_scale(3, 1, lower_left, 1.0 / m_focal, lower_left);
    matrix_scale(3, 1, upper_right, 1.0 / m_focal, upper_right);
    matrix_scale(3, 1, lower_right, 1.0 / m_focal, lower_right);

    /* Convert to world space coordinate */
    double rigid[16];
    GetRigid4x4(rigid);

    double rigid_inv[16];
    matrix_invert(4, rigid, rigid_inv);

    double upper_left_world[4], lower_left_world[4];
    double upper_right_world[4], lower_right_world[4];

    matrix_product(4, 4, 4, 1, rigid_inv, upper_left, upper_left_world);
    matrix_product(4, 4, 4, 1, rigid_inv, lower_left, lower_left_world);
    matrix_product(4, 4, 4, 1, rigid_inv, upper_right, upper_right_world);
    matrix_product(4, 4, 4, 1, rigid_inv, lower_right, lower_right_world);

    /* Compute planes for the "sides" of the viewing frustrum */
    double minuz[3] = { 0.0, 0.0, -1.0 }, dir[3];
    matrix_product(3, 3, 3, 1, R, minuz, dir);

    double infinity[4] = { dir[0], dir[1], dir[2], 10000.0 };
    double left_plane[4], right_plane[4], upper_plane[4], lower_plane[4];
    
    // matrix_cross4(pos4, upper_left_world, lower_left_world, left_plane);
    // matrix_cross4(pos4, upper_right_world, lower_right_world, right_plane);
    v3_t left_points[3] = { v3_new(pos[0], pos[1], pos[2]),
			    v3_new(upper_left_world[0],
				   upper_left_world[1],
				   upper_left_world[2]),
			    v3_new(lower_left_world[0], 
				   lower_left_world[1],
				   lower_left_world[2]) };
    
    v3_t right_points[3] = { v3_new(pos[0], pos[1], pos[2]),
			     v3_new(upper_right_world[0],
				    upper_right_world[1],
				    upper_right_world[2]),
			     v3_new(lower_right_world[0], 
				    lower_right_world[1],
				    lower_right_world[2]) };

    v3_t upper_points[3] = { v3_new(pos[0], pos[1], pos[2]),
			     v3_new(upper_left_world[0],
				    upper_left_world[1],
				    upper_left_world[2]),
			     v3_new(upper_right_world[0], 
				    upper_right_world[1],
				    upper_right_world[2]) };

    v3_t lower_points[3] = { v3_new(pos[0], pos[1], pos[2]),
			     v3_new(lower_left_world[0],
				    lower_left_world[1],
				    lower_left_world[2]),
			     v3_new(lower_right_world[0], 
				    lower_right_world[1],
				    lower_right_world[2]) };

    fit_3D_plane_orthogonal_regression(3, left_points, left_plane);
    fit_3D_plane_orthogonal_regression(3, right_points, right_plane);
    fit_3D_plane_orthogonal_regression(3, lower_points, upper_plane);
    fit_3D_plane_orthogonal_regression(3, upper_points, lower_plane);

    matrix_scale(4, 1, infinity, 1.0 / infinity[3], infinity);
    matrix_scale(4, 1, ground, 1.0 / ground[3], ground);
    matrix_scale(4, 1, left_plane, 1.0 / left_plane[3], left_plane);
    matrix_scale(4, 1, right_plane, 1.0 / right_plane[3], right_plane);
    matrix_scale(4, 1, upper_plane, 1.0 / upper_plane[3], upper_plane);
    matrix_scale(4, 1, lower_plane, 1.0 / lower_plane[3], lower_plane);

    /* Cross three planes to find the projection of the line at
     * infinity */
    v3_t left_planes[3] = 
	{ v3_new(infinity[0], infinity[1], infinity[2]),
	  v3_new(ground[0], ground[1], ground[2]),
	  v3_new(left_plane[0], left_plane[1], left_plane[2]) };

    v3_t right_planes[3] = 
	{ v3_new(infinity[0], infinity[1], infinity[2]),
	  v3_new(ground[0], ground[1], ground[2]),
	  v3_new(right_plane[0], right_plane[1], right_plane[2]) };

    v3_t upper_planes[3] = 
	{ v3_new(infinity[0], infinity[1], infinity[2]),
	  v3_new(ground[0], ground[1], ground[2]),
	  v3_new(upper_plane[0], upper_plane[1], upper_plane[2]) };

    v3_t lower_planes[3] = 
	{ v3_new(infinity[0], infinity[1], infinity[2]),
	  v3_new(ground[0], ground[1], ground[2]),
	  v3_new(lower_plane[0], lower_plane[1], lower_plane[2]) };
    
    double point_left[4], point_right[4], point_upper[4], point_lower[4];
    fit_3D_plane_orthogonal_regression(3, left_planes, point_left);
    fit_3D_plane_orthogonal_regression(3, right_planes, point_right);
    fit_3D_plane_orthogonal_regression(3, upper_planes, point_upper);
    fit_3D_plane_orthogonal_regression(3, lower_planes, point_lower);

    matrix_scale(4, 1, point_left, 1.0 / point_left[3], point_left);
    matrix_scale(4, 1, point_right, 1.0 / point_right[3], point_right);
    matrix_scale(4, 1, point_upper, 1.0 / point_upper[3], point_upper);
    matrix_scale(4, 1, point_lower, 1.0 / point_lower[3], point_lower);

    /* Project the points into the image */
    double proj_left[4], proj_right[4], proj_upper[4], proj_lower[4];

    matrix_product(3, 4, 4, 1, m_Pmatrix, point_left, proj_left);
    matrix_product(3, 4, 4, 1, m_Pmatrix, point_right, proj_right);
    matrix_product(3, 4, 4, 1, m_Pmatrix, point_upper, proj_upper);
    matrix_product(3, 4, 4, 1, m_Pmatrix, point_lower, proj_lower);

    matrix_scale(3, 1, proj_left, 1.0 / proj_left[2], proj_left);
    matrix_scale(3, 1, proj_right, 1.0 / proj_right[2], proj_right);
    matrix_scale(3, 1, proj_upper, 1.0 / proj_upper[2], proj_upper);
    matrix_scale(3, 1, proj_lower, 1.0 / proj_lower[2], proj_lower);

    /* Cross to get the horizon line */
    matrix_cross(proj_upper, proj_lower, m_horizon);

    /* Make sure the horizon has the correct orientation */
    double y_axis_img[3];
    matrix_product(3, 3, 3, 1, R, up, y_axis_img);

    double up_dir[2] = { y_axis_img[0], y_axis_img[1] };
    double norm = matrix_norm(2, 1, up_dir);
    matrix_scale(2, 1, up_dir, 1.0 / norm, up_dir);
    
    double up_dir3[3] = { up_dir[0], up_dir[1], 0.0 };
    double line_v[3];
    LineToUnitVector(m_horizon, line_v);

    double cross[3];
    matrix_cross(line_v, up_dir3, cross);
    
    if (cross[2] < 0.0) {
	matrix_scale(3, 1, m_horizon, -1.0, m_horizon);
    }
}

/* Returns true if the given 2D point is above the horizon */
bool CameraInfo::PointAboveHorizon(double *p) 
{
    double p3[3] = { p[0], p[1], 1.0 };
    double dot;
    matrix_product(1, 3, 3, 1, m_horizon, p3, &dot);
    
    return (dot > 0.0);
}

/* Returns true if the given point is in front of the camera */
bool CameraInfo::PointInFront(double *p) 
{
    /* Convert point to camera coordinates */
    double p_rot[3];
    matrix_product(3, 3, 3, 1, m_R, p, p_rot);
    matrix_sum(3, 1, 3, 1, p_rot, m_t, p_rot);
    
    /* Check z-coordinate */
    return (p_rot[2] < 0.0);
}

/* Interpolate between two camera views */
CameraInfo InterpolateCameras(const CameraInfo &cam1, 
                              const CameraInfo &cam2, double t)
{
    double pos1[3], pos2[3];

    cam1.GetPosition(pos1);
    cam2.GetPosition(pos2);

    double pos_new[3];
    pos_new[0] = (1.0 - t) * pos1[0] + t * pos2[0];
    pos_new[1] = (1.0 - t) * pos1[1] + t * pos2[1];
    pos_new[2] = (1.0 - t) * pos1[2] + t * pos2[2];

    /* Get quaternions */
    double R1[9], R2[9];
    cam1.GetPose(R1);
    cam2.GetPose(R2);

    double q1[4], q2[4];
    matrix_to_quaternion(R1, q1);
    matrix_to_quaternion(R2, q2);

    double dot;
    matrix_product(1, 3, 3, 1, q1, q2, &dot);

    if (dot < 0.0) {
        matrix_scale(4, 1, q2, -1.0, q2);
    }

    double q_new[4];
    q_new[0] = (1.0 - t) * q1[0] + t * q2[0];
    q_new[1] = (1.0 - t) * q1[1] + t * q2[1];
    q_new[2] = (1.0 - t) * q1[2] + t * q2[2];
    q_new[3] = (1.0 - t) * q1[3] + t * q2[3];

    /* Normalize */
    double mag = matrix_norm(4, 1, q_new);
    matrix_scale(4, 1, q_new, 1.0 / mag, q_new);

    double Rnew[9];
    quaternion_to_matrix(q_new, Rnew);

    CameraInfo cam_new;
    matrix_transpose(3, 3, Rnew, cam_new.m_R);

    matrix_product(3, 3, 3, 1, cam_new.m_R, pos_new, cam_new.m_t);
    matrix_scale(3, 1, cam_new.m_t, -1.0, cam_new.m_t);

    cam_new.m_focal = 1.0;

    return cam_new;
}

CameraInfo InterpolateCamerasThetaPhi(const CameraInfo &cam1, 
                                      const CameraInfo &cam2, double t,
                                      bool interp_fov)
{
    double pos1[3], pos2[3];

    cam1.GetPosition(pos1);
    cam2.GetPosition(pos2);

    double pos_new[3];
    pos_new[0] = (1.0 - t) * pos1[0] + t * pos2[0];
    pos_new[1] = (1.0 - t) * pos1[1] + t * pos2[1];
    pos_new[2] = (1.0 - t) * pos1[2] + t * pos2[2];

    /* Get quaternions */
    double v1[3], v2[3];
    cam1.GetViewDirection(v1);
    cam2.GetViewDirection(v2);

    double vt[3];
    slerp(v1, v2, t, vt);

    /* Scale by -1 to point in the right direction */
    matrix_scale(3, 1, vt, -1.0, vt);

    /* Project the y-axis onto the plane normal to vt */
    double dot = vt[1];
    double up[3] = { 0.0 - dot * vt[0],
                     1.0 - dot * vt[1],
                     0.0 - dot * vt[2] };

    double norm = matrix_norm(3, 1, up);
    matrix_scale(3, 1, up, 1.0 / norm, up);

    double xaxis[3];
    matrix_cross(up, vt, xaxis);

    double Rnew[9];
    memcpy(Rnew + 0, xaxis, 3 * sizeof(double));
    memcpy(Rnew + 3, up, 3 * sizeof(double));
    memcpy(Rnew + 6, vt, 3 * sizeof(double));

    CameraInfo cam_new;
    // matrix_transpose(3, 3, Rnew, cam_new.m_R);
    memcpy(cam_new.m_R, Rnew, 9 * sizeof(double));

    matrix_product(3, 3, 3, 1, cam_new.m_R, pos_new, cam_new.m_t);
    matrix_scale(3, 1, cam_new.m_t, -1.0, cam_new.m_t);

    if (!interp_fov) {
        cam_new.m_width = 1024;
        cam_new.m_height = 768;
        cam_new.SetFOV(70.0); // 60.0
    } else {
        cam_new.m_width = cam1.m_width;
        cam_new.m_height = cam1.m_height;
        double fov = (1.0 - t) * cam1.GetFOV() + t * cam2.GetFOV();
        cam_new.SetFOV(RAD2DEG(fov));
    }

    return cam_new;
}

#ifndef __BUNDLER__
/* Return a bezier curve given another camera */
Bezier CameraInfo::ComputeBezier(const CameraInfo &cam) const {
    double pos1[3], pos2[3];
    double R1[9], R2[9];

    GetPosition(pos1);
    cam.GetPosition(pos2);
    
    GetPose(R1);
    cam.GetPose(R2);
    
    /* Compute the x-axis direction in the camera coordinate frames */
    double xaxis[3] = { 1, 0, 0 };
    double yaxis[3] = { 0, 1, 0 };
    double zaxis[3] = { 0, 0, -1 };

    double x1[3], x2[3];
    double y1[3], y2[3];
    double z1[3], z2[3];

#if 1
    matrix_product(3, 3, 3, 1, R1, xaxis, x1);
    matrix_product(3, 3, 3, 1, R2, xaxis, x2);

    matrix_product(3, 3, 3, 1, R1, yaxis, y1);
    matrix_product(3, 3, 3, 1, R2, yaxis, y2);

    matrix_product(3, 3, 3, 1, R1, zaxis, z1);
    matrix_product(3, 3, 3, 1, R2, zaxis, z2);
#else
    matrix_transpose_product(3, 3, 3, 1, R1, xaxis, x1);
    matrix_transpose_product(3, 3, 3, 1, R2, xaxis, x2);

    matrix_transpose_product(3, 3, 3, 1, R1, yaxis, y1);
    matrix_transpose_product(3, 3, 3, 1, R2, yaxis, y2);

    matrix_transpose_product(3, 3, 3, 1, R1, zaxis, z1);
    matrix_transpose_product(3, 3, 3, 1, R2, zaxis, z2);
#endif

    double disp[3], dist;
    matrix_diff(3, 1, 3, 1, pos2, pos1, disp);
    dist = matrix_norm(3, 1, disp);

    /* Check the direction of the xaxes */
    double x_dot, y_dot, z_dot;
    matrix_product(1, 3, 3, 1, disp, x1, &x_dot);
    matrix_product(1, 3, 3, 1, disp, y1, &y_dot);
    matrix_product(1, 3, 3, 1, disp, z1, &z_dot);
    
    double dir1[3], dir2[3];
    if (fabs(x_dot) >= fabs(y_dot) && fabs(x_dot) >= fabs(z_dot)) {
	if (x_dot < 0)
	    matrix_scale(3, 1, x1, -1.0, dir1);
	else
	    matrix_scale(3, 1, x1, 1.0, dir1);
    } else if (fabs(y_dot) >= fabs(x_dot) && fabs(y_dot) >= fabs(z_dot)) {
	if (y_dot < 0)
	    matrix_scale(3, 1, y1, -1.0, dir1);
	else
	    matrix_scale(3, 1, y1, 1.0, dir1);
    } else if (fabs(z_dot) >= fabs(x_dot) && fabs(z_dot) >= fabs(y_dot)) {
	if (z_dot < 0)
	    matrix_scale(3, 1, z1, -1.0, dir1);
	else
	    matrix_scale(3, 1, z1, 1.0, dir1);
    }

#if 0    
    if (dot < 0)
	matrix_scale(3, 1, x1, -1.0, x1);
#endif    

    matrix_product(1, 3, 3, 1, disp, x2, &x_dot);
    matrix_product(1, 3, 3, 1, disp, y2, &y_dot);
    matrix_product(1, 3, 3, 1, disp, z2, &z_dot);

    if (fabs(x_dot) >= fabs(y_dot) && fabs(x_dot) >= fabs(z_dot)) {
	if (x_dot > 0)
	    matrix_scale(3, 1, x2, -1.0, dir2);
	else
	    matrix_scale(3, 1, x2, 1.0, dir2);
    } else if (fabs(y_dot) >= fabs(x_dot) && fabs(y_dot) >= fabs(z_dot)) {
	if (y_dot > 0)
	    matrix_scale(3, 1, y2, -1.0, dir2);
	else
	    matrix_scale(3, 1, y2, 1.0, dir2);
    } else if (fabs(z_dot) >= fabs(x_dot) && fabs(z_dot) >= fabs(y_dot)) {
	if (z_dot > 0)
	    matrix_scale(3, 1, z2, -1.0, dir2);
	else
	    matrix_scale(3, 1, z2, 1.0, dir2);
    }

#if 0
    if (dot > 0)
	matrix_scale(3, 1, x2, -1.0, x2);
#endif

    double disp1[3], disp2[3];
    double mid1[3], mid2[3];
    matrix_scale(3, 1, dir1, 0.1 * dist, disp1);
    matrix_scale(3, 1, dir2, 0.1 * dist, disp2);

    matrix_sum(3, 1, 3, 1, pos1, disp1, mid1);
    matrix_sum(3, 1, 3, 1, pos2, disp2, mid2);

    return Bezier(pos1, mid1, mid2, pos2);
}
#endif

/* Convert a pixel position to a ray direction */
void CameraInfo::PixelToCameraRay(double x, double y, double *ray) 
{
    double ray_cam[3] = { x, y, -m_focal };

#if 0
    matrix_transpose_product(3, 3, 3, 1, m_R, ray_cam, ray);
#endif

    double norm = matrix_norm(3, 1, ray_cam);
    matrix_scale(3, 1, ray_cam, 1.0 / norm, ray);
}

/* Convert a pixel position to an (absolute) ray direction */
void CameraInfo::PixelToCameraRayAbsolute(double x, double y, double *ray) 
{
    double ray_cam[3] = { x, y, -m_focal };

    matrix_transpose_product(3, 3, 3, 1, m_R, ray_cam, ray);

    double norm = matrix_norm(3, 1, ray);
    matrix_scale(3, 1, ray, 1.0 / norm, ray);
}

static void SphToRot(double theta, double phi, double *R)
{
    double v[3];
    
    v[0] = -cos(theta) * sin(phi);
    v[1] = -cos(phi);
    v[2] = -sin(theta) * sin(phi);

    /* Compute the new up vector */
    double phi_up = phi - 0.5 * M_PI;
    double theta_up = theta;
    double up[3];
    
    up[0] = cos(theta_up) * sin(phi_up);
    up[1] = cos(phi_up);
    up[2] = sin(theta_up) * sin(phi_up);

    double x_axis[3];
    matrix_cross(up, v, x_axis);
    
    memcpy(R + 0, x_axis, 3 * sizeof(double));
    memcpy(R + 3, up, 3 * sizeof(double));
    memcpy(R + 6, v, 3 * sizeof(double));
}

/* Point the camera in a different direction */
void CameraInfo::PointAt(double *ray)
{
    double r = matrix_norm(3, 1, ray);

    /* Convert the ray to spherical coordinates */
    double theta = atan2(ray[2], ray[0]);
    double phi = acos(ray[1] / r);

    double R[9];
    SphToRot(theta, phi, R);

    double pos[3];
    GetPosition(pos);
    
    double Rnew[9];
    matrix_product(3, 3, 3, 3, R, m_R, Rnew);

    memcpy(m_R, Rnew, 9 * sizeof(double));
    SetPosition(pos);
    
    Finalize();
}

/* Point the camera in a different direction */
void CameraInfo::PointAtAbsolute(double *ray)
{
    double r = matrix_norm(3, 1, ray);

    /* Convert the ray to spherical coordinates */
    double theta = atan2(ray[2], ray[0]);
    double phi = acos(ray[1] / r);

    double R[9];
    SphToRot(theta, phi, R);

    double pos[3];
    GetPosition(pos);
    
    memcpy(m_R, R, 9 * sizeof(double));
    SetPosition(pos);
    
    Finalize();
}

/* Get the bounding box for this camera */
BoundingBox CameraInfo::GetBoundingBox() const
{
    int w = m_width, h = m_height;

    return BoundingBox(-0.5 * w, -0.5 * h, 0.5 * w, 0.5 * h);
}

/* Return the view direction */
void CameraInfo::GetViewDirection(double *view) const
{
    // double R[9];
    // GetPose(R);
    // double minuz[3] = { 0.0, 0.0, -1.0 };
    // matrix_product(3, 3, 3, 1, R, minuz, view);

    view[0] = -m_R[6];
    view[1] = -m_R[7];
    view[2] = -m_R[8];
}

/* Get the twist angle of the camera */
double CameraInfo::GetTwistAngleRadians() const
{
    double R[9];
    GetPose(R);

    double c_twist = 
        (R[0] * R[8] - R[6] * R[2]) / sqrt(1 - R[5] * R[5]);
    
    c_twist = CLAMP(c_twist, -1.0 + 1.0e-8, 1.0 - 1.0e-8);

    double angle = acos(c_twist);  

    if (R[3] < 0.0)
        return -angle;
    else 
        return angle;
}

/* Return the halfspace in front of the camera */
void CameraInfo::GetFrontHalfspace(double *plane) const
{
    double v[3];
    GetViewDirection(v);

    double pos[3];
    GetPosition(pos);

    /* Set the position slightly in front of the real position */
    pos[0] += 1.0e-6 * v[0];
    pos[1] += 1.0e-6 * v[1];
    pos[2] += 1.0e-6 * v[2];

    double dot;
    matrix_product(1, 3, 3, 1, v, pos, &dot);

    plane[0] = v[0];
    plane[1] = v[1];
    plane[2] = v[2];
    plane[3] = -dot;
}

bool CameraInfo::PointInsideImage(const double *p) const
{
	double proj[2];
	bool in_front = Project(p, proj);
	if (!in_front) return false;
	return (proj[0]>-0.5*m_width && proj[0]<0.5*m_width && proj[1]<0.5*m_height && proj[1]>-0.5*m_height);	
}

/* Get a camera whose up vector points in the right direction */
CameraInfo CameraInfo::GetUpCamera(double *up) const 
{
    double pos[3];
    GetPosition(pos);

    double up_image[3];
    matrix_product(3, 3, 3, 1, (double *) m_R, up, up_image);
    
    /* Project up-image onto the xy-plane */
    double up_image_proj[3] = { up_image[0], up_image[1], 0.0 };
    double mag = matrix_norm(3, 1, up_image_proj);
    matrix_scale(3, 1, up_image_proj, 1.0 / mag, up_image_proj);

    double yaxis[3] = { 0.0, 1.0, 0.0 };
    double dot;
    matrix_product(1, 3, 3, 1, up_image_proj, yaxis, &dot);

    double angle = acos(dot);

    double cross[3];
    matrix_cross(up_image_proj, yaxis, cross);

    mag = matrix_norm(3, 1, cross);
    matrix_scale(3, 1, cross, 1.0 / mag, cross);

    /* Rotate '-angle' radians around the cross axis */
    double Rroll[9];
    axis_angle_to_matrix(cross, -angle, Rroll);

    /* The world to camera transformation becomes Rroll * R */
#if 0
    double Rz[9] = { cos(angle), -sin(angle), 0.0,
		     sin(angle),  cos(angle), 0.0,
		     0.0,         0.0,        1.0   };

    double Rnew[9];
    matrix_product(3, 3, 3, 3, Rz, (double *) m_R, Rnew);
#endif

    double Rnew[9];
    matrix_transpose_product(3, 3, 3, 3, Rroll, (double *) m_R, Rnew);
    
    double RnewT[9];
    matrix_transpose(3, 3, Rnew, RnewT);

    CameraInfo up_camera = *this;
    up_camera.SetPose(RnewT);
    up_camera.SetPosition(pos);
    up_camera.Finalize();

    matrix_product(3, 3, 3, 1, Rnew, up, up_image_proj);
    if (up_image_proj[0] > 1.0e-4) 
	printf("Error in up camera computation\n");

    return up_camera;
}

#ifndef __BUNDLER__
/* Read/write the links to a file */
void CameraInfo::ReadLinks(FILE *f)
{
    fscanf(f, "%d %d %d %d %d %d %d %d\n",
	   m_links + DIRECTION_LEFT,
	   m_links + DIRECTION_RIGHT,
	   m_links + DIRECTION_FORWARD,
	   m_links + DIRECTION_BACKWARD,
	   m_links + DIRECTION_UP,
	   m_links + DIRECTION_DOWN,
	   m_links + DIRECTION_ZOOM_OUT,
	   m_links + DIRECTION_ZOOM_IN);    
}

void CameraInfo::WriteLinks(FILE *f)
{
    fprintf(f, "%d %d %d %d %d %d %d %d\n",
	    m_links[DIRECTION_LEFT],
	    m_links[DIRECTION_RIGHT],
	    m_links[DIRECTION_FORWARD],
	    m_links[DIRECTION_BACKWARD],
	    m_links[DIRECTION_UP],
	    m_links[DIRECTION_DOWN],
	    m_links[DIRECTION_ZOOM_OUT],
	    m_links[DIRECTION_ZOOM_IN]);
}
#endif /* __BUNDLER__ */

void CameraInfo::WriteFile(FILE *f)
{
    if (m_adjusted) {
        double temp[3];
        GetPosition(temp);
        fprintf(f,"%f %f %f\n",temp[0],temp[1],temp[2]);
    } else {
        fprintf(f, "0.0 0.0 0.0\n");
    }
}

void CameraInfo::WriteXML(FILE *f)
{
    static const char *spacer = "    ";

    fprintf(f, "%s<focal> %0.8e </focal>\n", spacer, m_focal);
    fprintf(f, "%s<rot> %0.8e %0.8e %0.8e "
	               "%0.8e %0.8e %0.8e "
	               "%0.8e %0.8e %0.8e </rot>\n", 
	    spacer, 
	    m_R[0], m_R[1], m_R[2], 
	    m_R[3], m_R[4], m_R[5],
	    m_R[6], m_R[7], m_R[8]);

    fprintf(f, "%s<t> %0.8e %0.8e %0.8e </t>\n", spacer, 
	    m_t[0], m_t[1], m_t[2]);    

    // fprintf(f, "%s<back> %d </back>\n", spacer, m_links[DIRECTION_ZOOM_OUT]);
}


