/* 
 *  Copyright (c) 2008  Noah Snavely (snavely (at) cs.washington.edu)
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

/* sfm-driver.c */
/* Driver for sfm routines */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sba.h"

#include "matrix.h"
#include "vector.h"
#include "sfm.h"

#ifdef WIN32
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif /* WIN32 */

// #define COLIN_HACK
#define TEST_FOCAL

typedef struct {
    int num_cameras;               /* Number of cameras */
    int num_points;                /* Number of points */
    int num_params_per_camera;     /* Number of parameters for each camera */

    int est_focal_length;          /* Should the focal length be estimated? */
    int const_focal_length;        /* Is the focal length constant for all
			            * cameras? */
    int explicit_camera_centers;   /* Are the camera centers explicit? */
    int estimate_distortion;       /* Apply undistortion? */
    
    camera_params_t global_params;
    camera_params_t *init_params;  /* Initial camera parameters */

    v3_t *points;
} sfm_global_t;

static void *safe_malloc(int n, char *where)
{
    void *mem = malloc(n);
    
    if (mem == NULL) {
	if (where) {
	    printf("[safe_malloc] Error allocating %d bytes "
		   "of memory at %s\n", n, where);
	} else {
	    printf("[safe_malloc] Error allocating %d bytes of memory\n", n);
	}

	fflush(stdout);
	exit(1);
    }

    return mem;
}

/* Compute an updated rotation matrix given the initial rotation (R)
 * and the correction (w) */
void rot_update(double *R, double *w, double *Rnew) 
{
    double theta, sinth, costh, n[3];
    double nx[9], nxsq[9];
    double term2[9], term3[9];
    double tmp[9], dR[9];

    double ident[9] = 
	{ 1.0, 0.0, 0.0,
	  0.0, 1.0, 0.0,
	  0.0, 0.0, 1.0 };

    theta = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);

    if (theta == 0.0) {
	memcpy(Rnew, R, sizeof(double) * 9);
	return;
    }

    n[0] = w[0] / theta;
    n[1] = w[1] / theta;
    n[2] = w[2] / theta;

    nx[0] = 0.0;   nx[1] = -n[2];  nx[2] = n[1];
    nx[3] = n[2];  nx[4] = 0.0;    nx[5] = -n[0];
    nx[6] = -n[1]; nx[7] = n[0];   nx[8] = 0.0;

    matrix_product33(nx, nx, nxsq);

    sinth = sin(theta);
    costh = cos(theta);

    matrix_scale(3, 3, nx, sinth, term2);
    matrix_scale(3, 3, nxsq, 1.0 - costh, term3);

    matrix_sum(3, 3, 3, 3, ident, term2, tmp);
    matrix_sum(3, 3, 3, 3, tmp, term3, dR);

    matrix_product33(dR, R, Rnew);
}

v2_t sfm_project_final(camera_params_t *params, v3_t pt, 
		       int explicit_camera_centers, int undistort)
{
    double b_cam[3], b_proj[3];
    double b[3] = { Vx(pt), Vy(pt), Vz(pt) };
    v2_t proj;

    /* Project! */
    if (!explicit_camera_centers) {
	matrix_product331(params->R, b, b_cam);
	b_cam[0] += params->t[0];
	b_cam[1] += params->t[1];
	b_cam[2] += params->t[2];
    } else {
	double b2[3];
	b2[0] = b[0] - params->t[0];
	b2[1] = b[1] - params->t[1];
	b2[2] = b[2] - params->t[2];
	matrix_product331(params->R, b2, b_cam);
    }
    
    if (!params->known_intrinsics) {
        b_proj[0] = b_cam[0] * params->f;
        b_proj[1] = b_cam[1] * params->f;
        b_proj[2] = b_cam[2];

        b_proj[0] /= -b_proj[2];
        b_proj[1] /= -b_proj[2];

        if (undistort) {
            double rsq = 
                (b_proj[0] * b_proj[0] + b_proj[1] * b_proj[1]) / 
                (params->f * params->f);
            double factor = 1.0 + params->k[0] * rsq + params->k[1] * rsq * rsq;

            b_proj[0] *= factor;
            b_proj[1] *= factor;
        }
    } else {
        /* Apply intrinsics */
        double x_n = -b_cam[0] / b_cam[2];
        double y_n = -b_cam[1] / b_cam[2];

        double *k = params->k_known;
	double rsq = x_n * x_n + y_n * y_n;
	double factor = 1.0 + k[0] * rsq + 
            k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

        double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
        double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

	double x_d = x_n * factor + dx_x;
	double y_d = y_n * factor + dx_y;

        double *K = params->K_known;
        b_proj[0] = K[0] * x_d + K[1] * y_d + K[2];
        b_proj[1] = K[4] * y_d + K[5];        
    }

    Vx(proj) = b_proj[0];
    Vy(proj) = b_proj[1];

    return proj;
}

void sfm_project(camera_params_t *init, double *K, 
		 double *w, double *dt, double *b, double *p,
		 int explicit_camera_centers)
{
    double *R, *t;

    double Rnew[9];
    double tnew[3];

    double b_cam[3], b_proj[3];

    R = init->R;
    t = init->t;

    rot_update(R, w, Rnew);

    tnew[0] = dt[0]; // t[0] + dt[0];
    tnew[1] = dt[1]; // t[1] + dt[1];  // 0.0
    tnew[2] = dt[2]; // t[2] + dt[2];  // 0.0

    /* Project! */
    if (!explicit_camera_centers) {
	matrix_product331(Rnew, b, b_cam);
	b_cam[0] += tnew[0];
	b_cam[1] += tnew[1];
	b_cam[2] += tnew[2];
    } else {
	double b2[3];
	b2[0] = b[0] - tnew[0];
	b2[1] = b[1] - tnew[1];
	b2[2] = b[2] - tnew[2];
	matrix_product331(Rnew, b2, b_cam);
    }
    
    if (!init->known_intrinsics) {
        matrix_product331(K, b_cam, b_proj);
        p[0] = -b_proj[0] / b_proj[2];
        p[1] = -b_proj[1] / b_proj[2];
    } else {
        /* Apply intrinsics */
        double x_n = -b_cam[0] / b_cam[2];
        double y_n = -b_cam[1] / b_cam[2];

        double *k = init->k_known;
	double rsq = x_n * x_n + y_n * y_n;
	double factor = 1.0 + k[0] * rsq + 
            k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

        double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
        double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

	double x_d = x_n * factor + dx_x;
	double y_d = y_n * factor + dx_y;

        double *K = init->K_known;
        p[0] = K[0] * x_d + K[1] * y_d + K[2];
        p[1] = K[4] * y_d + K[5];
    }
}

void sfm_project2(camera_params_t *init, double f,
                  double *R, double *dt, double *b, double *p,
                  int explicit_camera_centers)
{
    double *t;

    double tnew[3];
    double b_cam[3];

    t = init->t;

    tnew[0] = dt[0];
    tnew[1] = dt[1];
    tnew[2] = dt[2];

    /* Project! */
    if (!explicit_camera_centers) {
	matrix_product331(R, b, b_cam);
	b_cam[0] += tnew[0];
	b_cam[1] += tnew[1];
	b_cam[2] += tnew[2];
    } else {
	double b2[3];
	b2[0] = b[0] - tnew[0];
	b2[1] = b[1] - tnew[1];
	b2[2] = b[2] - tnew[2];

	// matrix_product(3, 3, 3, 1, R, b2, b_cam);
        b_cam[0] = R[0] * b2[0] + R[1] * b2[1] + R[2] * b2[2];
        b_cam[1] = R[3] * b2[0] + R[4] * b2[1] + R[5] * b2[2];
        b_cam[2] = R[6] * b2[0] + R[7] * b2[1] + R[8] * b2[2];
    }
    
    if (!init->known_intrinsics) {
        p[0] = -b_cam[0] * f / b_cam[2];
        p[1] = -b_cam[1] * f / b_cam[2];
    } else {
        /* Apply intrinsics */
        double x_n = -b_cam[0] / b_cam[2];
        double y_n = -b_cam[1] / b_cam[2];

        double *k = init->k_known;
	double rsq = x_n * x_n + y_n * y_n;
	double factor = 1.0 + k[0] * rsq + 
            k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

        double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
        double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

	double x_d = x_n * factor + dx_x;
	double y_d = y_n * factor + dx_y;

        double *K = init->K_known;
        p[0] = K[0] * x_d + K[1] * y_d + K[2];
        p[1] = K[4] * y_d + K[5];
    }
}


void sfm_project_rd(camera_params_t *init, double *K, double *k,
                    double *R, double *dt, double *b, double *p,
                    int undistort, int explicit_camera_centers)
{
    double *t;

    double tnew[3];
    double b_cam[3];

    t = init->t;

    tnew[0] = dt[0];
    tnew[1] = dt[1];
    tnew[2] = dt[2];

    /* Project! */
    if (!explicit_camera_centers) {
	matrix_product331(R, b, b_cam);
	b_cam[0] += tnew[0];
	b_cam[1] += tnew[1];
	b_cam[2] += tnew[2];
    } else {
	double b2[3];
	b2[0] = b[0] - tnew[0];
	b2[1] = b[1] - tnew[1];
	b2[2] = b[2] - tnew[2];
	matrix_product331(R, b2, b_cam);
    }
    
    if (!init->known_intrinsics) {
        p[0] = -b_cam[0] * K[0] / b_cam[2];
        p[1] = -b_cam[1] * K[0] / b_cam[2];
    } else {
        /* Apply intrinsics */
        double x_n = -b_cam[0] / b_cam[2];
        double y_n = -b_cam[1] / b_cam[2];

        double *k = init->k_known;
	double rsq = x_n * x_n + y_n * y_n;
	double factor = 1.0 + k[0] * rsq + 
            k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

        double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
        double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

	double x_d = x_n * factor + dx_x;
	double y_d = y_n * factor + dx_y;

        double *K = init->K_known;
        p[0] = K[0] * x_d + K[1] * y_d + K[2];
        p[1] = K[4] * y_d + K[5];
    }

    // p[0] = b_cam[0] * K[0] / b_cam[2];
    // p[1] = b_cam[1] * K[0] / b_cam[2];

    /* Apply radial distortion */
    if (undistort) {
#ifndef TEST_FOCAL
        double k1 = k[0], k2 = k[1];
#else
        double k1 = k[0] / init->k_scale;
        double k2 = k[1] / init->k_scale;
#endif

	double rsq = (p[0] * p[0] + p[1] * p[1]) / (K[0] * K[0]);
	double factor = 1.0 + k1 * rsq + k2 * rsq * rsq;

	p[0] *= factor;
	p[1] *= factor;
    }
}

static double *global_last_ws = NULL;
static double *global_last_Rs = NULL;

static void sfm_project_point(int j, int i, double *aj, double *bi, 
                              double *xij, void *adata)
{
    sfm_global_t *globs = (sfm_global_t *) adata;

    double K[9] = 
	{ 1.0, 0.0, 0.0,
	  0.0, 1.0, 0.0,
	  0.0, 0.0, 1.0 };

    double *w, *dt;

    /* Compute intrinsics */
    if (!globs->est_focal_length) {
	K[0] = K[4] = globs->init_params[j].f; // globs->global_params.f;
    } else if (globs->const_focal_length) {
	printf("Error: case of constant focal length "
	       "has not been implemented.\n");
	K[0] = K[4] = globs->global_params.f;
    } else {
	K[0] = K[4] = aj[6];
    }
    
    /* Compute translation, rotation update */
    dt = aj + 0;
    w = aj + 3;

#ifdef COLIN_HACK
    w[0] = w[1] = w[2] = 0.0;
    dt[2] = 0.0;
#endif

    if (globs->estimate_distortion == 0) {
        sfm_project(globs->init_params + j, K, w, dt, bi, xij, 
                    globs->explicit_camera_centers);
    } else {
        double Rnew[9];
        rot_update(globs->init_params->R, w, Rnew);
        sfm_project_rd(globs->init_params + j, K, aj + 7,
                       Rnew, dt, bi, xij, 1, globs->explicit_camera_centers);
    }
}

// k_scale 100.0
// focal_scale 0.001

static void sfm_fisheye_distort(camera_params_t *cam,
                                double *x_u, double *x_d)
{
    double xn, yn, r, angle, rnew;

    if (cam->fisheye == 0) {
        x_d[0] = x_u[0];
        x_d[1] = x_u[1];
        return;
    }
    
    xn = x_u[0];
    yn = x_u[1];
    
    r = sqrt(xn * xn + yn * yn);
    angle = 180.0 * atan(r / cam->f_focal) / M_PI;
    rnew = cam->f_rad * angle / (0.5 * cam->f_angle);
    
    x_d[0] = xn * (rnew / r) + cam->f_cx;
    x_d[1] = yn * (rnew / r) + cam->f_cy;    
}

static void sfm_project_point2_fisheye(int j, int i, double *aj, double *bi, 
                                       double *xij, void *adata)
{
    sfm_global_t *globs = (sfm_global_t *) adata;

    double f;

    double *w, *dt;
    double xij_tmp[2];

    /* Compute intrinsics */
    if (!globs->est_focal_length) {
	f = globs->init_params[j].f; // globs->global_params.f;
    } else if (globs->const_focal_length) {
	printf("Error: case of constant focal length "
	       "has not been implemented.\n");
	f = globs->global_params.f;
    } else {
#ifndef TEST_FOCAL
	f = aj[6];
#else
        f = aj[6] / globs->init_params[j].f_scale;
#endif
    }
    
    /* Compute translation, rotation update */
    dt = aj + 0;
    w = aj + 3;

#ifdef COLIN_HACK
    w[0] = w[1] = w[2] = 0.0;
    dt[2] = 0.0;
#endif

    if (w[0] != global_last_ws[3 * j + 0] ||
	w[1] != global_last_ws[3 * j + 1] ||
	w[2] != global_last_ws[3 * j + 2]) {

	// printf("updating w: %0.3f, %0.3f, %0.3f\n", w[0], w[1], w[2]);

	rot_update(globs->init_params[j].R, w, global_last_Rs + 9 * j);
	global_last_ws[3 * j + 0] = w[0];
	global_last_ws[3 * j + 1] = w[1];
	global_last_ws[3 * j + 2] = w[2];
    }
    
    sfm_project2(globs->init_params + j, f, global_last_Rs + 9 * j, 
		 dt, bi, xij_tmp, globs->explicit_camera_centers);

    /* Distort the point */
    sfm_fisheye_distort(globs->init_params + j, xij_tmp, xij);
}

static void sfm_project_point2_fisheye_mot(int j, int i, double *aj, 
                                           double *xij, void *adata)
{
    sfm_global_t *globs = (sfm_global_t *) adata;
    double *b = globs->points[i].p;

    sfm_project_point2_fisheye(j, i, aj, b, xij, adata);
}

static void sfm_project_point3(int j, int i, double *aj, double *bi, 
			       double *xij, void *adata)
{
    sfm_global_t *globs = (sfm_global_t *) adata;

    double K[9] = 
	{ 1.0, 0.0, 0.0,
	  0.0, 1.0, 0.0,
	  0.0, 0.0, 1.0 };

    double *w, *dt, *k;

    /* Compute intrinsics */
    if (!globs->est_focal_length) {
	K[0] = K[4] = globs->init_params[j].f; // globs->global_params.f;
    } else if (globs->const_focal_length) {
	printf("Error: case of constant focal length "
	       "has not been implemented.\n");
	K[0] = K[4] = globs->global_params.f;
    } else {
#ifndef TEST_FOCAL
	K[0] = K[4] = aj[6];
#else
	K[0] = K[4] = aj[6] / globs->init_params[j].f_scale;
#endif
    }
    
    /* Compute translation, rotation update */
    dt = aj + 0;
    w = aj + 3;

    if (globs->est_focal_length)
        k = aj + 7;
    else
        k = aj + 6;

    if (w[0] != global_last_ws[3 * j + 0] ||
	w[1] != global_last_ws[3 * j + 1] ||
	w[2] != global_last_ws[3 * j + 2]) {

	rot_update(globs->init_params[j].R, w, global_last_Rs + 9 * j);
	global_last_ws[3 * j + 0] = w[0];
	global_last_ws[3 * j + 1] = w[1];
	global_last_ws[3 * j + 2] = w[2];
    }
    
    sfm_project_rd(globs->init_params + j, K, k, global_last_Rs + 9 * j, 
                   dt, bi, xij, globs->estimate_distortion, 
                   globs->explicit_camera_centers);
}

static void sfm_project_point3_mot(int j, int i, double *aj, 
                                   double *xij, void *adata)
{
    sfm_global_t *globs = (sfm_global_t *) adata;
    double *b = globs->points[i].p;

    sfm_project_point3(j, i, aj, b, xij, adata);    
}

static void sfm_mot_project_point(int j, int i, double *bi, 
				  double *xij, void *adata)
{
    sfm_global_t *globs = (sfm_global_t *) adata;

    double K[9] = 
	{ 1.0, 0.0, 0.0,
	  0.0, 1.0, 0.0,
	  0.0, 0.0, 1.0 };

    double w[3] = { 0.0, 0.0, 0.0 }, dt[3] = { 0.0, 0.0, 0.0 };

    /* Compute intrinsics */
    K[0] = K[4] = globs->init_params[j].f; // globs->global_params.f;
    
    if (1 || j > 0) {
	sfm_project(globs->init_params + j, K, w, dt, bi, xij, 
		    globs->explicit_camera_centers);
    } else {
	/* Keep first camera at origin */
	double w0[3] = { 0.0, 0.0, 0.0 };
	double dt0[3] = { 0.0, 0.0, 0.0 };
	sfm_project(globs->init_params + j, K, w0, dt0, bi, xij,
		    globs->explicit_camera_centers);
    }
}

#define SBA_V121

void run_sfm(int num_pts, int num_cameras, int ncons,
             char *vmask,
             double *projections,
             int est_focal_length,
             int const_focal_length,
             int undistort,
             int explicit_camera_centers,
             camera_params_t *init_camera_params,
             v3_t *init_pts, 
             int use_constraints, 
             int use_point_constraints,
             v3_t *pt_constraints,
             double pt_constraint_weight,
             int fix_points,
             int optimize_for_fisheye,
             double eps2,
             double *Vout, 
             double *Sout,
             double *Uout, double *Wout
             /* size num_cameras ** 2 * cnp * cnp */)
{
    int cnp;
    double *params;

#ifdef SBA_V121
    double opts[6]; // opts[5];
#else
    double opts[3];
#endif
    double info[10];

    int i, j, idx, base;
    int num_camera_params, num_pt_params, num_params;

    sfm_global_t global_params;

    // #ifndef SBA_V121
    camera_constraints_t *constraints = NULL;
    // #endif

    point_constraints_t *point_constraints = NULL;

    const double f_scale = 0.001;
    const double k_scale = 5.0;

    if (est_focal_length)
	cnp = 7;
    else
	cnp = 6;

    if (undistort)
        cnp += 2;

    num_camera_params = cnp * num_cameras;
    num_pt_params = 3 * num_pts;
    num_params = num_camera_params + num_pt_params;

    params = (double *) safe_malloc(sizeof(double) * num_params, "params");
    
    /* Fill parameters */
    for (j = 0; j < num_cameras; j++) {
        int c = 0;

#ifdef TEST_FOCAL
        init_camera_params[j].f_scale = f_scale;
        init_camera_params[j].k_scale = k_scale;
#else
        init_camera_params[j].f_scale = 1.0;
        init_camera_params[j].k_scale = 1.0;
#endif

	/* Translation is zero */
	params[cnp * j + 0] = init_camera_params[j].t[0]; // 0.0;
	params[cnp * j + 1] = init_camera_params[j].t[1]; // 0.0;
	params[cnp * j + 2] = init_camera_params[j].t[2]; // 0.0;
	
	/* Rotation is zero */
	params[cnp * j + 3] = 0.0;
	params[cnp * j + 4] = 0.0;
	params[cnp * j + 5] = 0.0;

	if (est_focal_length) {
	    /* Focal length is initial estimate */
#ifndef TEST_FOCAL
	    params[cnp * j + 6] = init_camera_params[j].f;
#else
	    params[cnp * j + 6] = 
                init_camera_params[j].f * init_camera_params[j].f_scale;
#endif
            c = 7;
	} else {
            c = 6;
        }
        
        if (undistort) {
#ifndef TEST_FOCAL
            params[cnp * j + c] = init_camera_params[j].k[0];
            params[cnp * j + c+1] = init_camera_params[j].k[1];
#else
            double scale = init_camera_params[j].k_scale;
            params[cnp * j + c] = init_camera_params[j].k[0] * scale;
            params[cnp * j + c+1] = init_camera_params[j].k[1] * scale;
#endif
        }
    }

    base = num_camera_params;
    for (i = 0; i < num_pts; i++) {
	params[base + 3 * i + 0] = Vx(init_pts[i]);
	params[base + 3 * i + 1] = Vy(init_pts[i]);
	params[base + 3 * i + 2] = Vz(init_pts[i]);
    }

    opts[0] = 1.0e-3;
    opts[1] = 1.0e-10; // 1.0e-15;
    opts[2] = eps2; // 0.0;  // 1.0e-10; // 1.0e-15;

#ifdef SBA_V121
    opts[3] = 1.0e-12;
    // opts[4] = 4.0e-2;
    opts[4] = 0.0;
    opts[5] = 4.0e-2; // change this back to opts[4] for sba v1.2.1
#endif

    // opts[1] = 1.0e-8;
    // opts[2] = 1.0e-8;

    // #ifndef SBA_V121
    /* Create the constraints */
    if (use_constraints) {
	constraints = 
	    (camera_constraints_t *) 
	       malloc(num_cameras * sizeof(camera_constraints_t));	

	for (i = 0; i < num_cameras; i++) {
	    constraints[i].constrained = (char *) malloc(cnp);
	    constraints[i].constraints = 
		(double *) malloc(sizeof(double) * cnp);
	    constraints[i].weights = (double *) malloc(sizeof(double) * cnp);

	    memcpy(constraints[i].constrained, 
		   init_camera_params[i].constrained, cnp * sizeof(char));
	    memcpy(constraints[i].constraints, 
		   init_camera_params[i].constraints, cnp * sizeof(double));
	    memcpy(constraints[i].weights,
		   init_camera_params[i].weights, cnp * sizeof(double));

#ifdef TEST_FOCAL
            if (est_focal_length) {
                constraints[i].constraints[6] *= f_scale;
                constraints[i].weights[6] *= (1.0 / (f_scale * f_scale));
            }
            
            if (undistort) {
                constraints[i].constraints[7] *= k_scale;
                constraints[i].weights[7] *= (1.0 / (k_scale * k_scale));

                constraints[i].constraints[8] *= k_scale;
                constraints[i].weights[8] *= (1.0 / (k_scale * k_scale));
            }
#endif
	}
    }
    // #endif

    if (use_point_constraints) {
	point_constraints = 
	    (point_constraints_t *) 
	        malloc(num_pts * sizeof(point_constraints_t));
	
	for (i = 0; i < num_pts; i++) {
	    if (Vx(pt_constraints[i]) == 0.0 &&
		Vy(pt_constraints[i]) == 0.0 &&
		Vz(pt_constraints[i]) == 0.0) {
		
		point_constraints[i].constrained = 0;
		point_constraints[i].constraints[0] = 0.0;
		point_constraints[i].constraints[1] = 0.0;
		point_constraints[i].constraints[2] = 0.0;
		point_constraints[i].weight = 0.0;
	    } else {
		// printf("[run_sfm] Constraining point %d\n", i);
		point_constraints[i].constrained = 1;
		point_constraints[i].weight = pt_constraint_weight;
		point_constraints[i].constraints[0] = Vx(pt_constraints[i]);
		point_constraints[i].constraints[1] = Vy(pt_constraints[i]);
		point_constraints[i].constraints[2] = Vz(pt_constraints[i]);
	    }
	}
    }

    /* Fill global param struct */
    global_params.num_cameras = num_cameras;
    global_params.num_points = num_pts;
    global_params.num_params_per_camera = cnp;

    global_params.est_focal_length = est_focal_length;
    global_params.const_focal_length = const_focal_length;
    global_params.estimate_distortion = undistort;
    global_params.explicit_camera_centers = explicit_camera_centers,
    
    global_params.global_params.f = 1.0;
    global_params.init_params = init_camera_params;

    global_last_ws = 
	safe_malloc(3 * num_cameras * sizeof(double), "global_last_ws");

    global_last_Rs = 
	safe_malloc(9 * num_cameras * sizeof(double), "global_last_ws");

    global_params.points = init_pts;

    for (i = 0; i < num_cameras; i++) {
	global_last_ws[3 * i + 0] = 0.0;
	global_last_ws[3 * i + 1] = 0.0;
	global_last_ws[3 * i + 2] = 0.0;

	memcpy(global_last_Rs + 9 * i, 
	       init_camera_params[i].R, 9 * sizeof(double));
    }

    /* Run sparse bundle adjustment */
#define MAX_ITERS 150 // 256
#define VERBOSITY 3

#ifdef SBA_V121
    if (fix_points == 0) {
        if (optimize_for_fisheye == 0) {
            sba_motstr_levmar(num_pts, num_cameras, ncons, 
                              vmask, params, cnp, 3, projections, NULL, 2, 
                              //remove NULL in prev line for sba v1.2.1
                              sfm_project_point3, NULL, 
                              (void *) (&global_params),
                              MAX_ITERS, VERBOSITY, opts, info,
                              use_constraints, constraints,
                              use_point_constraints,
                              point_constraints, Vout, Sout, Uout, Wout);
        } else {
            sba_motstr_levmar(num_pts, num_cameras, ncons, 
                              vmask, params, cnp, 3, projections, NULL, 2,
                              sfm_project_point2_fisheye, NULL, 
                              (void *) (&global_params),
                              MAX_ITERS, VERBOSITY, opts, info,
                              use_constraints, constraints,
                              use_point_constraints,
                              point_constraints, Vout, Sout, Uout, Wout); 
        }
    } else {
        if (optimize_for_fisheye == 0) {
            sba_mot_levmar(num_pts, num_cameras, ncons, 
                           vmask, params, cnp, projections, NULL, 2,
                           sfm_project_point3_mot, NULL, 
                           (void *) (&global_params),
                           MAX_ITERS, VERBOSITY, opts, info,
                           use_constraints, constraints);
        } else {
            sba_mot_levmar(num_pts, num_cameras, ncons, 
                           vmask, params, cnp, projections, NULL, 2,
                           sfm_project_point2_fisheye_mot, NULL, 
                           (void *) (&global_params),
                           MAX_ITERS, VERBOSITY, opts, info,
                           use_constraints, constraints);
        }
    }
#else
    if (fix_points == 0) {
	sba_motstr_levmar(num_pts, num_cameras, ncons, 
			  vmask, params, cnp, 3, projections, 2,
			  sfm_project_point2, NULL, (void *) (&global_params),
			  MAX_ITERS, VERBOSITY, opts, info, 
			  use_constraints, constraints, 
                          Vout, Sout, Uout, Wout);
    } else {
	sba_mot_levmar(num_pts, num_cameras, ncons, 
		       vmask, params, cnp, projections, 2,
		       sfm_mot_project_point, NULL, (void *) (&global_params),
		       MAX_ITERS, VERBOSITY, opts, info);
    }
#endif
    
    printf("[run_sfm] Number of iterations: %d\n", (int) info[5]);
    printf("info[6] = %0.3f\n", info[6]);
    
    /* Copy out the params */
    for (j = 0; j < num_cameras; j++) {
	double *dt = params + cnp * j + 0;
	double *w = params + cnp * j + 3;
	double Rnew[9];
        int c;

	/* Translation */
	init_camera_params[j].t[0] = dt[0];
	init_camera_params[j].t[1] = dt[1];
	init_camera_params[j].t[2] = dt[2];
	// init_camera_params[j].t[0] += dt[0];
	// init_camera_params[j].t[1] += dt[1];
	// init_camera_params[j].t[2] += dt[2];

	/* Rotation */
	rot_update(init_camera_params[j].R, w, Rnew);
	memcpy(init_camera_params[j].R, Rnew, 9 * sizeof(double));

	/* Focal length */
	if (est_focal_length) {
            c = 7;
#ifndef TEST_FOCAL
	    init_camera_params[j].f = params[cnp * j + 6];
#else
	    init_camera_params[j].f = 
                params[cnp * j + 6] / init_camera_params[j].f_scale;
#endif
        } else {
            c = 6;
        }
        
        if (undistort) {
#ifndef TEST_FOCAL
            init_camera_params[j].k[0] = params[cnp * j + c];
            init_camera_params[j].k[1] = params[cnp * j + c+1];
#else
            double scale = init_camera_params[j].k_scale;
            init_camera_params[j].k[0] = params[cnp * j + c] / scale;
            init_camera_params[j].k[1] = params[cnp * j + c+1] / scale;
#endif
        }

#ifdef TEST_FOCAL
        init_camera_params[j].f_scale = 1.0;
        init_camera_params[j].k_scale = 1.0;
#endif
    }

    base = num_camera_params;
    for (i = 0; i < num_pts; i++) {
	Vx(init_pts[i]) = params[base + 3 * i + 0];
	Vy(init_pts[i]) = params[base + 3 * i + 1];
	Vz(init_pts[i]) = params[base + 3 * i + 2];
    }

    // #define DEBUG_SFM
#ifdef DEBUG_SFM
    for (i = 0; i < num_cameras; i++) {
	int num_projs = 0;
	double error = 0.0;
	double error_max = 0.0;
	int idx_max = 0;
	double px_max = 0.0, py_max = 0.0;

	double K[9] = { init_camera_params[i].f, 0.0, 0.0, 
			0.0, init_camera_params[i].f, 0.0,
			0.0, 0.0, 1.0 };
	double w[3] = { 0.0, 0.0, 0.0 };
	double dt[3] = { init_camera_params[i].t[0],
			 init_camera_params[i].t[1],
			 init_camera_params[i].t[2] };

	// double dt[3] = { 0.0, 0.0, 0.0 };

	for (j = 0; j < num_pts; j++) {
	    double b[3], pr[2];
	    double dx, dy, dist;

	    if (!vmask[j * num_cameras + i])
		continue;

	    b[0] = Vx(init_pts[j]);
	    b[1] = Vy(init_pts[j]);
	    b[2] = Vz(init_pts[j]);

	    sfm_project(&(init_camera_params[i]), K, w, dt, b, pr,
			global_params.explicit_camera_centers);

	    dx = pr[0] - Vx(projections[j * num_cameras + i]);
	    dy = pr[1] - Vy(projections[j * num_cameras + i]);

	    dist = dx * dx + dy * dy;
	    error += dist;
	    
	    if (dist > error_max) {
		idx_max = j;
		error_max = dist;
		px_max = Vx(projections[j * num_cameras + i]);
		py_max = Vy(projections[j * num_cameras + i]);
	    }
	    
	    num_projs++;
	}

	printf("Camera %d:  error = %0.3f (%0.3f)\n", i,
	       error, sqrt(error / num_projs));
	printf("           error_max = %0.3f (%d)\n", sqrt(error_max), idx_max);
	printf("           proj = %0.3f, %0.3f\n", px_max, py_max);
    }
#endif /* DEBUG_SFM */

    free(params);

    // #ifndef SBA_V121
    if (use_constraints) {
	for (i = 0; i < num_cameras; i++) {
	    free(constraints[i].constraints);
	    free(constraints[i].constrained);
	    free(constraints[i].weights);
	}
	free(constraints);
    }

    free(global_last_ws);
    free(global_last_Rs);

    // #endif
}


static int global_num_points = 0;
static sfm_global_t *global_params = NULL;
static v3_t *global_points = NULL;
static v2_t *global_projections = NULL;
static int global_constrain_focal = 0;
static double global_init_focal = 0.0;
static double global_constrain_focal_weight = 0.0;
static double global_constrain_rd_weight = 0.0;
static int global_round = 0;

void camera_refine_residual(const int *m, const int *n, 
			    double *x, double *fvec, int *iflag) 
{
    int i;
    double error = 0.0, error2 = 0.0;

    for (i = 0; i < global_num_points; i++) {
	double pt[3] = { Vx(global_points[i]), 
			 Vy(global_points[i]),
			 Vz(global_points[i]) };

	double proj[2], dx, dy;
        
	sfm_project_point(0, i, x, pt, proj, (void *) global_params);

	dx = Vx(global_projections[i]) - proj[0];
	dy = Vy(global_projections[i]) - proj[1];

	fvec[2 * i + 0] = dx;
	fvec[2 * i + 1] = dy;

	if (*iflag == 0) {
	    error += dx * dx + dy * dy;
	    error2 += sqrt(dx * dx + dy * dy);
	}
    }

    if (global_constrain_focal == 1) {
	double focal_diff = global_init_focal - x[6];
	fvec[2 * global_num_points] = 
	    global_constrain_focal_weight * focal_diff;

        if (global_params->estimate_distortion) {
            fvec[2 * global_num_points + 1] = 
                -global_constrain_rd_weight * x[7];
            fvec[2 * global_num_points + 2] = 
                -global_constrain_rd_weight * x[8];
        }
    } else if (global_params->estimate_distortion) {
        fvec[2 * global_num_points + 0] = 
            -global_constrain_rd_weight * x[7];
        fvec[2 * global_num_points + 1] = 
            -global_constrain_rd_weight * x[8];
    }

    if (*iflag == 0) {
        if (global_params->estimate_distortion) {
            printf("  Round[%d]: RMS error = %0.8f [%0.8f], "
                   "f = %0.3f; %0.3e %0.3e\n", 
                   global_round, sqrt(error / global_num_points), 
                   error2 / global_num_points, x[6], x[7], x[8]);
        } else {
            if (global_params->est_focal_length) {
                printf("  Round[%d]: RMS error = %0.8f [%0.8f], f = %0.3f\n", 
                       global_round, sqrt(error / global_num_points),
                       error2 / global_num_points, x[6]);
            } else {
                printf("  Round[%d]: RMS error = %0.8f [%0.8f]\n", 
                       global_round, sqrt(error / global_num_points),
                       error2 / global_num_points);        
            }
        }
        
	if (global_constrain_focal == 1) {
	    printf("  Round[%d]: df = %0.3f\n", 
		   global_round, global_init_focal - x[6]);
	}
	
	global_round++;
    }
}

/* Refine the position of a single camera */
void camera_refine(int num_points, v3_t *points, v2_t *projs, 
		   camera_params_t *params, int adjust_focal, 
                   int estimate_distortion)
{
    if (adjust_focal) {
        int num_camera_params = 7;
	double x[9] = { params->t[0], params->t[1], params->t[2], 
			0.0, 0.0, 0.0, params->f, params->k[0], params->k[1] };
	double Rnew[9];

	sfm_global_t globs;
	int focal_constraint = 0;
    
        if (estimate_distortion)
            num_camera_params += 2;

	globs.num_cameras = 1;
	globs.num_points = num_points;
	globs.est_focal_length = 1;
	globs.const_focal_length = 0;
	globs.explicit_camera_centers = 1;
	globs.global_params.f = params->f;
	globs.init_params = params;
        globs.estimate_distortion = estimate_distortion;

	global_num_points = num_points;
	global_params = &globs;
	global_points = points;
	global_projections = projs;
	global_round = 0;

	if (params->constrained[6]) {
	    printf("[camera_refine] Constraining focal length to %0.3f "
		   "(weight: %0.3f)\n",
		   params->constraints[6], 
		   num_points * params->weights[6]);
	    focal_constraint = 1;
	    global_init_focal = params->constraints[6];
	    global_constrain_focal = 1;
	    global_constrain_focal_weight = 
		1.0e0 /*1.0e1*/ * num_points * params->weights[6];
	} else {
	    focal_constraint = 0;
	    global_init_focal = 0.0;
	    global_constrain_focal = 0;
	    global_constrain_focal_weight = 0.0;	    
	}

        if (estimate_distortion) {       
            global_constrain_rd_weight = 0.05 * num_points; 
                // 1.0e-1 * num_points;
        }
    
	lmdif_driver2(camera_refine_residual, 
                      2 * num_points + focal_constraint + 
                      2 * estimate_distortion, 
                      num_camera_params, x, 1.0e-12);

	/* Copy out the parameters */
	memcpy(params->t, x + 0, 3 * sizeof(double));
	rot_update(params->R, x + 3, Rnew);
	memcpy(params->R, Rnew, 9 * sizeof(double));    
	params->f = x[6];

        if (estimate_distortion) {
            params->k[0] = x[7];
            params->k[1] = x[8];
        }
    } else {
	double x[6] = { params->t[0], params->t[1], params->t[2], 
			0.0, 0.0, 0.0 };
	double Rnew[9];

	sfm_global_t globs;
    
	globs.num_cameras = 1;
	globs.num_points = num_points;
	globs.est_focal_length = 0;
	globs.const_focal_length = 1;
	globs.explicit_camera_centers = 1;
	globs.global_params.f = params->f;
	globs.init_params = params;
        globs.estimate_distortion = estimate_distortion;

	global_num_points = num_points;
	global_params = &globs;
	global_points = points;
	global_projections = projs;
	global_round = 0;

	global_init_focal = 0.0;
	global_constrain_focal = 0;
	global_constrain_focal_weight = 0.0;	    
    
	lmdif_driver2(camera_refine_residual, 2 * num_points, 6, x, 1.0e-12);

	/* Copy out the parameters */
	memcpy(params->t, x + 0, 3 * sizeof(double));
	rot_update(params->R, x + 3, Rnew);
	memcpy(params->R, Rnew, 9 * sizeof(double));    
    }
}


