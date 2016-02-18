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

/* triangulate.c */
/* Triangulate two image points */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "triangulate.h"
#include "vector.h"

static v2_t global_p, global_q;
static double *global_R0, *global_t0, *global_R1, *global_t1;

void quick_svd(double *E, double *U, double *S, double *VT) {
    double e1[3] = { E[0], E[3], E[6] };
    double e2[3] = { E[1], E[4], E[7] };
    double e3[3] = { E[2], E[5], E[8] };

    double e1e2[3], e1e3[3], e2e3[3];
    double e1e2_mag, e1e3_mag, e2e3_mag;

    double u_a[3], u_b[3], u_c[3];
    double v_a[3], v_b[3], v_c[3];
    double tmp[3], mag;
    
    double e1_mag = matrix_norm(3, 1, e1);
    double e2_mag = matrix_norm(3, 1, e2);

    matrix_cross(e1, e2, e1e2);
    matrix_cross(e1, e3, e1e3);
    matrix_cross(e2, e3, e2e3);

    e1e2_mag = matrix_norm(3, 1, e1e2);
    e1e3_mag = matrix_norm(3, 1, e1e3);
    e2e3_mag = matrix_norm(3, 1, e2e3);

    if (e1e2_mag >= e1e3_mag && e1e2_mag >= e1e3_mag) {
	matrix_scale(3, 1, e1, 1.0 / e1_mag, v_a);
	matrix_scale(3, 1, e1e2, 1.0 / e1e2_mag, v_c);
	matrix_cross(v_c, v_a, v_b);
    } else if (e1e3_mag > e1e2_mag && e1e3_mag > e2e3_mag) {
	matrix_scale(3, 1, e1, 1.0 / e1_mag, v_a);
	matrix_scale(3, 1, e1e3, 1.0 / e1e3_mag, v_b);
	matrix_cross(v_b, v_a, v_c);
    } else if (e2e3_mag > e1e2_mag && e2e3_mag > e1e3_mag) {
	matrix_scale(3, 1, e2, 1.0 / e2_mag, v_a);
	matrix_scale(3, 1, e2e3, 1.0 / e2e3_mag, v_a);
	matrix_cross(v_b, v_c, v_a);
    }

    matrix_product331(E, v_a, tmp);
    mag = matrix_norm(3, 1, tmp);
    matrix_scale(3, 1, tmp, 1.0 / mag, u_a);

    matrix_product331(E, v_b, tmp);
    mag = matrix_norm(3, 1, tmp);
    matrix_scale(3, 1, tmp, 1.0 / mag, u_b);

    matrix_cross(u_a, u_b, u_c);
}

/* Project a point onto an image */
v2_t project(double *R, double *t0, double *P) {
    double tmp[3], tmp2[3];
    v2_t result;

    /* Rigid transform */
    matrix_product331(R, P, tmp);
    matrix_sum(3, 1, 3, 1, tmp, t0, tmp2);
    
    /* Perspective division */
    Vx(result) = tmp2[0] / tmp2[2];
    Vy(result) = tmp2[1] / tmp2[2];

    return result;
}

void triangulation_residual(const int *m, const int *n, double *x, 
			    double *fvec, double *iflag) 
{
    /* Project the point into the two views */
    v2_t p = project(global_R0, global_t0, x);
    v2_t q = project(global_R1, global_t1, x);
    
    fvec[0] = Vx(global_p) - Vx(p);
    fvec[1] = Vy(global_p) - Vy(p);
    fvec[2] = Vx(global_q) - Vx(q);
    fvec[3] = Vy(global_q) - Vy(q);
}

static int global_num_points;
static double *global_Rs = NULL;
static double *global_ts = NULL;
static v2_t *global_ps;

void triangulate_n_residual(const int *m, const int *n, 
			    double *x, double *fvec, double *iflag) 
{
    int i;

    for (i = 0; i < global_num_points; i++) {
	int Roff = 9 * i;
	int toff = 3 * i;

	/* Project the point into the view */
	v2_t p = project(global_Rs + Roff, global_ts + toff, x);
    
	fvec[2 * i + 0] = Vx(global_ps[i]) - Vx(p);
	fvec[2 * i + 1] = Vy(global_ps[i]) - Vy(p);
    }
}


/* Find the point with the smallest squared projection error */
v3_t triangulate_n_refine(v3_t pt, int num_points, 
			  v2_t *p, double *R, double *t, double *error_out) 
{
    int num_eqs = 2 * num_points;
    int num_vars = 3;

    double x[3] = { Vx(pt), Vy(pt), Vz(pt) };
    double error;

    int i;

    /* Run a non-linear optimization to polish the result */
    global_num_points = num_points;
    global_ps = p;
    global_Rs = R;  global_ts = t;
    lmdif_driver(triangulate_n_residual, num_eqs, num_vars, x, 1.0e-5);

    error = 0.0;
    for (i = 0; i < num_points; i++) {
	double dx, dy;
	int Roff = 9 * i;
	int toff = 3 * i;
	double pp[3];

	/* Compute projection error */
	matrix_product331(R + Roff, x, pp);
	pp[0] += t[toff + 0];
	pp[1] += t[toff + 1];
	pp[2] += t[toff + 2];

	dx = pp[0] / pp[2] - Vx(p[i]);
	dy = pp[1] / pp[2] - Vy(p[i]);
	error += dx * dx + dy * dy;
    }

    error = sqrt(error / num_points);

    // printf("[triangulate_n] Error [after polishing]: %0.3e\n", error);

    if (error_out != NULL) {
	*error_out = error;
    }

    return v3_new(x[0], x[1], x[2]);
}


/* Find the point with the smallest squared projection error */
v3_t triangulate_n(int num_points, 
		   v2_t *p, double *R, double *t, double *error_out) 
{
    int num_eqs = 2 * num_points;
    int num_vars = 3;

    double *A = (double *) malloc(sizeof(double) * num_eqs * num_vars);
    double *b = (double *) malloc(sizeof(double) * num_eqs);
    double *x = (double *) malloc(sizeof(double) * num_vars);

    int i;
    double error;

    v3_t r;

    for (i = 0; i < num_points; i++) {
	int Roff = 9 * i;
	int row = 6 * i;
	int brow = 2 * i;
	int toff = 3 * i;

	A[row + 0] = R[Roff + 0] - Vx(p[i]) * R[Roff + 6];  
	A[row + 1] = R[Roff + 1] - Vx(p[i]) * R[Roff + 7];  
	A[row + 2] = R[Roff + 2] - Vx(p[i]) * R[Roff + 8];

	A[row + 3] = R[Roff + 3] - Vy(p[i]) * R[Roff + 6];  
	A[row + 4] = R[Roff + 4] - Vy(p[i]) * R[Roff + 7];  
	A[row + 5] = R[Roff + 5] - Vy(p[i]) * R[Roff + 8];

	b[brow + 0] = t[toff + 2] * Vx(p[i]) - t[toff + 0];
	b[brow + 1] = t[toff + 2] * Vy(p[i]) - t[toff + 1];
    }
    
    /* Find the least squares result */
    dgelsy_driver(A, b, x, num_eqs, num_vars, 1);

    error = 0.0;
    for (i = 0; i < num_points; i++) {
	double dx, dy;
	int Roff = 9 * i;
	int toff = 3 * i;
	double pp[3];

	/* Compute projection error */
	matrix_product331(R + Roff, x, pp);
	pp[0] += t[toff + 0];
	pp[1] += t[toff + 1];
	pp[2] += t[toff + 2];

	dx = pp[0] / pp[2] - Vx(p[i]);
	dy = pp[1] / pp[2] - Vy(p[i]);
	error += dx * dx + dy * dy;
    }

    error = sqrt(error / num_points);

    // printf("[triangulate_n] Error [before polishing]: %0.3e\n", error);

    /* Run a non-linear optimization to refine the result */
    global_num_points = num_points;
    global_ps = p;
    global_Rs = R;  global_ts = t;
    lmdif_driver(triangulate_n_residual, num_eqs, num_vars, x, 1.0e-5);

    error = 0.0;
    for (i = 0; i < num_points; i++) {
	double dx, dy;
	int Roff = 9 * i;
	int toff = 3 * i;
	double pp[3];

	/* Compute projection error */
	matrix_product331(R + Roff, x, pp);
	pp[0] += t[toff + 0];
	pp[1] += t[toff + 1];
	pp[2] += t[toff + 2];

	dx = pp[0] / pp[2] - Vx(p[i]);
	dy = pp[1] / pp[2] - Vy(p[i]);
	error += dx * dx + dy * dy;
    }

    error = sqrt(error / num_points);

    // printf("[triangulate_n] Error [after polishing]: %0.3e\n", error);

    if (error_out != NULL) {
	*error_out = error;
    }

    r = v3_new(x[0], x[1], x[2]);

    free(A);
    free(b);
    free(x);

    return r;
}


/* Find the point with the smallest squared projection error */
v3_t triangulate(v2_t p, v2_t q, 
		 double *R0, double *t0, 
		 double *R1, double *t1, double *error) 
{
    double A[12];
    double b[4];
    double x[3];

	double dx1, dx2, dy1, dy2;

    A[0] = R0[0] - Vx(p) * R0[6];  
    A[1] = R0[1] - Vx(p) * R0[7];  
    A[2] = R0[2] - Vx(p) * R0[8];
    
    A[3] = R0[3] - Vy(p) * R0[6];  
    A[4] = R0[4] - Vy(p) * R0[7];  
    A[5] = R0[5] - Vy(p) * R0[8];

    A[6] = R1[0] - Vx(q) * R1[6];  
    A[7] = R1[1] - Vx(q) * R1[7];  
    A[8] = R1[2] - Vx(q) * R1[8];

    A[9] = R1[3] - Vy(q) * R1[6];  
    A[10] = R1[4] - Vy(q) * R1[7];  
    A[11] = R1[5] - Vy(q) * R1[8];

    b[0] = t0[2] * Vx(p) - t0[0];
    b[1] = t0[2] * Vy(p) - t0[1];
    b[2] = t1[2] * Vx(q) - t1[0];
    b[3] = t1[2] * Vy(q) - t1[1];

    /* Find the least squares result */
    dgelsy_driver(A, b, x, 4, 3, 1);

    /* Run a non-linear optimization to refine the result */
    global_p = p;
    global_q = q;
    global_R0 = R0;  global_t0 = t0;
    global_R1 = R1;  global_t1 = t1;
    lmdif_driver(triangulation_residual, 4, 3, x, 1.0e-10);

    if (error != NULL) {
	double pp[3], qp[3];

	/* Compute projection error */
	matrix_product331(R0, x, pp);
	pp[0] += t0[0];
	pp[1] += t0[1];
	pp[2] += t0[2];
		
	matrix_product331(R1, x, qp);
	qp[0] += t1[0];
	qp[1] += t1[1];
	qp[2] += t1[2];	

	dx1 = pp[0] / pp[2] - Vx(p);
	dy1 = pp[1] / pp[2] - Vy(p);

	dx2 = qp[0] / qp[2] - Vx(q);
	dy2 = qp[1] / qp[2] - Vy(q);

	*error = dx1 * dx1 + dy1 * dy1 + dx2 * dx2 + dy2 * dy2;
    }

    return v3_new(x[0], x[1], x[2]);
}

/* Given an F matrix, two calibration matrices, and a point correspondence, find R and t */
void find_extrinsics(double *F, double *K1, double *K2, 
                     v2_t p1, v2_t p2, double *R, double *t) {
    double E[9];
    double K1_inv[9], K2_inv[9], tmp[9];
    
    /* Find the essential matrix */
    matrix_invert(3, K1, K1_inv);
    matrix_invert(3, K2, K2_inv);
    
    matrix_product33(F, K1, tmp);
    matrix_transpose_product(3, 3, 3, 3, K2, tmp, E);

    find_extrinsics_essential(E, p1, p2, R, t);
}

/* Given an E matrix and a point correspondence, find R and t */
int find_extrinsics_essential(double *E, v2_t p1, v2_t p2, 
                              double *R, double *t)
{
    double tmp[9], tmp2[3], Qv[3];
    double U[9], S[3], VT[9];
    double tu[3], Ra[9], Rb[9];

    double D[9] = 
	{  0.0, 1.0, 0.0,
	  -1.0, 0.0, 0.0,
	   0.0, 0.0, 1.0 };

    double DT[9] = 
	{  0.0, -1.0, 0.0,
	   1.0,  0.0, 0.0,
	   0.0,  0.0, 1.0 };

    double I[9] = 
	{  1.0, 0.0, 0.0, 
	   0.0, 1.0, 0.0,
	   0.0, 0.0, 1.0 };
    
    double t0[3] = { 0.0, 0.0, 0.0 };

    v3_t Q, PQ;
    double c1, c2;
    double error;

    /* Now find the SVD of E */
    dgesvd_driver(3, 3, E, U, S, VT);

#if 0
    printf("S[0] = %0.3e\n", S[0]);
    printf("S[1] = %0.3e\n", S[1]);
    printf("S[2] = %0.3e\n", S[2]);
#endif

    /* Now find R and t */
    tu[0] = U[2];  tu[1] = U[5];  tu[2] = U[8];

    matrix_product33(U, D, tmp);
    matrix_product33(tmp, VT, Ra);
    matrix_product33(U, DT, tmp);

    matrix_product33(tmp, VT, Rb);

    if (matrix_determinant3(Ra) < 0.0) {
	// printf("flipping...\n");
	matrix_scale(3, 3, Ra, -1.0, Ra);
    }

    if (matrix_determinant3(Rb) < 0.0) {
	// printf("flopping...\n");
	matrix_scale(3, 3, Rb, -1.0, Rb);
    }

    /* Figure out which configuration is correct using the supplied
     * point */

    Q = triangulate(p1, p2, I, t0, Ra, tu, &error);
    Qv[0] = Vx(Q), Qv[1] = Vy(Q), Qv[2] = Vz(Q);
    matrix_product331(Ra, Qv, tmp);
    matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
    PQ = v3_new(tmp2[0], tmp2[1], tmp2[2]);

    c1 = Vz(Q);
    c2 = Vz(PQ);

    if (c1 < 0 && c2 < 0) {
	memcpy(R, Ra, 9 * sizeof(double));
	t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
    } else if (c1 > 0 && c2 > 0) {
	memcpy(R, Ra, 9 * sizeof(double));
	t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
        Q = triangulate(p1, p2, I, t0, R, t, &error);
    } else {
        /* Triangulate again */
        Q = triangulate(p1, p2, I, t0, Rb, tu, &error);
        Qv[0] = Vx(Q), Qv[1] = Vy(Q), Qv[2] = Vz(Q);
        matrix_product331(Rb, Qv, tmp);
        matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
        PQ = v3_new(tmp2[0], tmp2[1], tmp2[2]);

        c1 = Vz(Q);
        c2 = Vz(PQ);

        if (c1 < 0 && c2 < 0) {
            memcpy(R, Rb, 9 * sizeof(double));
            t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
        } else if (c1 > 0 && c2 > 0) {
            memcpy(R, Rb, 9 * sizeof(double));
            t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
            Q = triangulate(p1, p2, I, t0, R, t, &error);
        } else {
            printf("[find_extrinsics] Error: no case found!\n");
            return 0;
        }
    }

    matrix_product331(R, Q.p, tmp2);
    matrix_sum(3, 1, 3, 1, tmp2, t, tmp2);

    printf("[find_extrinsics] error: %0.3e\n", error);
    // printf("  %0.3f => %0.3f\n", Vx(p2), tmp2[0] / tmp2[2]);
    // printf("  %0.3f => %0.3f\n", Vy(p2), tmp2[1] / tmp2[2]);
    // printf("  %0.3f\n", tmp2[2]);

    return 1;
}

/* Given an E matrix and a point correspondence, find R and t */
int find_extrinsics_essential_multipt(double *E, int n,
                                      v2_t *p1, v2_t *p2, 
                                      double *R, double *t)
{
    double tmp[9], tmp2[3], Qv[3];
    double U[9], S[3], VT[9];
    double tu[3], Ra[9], Rb[9];

    double D[9] = 
	{  0.0, 1.0, 0.0,
	  -1.0, 0.0, 0.0,
	   0.0, 0.0, 1.0 };

    double DT[9] = 
	{  0.0, -1.0, 0.0,
	   1.0,  0.0, 0.0,
	   0.0,  0.0, 1.0 };

    double I[9] = 
	{  1.0, 0.0, 0.0, 
	   0.0, 1.0, 0.0,
	   0.0, 0.0, 1.0 };
    
    double t0[3] = { 0.0, 0.0, 0.0 };

    v3_t Q, PQ;
    double c1, c2;
    double error;

    int i;
    int c1_pos = 0, c1_neg = 0;
    int c2_pos = 0, c2_neg = 0;

    /* Now find the SVD of E */
    dgesvd_driver(3, 3, E, U, S, VT);

#if 0
    printf("S[0] = %0.3e\n", S[0]);
    printf("S[1] = %0.3e\n", S[1]);
    printf("S[2] = %0.3e\n", S[2]);
#endif

    /* Now find R and t */
    tu[0] = U[2];  tu[1] = U[5];  tu[2] = U[8];

    matrix_product33(U, D, tmp);
    matrix_product33(tmp, VT, Ra);
    matrix_product33(U, DT, tmp);

    matrix_product33(tmp, VT, Rb);

    if (matrix_determinant3(Ra) < 0.0) {
	// printf("flipping...\n");
	matrix_scale(3, 3, Ra, -1.0, Ra);
    }

    if (matrix_determinant3(Rb) < 0.0) {
	// printf("flopping...\n");
	matrix_scale(3, 3, Rb, -1.0, Rb);
    }

    /* Figure out which configuration is correct using the supplied
     * points */

    for (i = 0; i < n; i++) {
        Q = triangulate(p1[i], p2[i], I, t0, Ra, tu, &error);
        Qv[0] = Vx(Q), Qv[1] = Vy(Q), Qv[2] = Vz(Q);
        matrix_product331(Ra, Qv, tmp);
        matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
        PQ = v3_new(tmp2[0], tmp2[1], tmp2[2]);

        c1 = Vz(Q);
        c2 = Vz(PQ);

        if (c1 > 0)
            c1_pos++;
        else
            c1_neg++;
        
        if (c2 > 0)
            c2_pos++;
        else
            c2_neg++;
    }

    if (c1_pos < c1_neg && c2_pos < c2_neg) {
	memcpy(R, Ra, 9 * sizeof(double));
	t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
    } else if (c1_pos > c1_neg && c2_pos > c2_neg) {
	memcpy(R, Ra, 9 * sizeof(double));
	t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
    } else {
        /* Triangulate again */
        c1_pos = c1_neg = c2_pos = c2_neg = 0;

        for (i = 0; i < n; i++) {
            Q = triangulate(p1[i], p2[i], I, t0, Rb, tu, &error);
            Qv[0] = Vx(Q), Qv[1] = Vy(Q), Qv[2] = Vz(Q);
            matrix_product331(Rb, Qv, tmp);
            matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
            PQ = v3_new(tmp2[0], tmp2[1], tmp2[2]);

            c1 = Vz(Q);
            c2 = Vz(PQ);

            if (c1 > 0)
                c1_pos++;
            else
                c1_neg++;
        
            if (c2 > 0)
                c2_pos++;
            else
                c2_neg++;
        }

        if (c1_pos < c1_neg && c2_pos < c2_neg) {
            memcpy(R, Rb, 9 * sizeof(double));
            t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
        } else if (c1_pos > c1_neg && c2_pos > c2_neg) {
            memcpy(R, Rb, 9 * sizeof(double));
            t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
        } else {
            fprintf(stderr, "[find_extrinsics] Error: no case found!\n");
            return 0;
        }
    }

    return 1;
}


static v3_t *condition_points_3D(int num_points, v3_t *pts, double *T) {
    v3_t *pts_new = (v3_t *) malloc(sizeof(v3_t) * num_points);
    v3_t *pts_zero_mean = (v3_t *) malloc(sizeof(v3_t) * num_points);
    double U[9], S[3], VT[9], U16[9], Sinv[16], Ttrans[16], action[16];

    v3_t mean = v3_mean(num_points, pts);

    double total_dist = 0.0;
    double avg_dist;
    double factor;
    int i;

    for (i = 0; i < num_points; i++) {
	double dx = Vx(pts[i]) - Vx(mean);
	double dy = Vy(pts[i]) - Vy(mean);
	double dz = Vz(pts[i]) - Vz(mean);
	
	pts_zero_mean[i] = v3_new(Vx(pts[i]) - Vx(mean),
				  Vy(pts[i]) - Vy(mean),
				  Vz(pts[i]) - Vz(mean));

	total_dist += sqrt(dx * dx + dy * dy + dz * dz);
    }

    v3_svd(num_points, pts_zero_mean, U, S, VT);

    Sinv[0] = sqrt(3.0) / sqrt(S[0]);  Sinv[1] = 0.0;  Sinv[2] = 0.0;  Sinv[3] = 0.0;
    Sinv[4] = 0.0;  Sinv[5] = sqrt(3.0) / sqrt(S[1]);  Sinv[6] = 0.0;  Sinv[7] = 0.0;
    Sinv[8] = 0.0;  Sinv[9] = 0.0;  Sinv[10] = sqrt(3.0) / sqrt(S[2]); Sinv[11] = 0.0;
    Sinv[12] = 0.0; Sinv[13] = 0.0; Sinv[14] = 0.0; Sinv[15] = 1.0;

    matrix_ident(4, U16);
    memcpy(U16 + 0, U + 0, 3 * sizeof(double));
    memcpy(U16 + 4, U + 3, 3 * sizeof(double));
    memcpy(U16 + 8, U + 6, 3 * sizeof(double));

    matrix_ident(4, Ttrans);
    Ttrans[3]  = -Vx(mean);
    Ttrans[7]  = -Vy(mean);
    Ttrans[11] = -Vz(mean);

    avg_dist = total_dist / num_points;
    factor = sqrt(3.0) / avg_dist;

    matrix_transpose_product2(4, 4, 4, 4, Sinv, U16, action);
    matrix_product44(action, Ttrans, T);

    for (i = 0; i < num_points; i++) {
	double x = factor * (Vx(pts[i]) - Vx(mean));
	double y = factor * (Vy(pts[i]) - Vy(mean));
	double z = factor * (Vz(pts[i]) - Vz(mean));

	double pt[4] = { Vx(pts[i]), Vy(pts[i]), Vz(pts[i]), 1.0 };
	double Tpt[4];
	
	matrix_product441(T, pt, Tpt);
	
	pts_new[i] = v3_new(Tpt[0], Tpt[1], Tpt[2]);
    }

    free(pts_zero_mean);

    return pts_new;
}

static v2_t *condition_points_2D(int num_points, v2_t *pts, double *T) {
    v2_t *pts_new = (v2_t *) malloc(sizeof(v2_t) * num_points);

    v2_t mean = v2_mean(num_points, pts);
    double total_dist = 0.0;
    double avg_dist;
    double factor;
    int i;

    for (i = 0; i < num_points; i++) {
	double dx = Vx(pts[i]) - Vx(mean);
	double dy = Vy(pts[i]) - Vy(mean);
	total_dist += sqrt(dx * dx + dy * dy);
    }

    avg_dist = total_dist / num_points;
    factor = sqrt(2.0) / avg_dist;

    for (i = 0; i < num_points; i++) {
	double x = factor * (Vx(pts[i]) - Vx(mean));
	double y = factor * (Vy(pts[i]) - Vy(mean));
	pts_new[i] = v2_new(x, y);
    }

    T[0] = factor;  T[1] = 0.0;     T[2] = -factor * Vx(mean);
    T[3] = 0.0;     T[4] = factor;  T[5] = -factor * Vy(mean);
    T[6] = 0.0;      T[7] = 0.0;    T[8] = 1.0;

    return pts_new;
}

/* Solve for a 3x4 projection matrix, given a set of 3D points and 2D
 * projections */
int find_projection_3x4(int num_pts, v3_t *points, v2_t *projs, double *P) {
    if (num_pts < 6) {
	printf("[find_projection_3x4] Need at least 6 points!\n");
	return -1;
    } else {

	// #define _CONDITION_
#ifdef _CONDITION_
	double Tpoints[16];
	v3_t *points_new = condition_points_3D(num_pts, points, Tpoints);
	
	double Tprojs[9];
	v2_t *projs_new = condition_points_2D(num_pts, projs, Tprojs);

	double Tprojs_inv[9];
	double Ptmp[12];
#else
	v3_t *points_new = points;
	v2_t *projs_new = projs;
#endif

	int num_eqns = 2 * num_pts;
	int num_vars = 11;

	double *A = malloc(sizeof(double) * num_eqns * num_vars);
	double *b = malloc(sizeof(double) * num_eqns);
	double X[11];

	double error = 0.0;
    
	int i;

	for (i = 0; i < num_pts; i++) {
	    double *row1 = A + 2 * i * num_vars;
	    double *row2 = A + (2 * i + 1) * num_vars;
	    
	    row1[0]  = Vx(points_new[i]);
	    row1[1]  = Vy(points_new[i]);
	    row1[2]  = Vz(points_new[i]);
	    row1[3]  = 1.0;
	
	    row1[4]  = 0.0;
	    row1[5]  = 0.0;
	    row1[6]  = 0.0;
	    row1[7]  = 0.0;
	
	    row1[8]  = Vx(projs_new[i]) * Vx(points_new[i]);
	    row1[9]  = Vx(projs_new[i]) * Vy(points_new[i]);
	    row1[10] = Vx(projs_new[i]) * Vz(points_new[i]);
	
	    b[2 * i] = -Vx(projs_new[i]);


	    row2[0]  = 0.0;
	    row2[1]  = 0.0;
	    row2[2]  = 0.0;
	    row2[3]  = 0.0;

	    row2[4]  = Vx(points_new[i]);
	    row2[5]  = Vy(points_new[i]);
	    row2[6]  = Vz(points_new[i]);
	    row2[7]  = 1.0;	

	    row2[8]  = Vy(projs_new[i]) * Vx(points_new[i]);
	    row2[9]  = Vy(projs_new[i]) * Vy(points_new[i]);
	    row2[10] = Vy(projs_new[i]) * Vz(points_new[i]);
	
	    b[2 * i + 1] = -Vy(projs_new[i]);
	}

	dgelsy_driver(A, b, X, num_eqns, num_vars, 1);

	memcpy(P, X, sizeof(double) * 11);
	P[11] = 1.0;

#ifdef _CONDITION_
	matrix_invert(3, Tprojs, Tprojs_inv);
	matrix_product(3, 3, 3, 4, Tprojs_inv, P, Ptmp);
	matrix_product(3, 4, 4, 4, Ptmp, Tpoints, P);
	
	matrix_scale(3, 4, P, 1.0 / P[11], P);
#endif

	for (i = 0; i < num_pts; i++) {
	    double pt[4] = { Vx(points[i]), 
			     Vy(points[i]), 
			     Vz(points[i]), 1.0 };
	    double pr[3];
	    double dx, dy, dist;

	    matrix_product341(P, pt, pr);
	    pr[0] /= -pr[2];
	    pr[1] /= -pr[2];
	    
	    dx = pr[0] - Vx(projs[i]);
	    dy = pr[1] - Vy(projs[i]);

	    dist = dx * dx + dy * dy;

	    error += dist;
	}

	// printf("[find_projection_3x4] Average error is %0.3f\n", 
	//       error / num_pts);

	free(A);
	free(b);

#ifdef _CONDITION_
	free(points_new);
	free(projs_new);
#endif

	return 0;
    }
}

static int global_num_pts;
static v3_t *global_points;
static v2_t *global_projs;

static void projection_residual(const int *m, const int *n, double *x, 
				double *fvec, double *iflag) 
{
    int i;
    
    double P[12];
    memcpy(P, x, sizeof(double) * 11);
    P[11] = 1.0;

    for (i = 0; i < global_num_pts; i++) {
	double pt[4] = { Vx(global_points[i]), 
			 Vy(global_points[i]), 
			 Vz(global_points[i]), 1.0 };

	double pr[3];
	double dx, dy;
	
	matrix_product341(P, pt, pr);
	pr[0] /= -pr[2];
	pr[1] /= -pr[2];
	    
	dx = pr[0] - Vx(global_projs[i]);
	dy = pr[1] - Vy(global_projs[i]);

	fvec[2 * i + 0] = dx;
	fvec[2 * i + 1] = dy;
    }
}

/* Solve for a 3x4 projection matrix, given a set of 3D points and 2D
 * projections using non-linear optimization */
int find_projection_3x4_nonlinear(int num_pts, v3_t *points, v2_t *projs, 
				  double *Pin, double *Pout) 
{
    if (num_pts < 6) {
	printf("[find_projection_3x4_nonlinear] Need at least 6 points!\n");
	return -1;
    } else {
	int num_eqns = 2 * num_pts;
	int num_vars = 11;
	double x[11];

	global_num_pts = num_pts;
	global_points = points;
	global_projs = projs;

	memcpy(x, Pin, sizeof(double) * 11);
	lmdif_driver(projection_residual, num_eqns, num_vars, x, 1.0e-5);

	memcpy(Pout, x, sizeof(double) * 11);
	Pout[11] = 1.0;

	return 0;
    }
}


/* Solve for a 3x4 projection matrix using RANSAC, given a set of 3D
 * points and 2D projections */
int find_projection_3x4_ransac(int num_pts, v3_t *points, v2_t *projs, 
			       double *P, 
			       int ransac_rounds, double ransac_threshold) 
{
    if (num_pts < 6) {
	printf("[find_projection_3x4_ransac] Error: need at least 6 points!\n");
	return -1;
    } else {
#define MIN_PTS 6
	// const int min_pts = 6;
	int *inliers = (int *) malloc(sizeof(int) * num_pts);
	int indices[MIN_PTS];
	int round, i, j;
	int max_inliers = 0;
	double max_error = 0.0;
	double Pbest[12];
	int num_inliers = 0, num_inliers_new = 0;
	v3_t *pts_final = NULL;
	v2_t *projs_final = NULL;
	double Plinear[12];

	double Rinit[9];
	double triangular[9], orthogonal[9];
	int neg, sign;

	double thresh_sq = ransac_threshold * ransac_threshold;
	double error = 0.0;

	int num_inliers_polished = 0;

	for (round = 0; round < ransac_rounds; round++) {
	    v3_t pts_inner[MIN_PTS];
	    v2_t projs_inner[MIN_PTS];
	    double Ptmp[12];

	    num_inliers = 0;
	    for (i = 0; i < MIN_PTS; i++) {
		int redo = 0;
		int idx;
                int redo_count = 0;

		do {
                    if (redo_count > 10000) {
                        free(inliers);
                        return -1;
                    }

		    idx = rand() % num_pts;
		    
		    redo = 0;
		    for (j = 0; j < i; j++) {
			if (idx == indices[j]) {
			    redo = 1;
			    break;
			} else if (Vx(projs[idx]) == Vx(projs[indices[j]]) && 
                                   Vy(projs[idx]) == Vy(projs[indices[j]])) {
                            redo = 1;
                        }
		    }

                    redo_count++;
		} while(redo);
	    
		indices[i] = idx;
		pts_inner[i] = points[idx];
		projs_inner[i] = projs[idx];
	    }

	    /* Solve for the parameters */
	    find_projection_3x4(MIN_PTS, pts_inner, projs_inner, Ptmp);

#if 1
	    /* Fix the sign on the P matrix */
            memcpy(Rinit + 0, Ptmp + 0, 3 * sizeof(double));
            memcpy(Rinit + 3, Ptmp + 4, 3 * sizeof(double));
            memcpy(Rinit + 6, Ptmp + 8, 3 * sizeof(double));

            dgerqf_driver(3, 3, Rinit, triangular, orthogonal);	    

	    /* Check the parity along the diagonal */
	    neg = 
		(triangular[0] < 0.0) + 
		(triangular[4] < 0.0) + 
		(triangular[8] < 0.0);

	    if ((neg % 2) == 1) {
		sign = -1;
	    } else {
		sign = 1;
	    }
#endif
	    
	    /* Count the number of inliers */
	    error = 0.0;
	    for (i = 0; i < num_pts; i++) {
		double pt[4] = { Vx(points[i]), 
				 Vy(points[i]), 
				 Vz(points[i]), 1.0 };
		double pr[3];
		double dx, dy, dist;

		matrix_product341(Ptmp, pt, pr);

		/* Check cheirality */
		if (sign * pr[2] > 0.0) 
		    continue;

		pr[0] /= -pr[2];
		pr[1] /= -pr[2];
	    
		dx = pr[0] - Vx(projs[i]);
		dy = pr[1] - Vy(projs[i]);

		dist = dx * dx + dy * dy;

		if (dist < thresh_sq) {
		    inliers[num_inliers] = i;
		    num_inliers++;
		    error += dist;
		}
	    }
	    
	    if (num_inliers > max_inliers) {
		memcpy(Pbest, Ptmp, sizeof(double) * 12);
		max_error = error;
		max_inliers = num_inliers;
	    }
	}
	
	memcpy(P, Pbest, sizeof(double) * 12);

	printf("[find_projection_3x4_ransac] num_inliers = %d (out of %d)\n",
	       max_inliers, num_pts);
	printf("[find_projection_3x4_ransac] error = %0.3f\n", 
	       sqrt(max_error / max_inliers));

        if (max_inliers < 6) {
            printf("[find_projection_3x4_ransac] "
                   "Too few inliers to continue.\n");
            
            free(inliers);

            return -1;
        }
	
	/* Do the final least squares minimization */

#if 1
	/* Fix the sign on the P matrix */
	memcpy(Rinit + 0, Pbest + 0, 3 * sizeof(double));
	memcpy(Rinit + 3, Pbest + 4, 3 * sizeof(double));
	memcpy(Rinit + 6, Pbest + 8, 3 * sizeof(double));

	dgerqf_driver(3, 3, Rinit, triangular, orthogonal);	    

	/* Check the parity along the diagonal */
	neg = 
	    (triangular[0] < 0.0) + 
	    (triangular[4] < 0.0) + 
	    (triangular[8] < 0.0);

	if ((neg % 2) == 1) {
	    sign = -1;
	} else {
	    sign = 1;
	}
#endif

	num_inliers = 0;
	pts_final = (v3_t *) malloc(sizeof(v3_t) * max_inliers);
	projs_final = (v2_t *) malloc(sizeof(v2_t) * max_inliers);
	
	for (i = 0; i < num_pts; i++) {
	    double pt[4] = { Vx(points[i]), 
			     Vy(points[i]), 
			     Vz(points[i]), 1.0 };

	    double pr[3];
	    double dx, dy, dist;
	    
	    matrix_product341(Pbest, pt, pr);

	    /* Check cheirality */
	    if (sign * pr[2] > 0.0) 
		continue;

	    pr[0] /= -pr[2];
	    pr[1] /= -pr[2];
	    
	    dx = pr[0] - Vx(projs[i]);
	    dy = pr[1] - Vy(projs[i]);

	    dist = dx * dx + dy * dy;

	    if (dist < thresh_sq) {
		pts_final[num_inliers] = points[i];
		projs_final[num_inliers] = projs[i];
		num_inliers++;
	    }
	}

	if (num_inliers != max_inliers) {
	    printf("[find_projection_3x4_ransac] Error! There was a miscount "
		   "somewhere: (%d != %d)\n", num_inliers, max_inliers);
	}

	find_projection_3x4(max_inliers, pts_final, projs_final, Plinear);

#if 1
	/* Fix the sign on the P matrix */
	memcpy(Rinit + 0, Plinear + 0, 3 * sizeof(double));
	memcpy(Rinit + 3, Plinear + 4, 3 * sizeof(double));
	memcpy(Rinit + 6, Plinear + 8, 3 * sizeof(double));

	dgerqf_driver(3, 3, Rinit, triangular, orthogonal);	    
	
	/* Check the parity along the diagonal */
	neg = 
	    (triangular[0] < 0.0) + 
	    (triangular[4] < 0.0) + 
	    (triangular[8] < 0.0);

	if ((neg % 2) == 1) {
	    sign = -1;
	} else {
	    sign = 1;
	}
#endif

	for (i = 0; i < num_pts; i++) {
	    double pt[4] = 
		{ Vx(points[i]), Vy(points[i]), Vz(points[i]), 1.0 };
	    double pr[3];
	    double dx, dy, dist;
	    
	    matrix_product341(Plinear, pt, pr);

	    if (sign * pr[2] > 0.0)
		continue;

	    pr[0] /= -pr[2];
	    pr[1] /= -pr[2];
	    
	    dx = pr[0] - Vx(projs[i]);
	    dy = pr[1] - Vy(projs[i]);

	    dist = dx * dx + dy * dy;

	    if (dist < thresh_sq) {
		num_inliers_new++;
	    }
	}

	if (num_inliers_new < max_inliers) {
	    printf("[find_projection_3x4_ransac] Reverting to old solution\n");
	    memcpy(Plinear, Pbest, 12 * sizeof(double));
	}
	
	printf("Best matrix (pre-opt):\n");
	matrix_print(3, 4, Plinear);

	error = 0.0;
	for (i = 0; i < max_inliers; i++) {
	    double pt[4] = 
		{ Vx(pts_final[i]), Vy(pts_final[i]), Vz(pts_final[i]), 1.0 };
	    double pr[3];
	    double dx, dy, dist;
	    
	    matrix_product341(Plinear, pt, pr);
	    pr[0] /= pr[2];
	    pr[1] /= pr[2];
	    
	    dx = pr[0] - Vx(projs_final[i]);
	    dy = pr[1] - Vy(projs_final[i]);

	    dist = dx * dx + dy * dy;

	    error += dist;
	}
	
	printf("Old error: %0.3e\n", sqrt(error / max_inliers));

	/* Polish the result */
	if (max_inliers >= 6) {
	    int num_inliers_polished = 0;
	    find_projection_3x4_nonlinear(max_inliers, pts_final, projs_final,
					  Plinear, P);

#if 1
            /* Fix the sign on the P matrix */
            memcpy(Rinit + 0, P + 0, 3 * sizeof(double));
            memcpy(Rinit + 3, P + 4, 3 * sizeof(double));
            memcpy(Rinit + 6, P + 8, 3 * sizeof(double));

            dgerqf_driver(3, 3, Rinit, triangular, orthogonal);	    
            
            /* Check the parity along the diagonal */
            neg = 
                (triangular[0] < 0.0) + 
                (triangular[4] < 0.0) + 
                (triangular[8] < 0.0);

            if ((neg % 2) == 1) {
                sign = -1;
            } else {
                sign = 1;
            }
#endif

            /* Check that the number of inliers hasn't gone down */
	    num_inliers_polished = 0;
	    for (i = 0; i < num_pts; i++) {
		double pt[4] = 
		    { Vx(points[i]), Vy(points[i]), Vz(points[i]), 1.0 };
		double pr[3];
		double dx, dy, dist;
	    
		matrix_product341(P, pt, pr);

		if (sign * pr[2] > 0.0)
		    continue;

		pr[0] /= -pr[2];
		pr[1] /= -pr[2];
	    
		dx = pr[0] - Vx(projs[i]);
		dy = pr[1] - Vy(projs[i]);

		dist = dx * dx + dy * dy;

		if (dist < thresh_sq) {
		    num_inliers_polished++;
		}
	    }

	    if (num_inliers_polished < max_inliers) {
		printf("Decreased number of inliers (%d < %d), reverting\n",
		       num_inliers_polished, max_inliers);

		memcpy(P, Plinear, sizeof(double) * 12);		
	    }
	} else {
	    memcpy(P, Plinear, sizeof(double) * 12);
	}

	printf("Best matrix (post-opt):\n");
	matrix_print(3, 4, P);

	error = 0.0;
	for (i = 0; i < max_inliers; i++) {
	    double pt[4] = 
		{ Vx(pts_final[i]), Vy(pts_final[i]), Vz(pts_final[i]), 1.0 };
	    double pr[3];
	    double dx, dy, dist;
	    
	    matrix_product341(P, pt, pr);

	    pr[0] /= -pr[2];
	    pr[1] /= -pr[2];
	    
	    dx = pr[0] - Vx(projs_final[i]);
	    dy = pr[1] - Vy(projs_final[i]);

	    dist = dx * dx + dy * dy;

	    error += dist;
	}
	
	printf("New error: %0.3e\n", sqrt(error / max_inliers));

	free(inliers);
	free(pts_final);
	free(projs_final);

	return max_inliers;
    }
#undef MIN_PTS
}

