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

/* matrix.h */
/* Various linear algebra routines */

#ifndef __matrix_h__
#define __matrix_h__

#ifdef __cplusplus
extern "C" {
#endif

/* Fill a given matrix with an n x n identity matrix */
void matrix_ident(int n, double *A);

/* Fill a given matrix with an m x n matrix of zeroes */
void matrix_zeroes(int m, int n, double *A);
    
/* Transpose the m x n matrix A and put the result in the n x m matrix AT */
void matrix_transpose(int m, int n, double *A, double *AT);

/* Compute the matrix product R = AB */
void matrix_product(int Am, int An, int Bm, int Bn, 
                    const double *A, const double *B, double *R);

void matrix_product_old(int Am, int An, int Bm, int Bn, 
                        const double *A, const double *B, double *R);
void matrix_transpose_product_old(int Am, int An, int Bm, int Bn, 
                                  double *A, double *B, double *R);
void matrix_transpose_product2_old(int Am, int An, int Bm, int Bn, 
                                   double *A, double *B, double *R);

void matrix_product_ipp(int Am, int An, int Bn, 
                        const double *A, const double *B, double *R);
void matrix_transpose_product_ipp(int Am, int An, int Bn, 
                                  const double *A, const double *B, double *R);
void matrix_transpose_product2_ipp(int Am, int An, int Bm, 
                                   const double *A, const double *B, 
                                   double *R);
void matrix_array_product_ipp(int count, int Am, int An, int Bn,
                              const double *A, const double *B, double *R);

static inline void matrix_product33(double *A, double *B, double *R)
{
    R[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    R[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    R[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];    
    
    R[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    R[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    R[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];    

    R[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    R[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    R[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];            
}       

static inline void matrix_product121(double *A, double *b, double *r)
{
    r[0] = A[0] * b[0] + A[1] * b[1];            
}
    
static inline void matrix_product131(double *A, double *b, double *r)
{
    r[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
}
       
static inline void matrix_product331(double *A, double *b, double *r)
{
    r[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
    r[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
    r[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}

void matrix_product341(double *A, double *b, double *r);    
void matrix_product44(double *A, double *B, double *R);
void matrix_product441(double *A, double *b, double *r);
    
/* Compute the power of a matrix */
void matrix_power(int n, double *A, int pow, double *R);

/* Compute the matrix product R = A B^T */
void matrix_transpose_product2(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

/* Compute the matrix sum R = A + B */
void matrix_sum(int Am, int An, int Bm, int Bn, 
                double *A, double *B, double *R);

/* Compute the matrix difference R = A - B */
void matrix_diff(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

/* Compute the determinant of a 3x3 matrix */
double matrix_determinant3(double *A);

/* Compute the matrix product R = A^T B */
void matrix_transpose_product(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

/* Compute (transpose of) LU decomposition of A */
void matrix_lu(int n, double *A, double *LU, int *ipiv);
void matrix_lu_no_transpose(int n, double *A, double *LU, int *ipiv);
    
/* Solve a system of equations using a precomputed LU decomposition */
void matrix_solve_lu(int n, double *LU, int *ipiv, double *b, double *x);

/* Invert the n-by-n matrix A, storing the result in Ainv */
void matrix_invert(int n, double *A, double *Ainv);
void matrix_invert_inplace(int n, double *A);
    
/* Get the norm of the matrix */
double matrix_norm(int m, int n, double *A);

/* Get the [squared] norm of the matrix */
double matrix_normsq(int m, int n, double *A);

/* Scale a matrix by a scalar */
void matrix_scale(int m, int n, double *A, double s, double *R);

/* Print the given m x n matrix */
void matrix_print(int m, int n, double *A);

/* Read a matrix from a file */
void matrix_read_file(int m, int n, double *matrix, char *fname);

/* Write a matrix to a file */
void matrix_write_file(int m, int n, double *matrix, char *fname);

/* Return the product x**T A x */
double matrix_double_product(int n, double *A, double *x);

/* Compute the cross product of two 3 x 1 vectors */
void matrix_cross(const double *u, const double *v, double *w);
void matrix_cross4(const double *u, const double *v, const double *w, 
		   double *x);
    
/* Create the 3x3 cross product matrix from a 3-vector */
void matrix_cross_matrix(double *v, double *v_cross);

/* Convert a rotation matrix to axis and angle representation */
void matrix_to_axis_angle(double *R, double *axis, double *angle);
void axis_angle_to_matrix(double *axis, double angle, double *R);
void axis_angle_to_matrix4(double *axis, double angle, double *R);

/* Convert a matrix to a normalize quaternion */
void matrix_to_quaternion(double *R, double *q);
/* Convert a normalized quaternion to a matrix */
void quaternion_to_matrix(double *q, double *R);
    
/* Decompose a square matrix into an orthogonal matrix and a symmetric
 * positive semidefinite matrix */
void matrix_polar_decomposition(int n, double *A, double *Q, double *S);
    
/* Driver for the minpack function lmdif, which uses
 * Levenberg-Marquardt for non-linear least squares minimization */
void lmdif_driver(void *fcn, int m, int n, double *xvec, double tol);
void lmdif_driver2(void *fcn, int m, int n, double *xvec, double tol);
void lmdif_driver3(void *fcn, int m, int n, double *xvec, double tol,
                   int maxfev, double *H);
    
/* Driver for the lapack function dgelss, which finds x to minimize
 * norm(b - A * x) */
void dgelss_driver(double *A, double *b, double *x, int m, int n, int nrhs);
void dgelsy_driver(double *A, double *b, double *x, int m, int n, int nrhs);

/* Version of above where matrix is already in column-major order */
void dgelsy_driver_transpose(double *A, double *b, double *x, 
			     int m, int n, int nrhs);

/* Solve an n x n system */
void dgesv_driver(int n, double *A, double *b, double *x);
    
/* n: the order of matrix A
 * A: matrix for which the eigenvectors/values are to be computed
 * evec: output array containing the eigenvectors
 * eval: output array containing the eigenvalues
 *
 * Note: Assumes the results are real! */
int dgeev_driver(int n, double *A, double *evec, double *eval);

/* Compute singular value decomposition of an m x n matrix A */
int dgesvd_driver(int m, int n, double *A, double *U, double *S, double *VT);
/* Compute singular value decomposition of an m x n matrix A 
 * (only compute S and VT) */
int dgesvd_driver_vt(int m, int n, double *A, double *S, double *VT);

/* Compute Cholesky decomposition of an nxn matrix */
void dpotrf_driver(int n, double *A, double *U);

/* Compute a QR factorization of an m by n matrix A */
void dgeqrf_driver(int m, int n, double *A, double *Q, double *R);

/* Compute an RQ factorization of an m by n matrix A */
void dgerqf_driver(int m, int n, double *A, double *R, double *Q);

/* Find the unit vector that minimizes ||Ax|| */
void matrix_minimum_unit_norm_solution(int m, int n, double *A, double *x);

void cblas_dgemm_driver_x(int m1, int n12, int k2,
                          int as, int bs, int cs,
                          double *a, double *b, double *c);

void cblas_dgemm_driver(int m1, int n12, int k2, 
                        double *a, double *b, double *c);
    
void cblas_dgemm_driver_transpose(int m1, int n12, int k2, 
                                  double *a, double *b, double *c);
    
void cblas_dgemm_driver_transpose2(int m1, int n12, int k2, 
                                   double *a, double *b, double *c);

void slerp(double *v1, double *v2, double t, double *v3);
    
#ifdef __cplusplus
}
#endif

#endif /* __matrix_h__ */
