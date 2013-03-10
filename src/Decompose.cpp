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

/* Decompose.cpp */
/* Decompose a homography into R and T */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BundleAdd.h"

#include "defines.h"
#include "matrix.h"

/* H is a map from image 2 to image 1 */
bool DecomposeHomography(double *H, double f1, double f2, 
                         double *Ra, double *ta, 
                         double *Rb, double *tb, v2_t p, v2_t q)
{
    double K2_inv[9] = { 1.0 / f2, 0.0,      0.0,
                         0.0,      1.0 / f2, 0.0,
                         0.0,      0.0,      1.0 };
    
    double K1[9] = { f1, 0.0, 0.0,
                     0.0, f1, 0.0,
                     0.0, 0.0, 1.0 };

    /* Reflect about y-axis */
#if 0
    double Yref[9] = { 1.0, 0.0, 0.0,
                       0.0, -1.0, 0.0,
                       0.0, 0.0, 1.0 };
#endif

    double Hnorm[9], tmp[9];
    matrix_product(3, 3, 3, 3, K2_inv, H, tmp);
    matrix_product(3, 3, 3, 3, tmp, K1, Hnorm);

    double U0[9], S0[3], VT0[9];
    dgesvd_driver(3, 3, Hnorm, U0, S0, VT0);
    matrix_scale(3, 3, Hnorm, 1.0 / S0[1], Hnorm);

    double HTH[9];
    matrix_transpose_product(3, 3, 3, 3, Hnorm, Hnorm, HTH);

    double U[9], S[3], VT[9];
    dgesvd_driver(3, 3, HTH, U, S, VT);
    S[1] = 1.0;

    double *v1 = VT + 0, *v2 = VT + 3, *v3 = VT + 6;

    double scale1 = sqrt(S[0] - 1.0);
    double scale3 = sqrt(1.0 - S[2]);
    
    double scale13 = sqrt(S[0] - S[2]);

    double v1s[3], v3s[3];
    matrix_scale(3, 1, v1, scale3, v1s);
    matrix_scale(3, 1, v3, scale1, v3s);
    
    double u1[3], u2[3];
    matrix_sum(3, 1, 3, 1, v1s, v3s, u1);
    matrix_diff(3, 1, 3, 1, v1s, v3s, u2);
    
    matrix_scale(3, 1, u1, 1.0 / scale13, u1);
    matrix_scale(3, 1, u2, 1.0 / scale13, u2);

    double Hu1[3], Hu2[3], Hv2[3];
    double v2_cross[9], Hv2_cross[9];

    matrix_product(3, 3, 3, 1, Hnorm, u1, Hu1);
    matrix_product(3, 3, 3, 1, Hnorm, u2, Hu2);
    matrix_product(3, 3, 3, 1, Hnorm, v2, Hv2);

    matrix_cross_matrix(v2, v2_cross);
    matrix_cross_matrix(Hv2, Hv2_cross);

    double U1T[9], U2T[9];
    double W1T[9], W2T[9];    
    double v2u1[3], v2u2[3], Hv2Hu1[3], Hv2Hu2[3];

    matrix_product(3, 3, 3, 1, v2_cross, u1, v2u1);
    matrix_product(3, 3, 3, 1, v2_cross, u2, v2u2);
    matrix_product(3, 3, 3, 1, Hv2_cross, Hu1, Hv2Hu1);
    matrix_product(3, 3, 3, 1, Hv2_cross, Hu2, Hv2Hu2);    

    memcpy(U1T + 0, v2, 3 * sizeof(double));
    memcpy(U1T + 3, u1, 3 * sizeof(double));
    memcpy(U1T + 6, v2u1, 3 * sizeof(double));
    
    memcpy(U2T + 0, v2, 3 * sizeof(double));
    memcpy(U2T + 3, u2, 3 * sizeof(double));
    memcpy(U2T + 6, v2u2, 3 * sizeof(double));

    memcpy(W1T + 0, Hv2, 3 * sizeof(double));
    memcpy(W1T + 3, Hu1, 3 * sizeof(double));
    memcpy(W1T + 6, Hv2Hu1, 3 * sizeof(double));

    memcpy(W2T + 0, Hv2, 3 * sizeof(double));
    memcpy(W2T + 3, Hu2, 3 * sizeof(double));
    memcpy(W2T + 6, Hv2Hu2, 3 * sizeof(double));

    double U1[9], U2[9], W1[9], W2[9];
    matrix_transpose(3, 3, U1T, U1);
    matrix_transpose(3, 3, U2T, U2);
    matrix_transpose(3, 3, W1T, W1);
    matrix_transpose(3, 3, W2T, W2);

    double R1[9], R2[9], R3[9], R4[9];
    double n1[3], n2[3], n3[3], n4[3];
    double t1[3], t2[3], t3[3], t4[3];

    /* First set of parameters */
    matrix_transpose_product2(3, 3, 3, 3, W1, U1, R1);
    memcpy(n1, v2u1, sizeof(double) * 3);
    matrix_diff(3, 3, 3, 3, Hnorm, R1, tmp);
    matrix_product(3, 3, 3, 1, tmp, n1, t1);

    /* Second set of parameters */
    matrix_transpose_product2(3, 3, 3, 3, W2, U2, R2);
    memcpy(n2, v2u2, sizeof(double) * 3);
    matrix_diff(3, 3, 3, 3, Hnorm, R2, tmp);
    matrix_product(3, 3, 3, 1, tmp, n2, t2);
    
    /* Third set of parameters */
    memcpy(R3, R1, sizeof(double) * 9);
    matrix_scale(3, 1, n1, -1.0, n3);
    matrix_scale(3, 1, t1, -1.0, t3);
    
    /* Fourth set of parameters */
    memcpy(R4, R2, sizeof(double) * 9);
    matrix_scale(3, 1, n2, -1.0, n4);
    matrix_scale(3, 1, t2, -1.0, t4);

    /* Choose the parameters by triangulating a point */
    camera_params_t camera0, camera1;
    matrix_ident(3, camera0.R);
    matrix_zeroes(3, 1, camera0.t);
    camera0.f = f1;

    camera1.f = f2;

    double error1, angle1;
    bool in_front1;
    memcpy(camera1.R, R1, sizeof(double) * 9);
    // memcpy(camera1.t, t1, sizeof(double) * 3);
    // Triangulate(p, q, camera0, camera1, error1, in_front1, angle1, false);
    matrix_transpose_product(3, 3, 3, 1, R1, t1, camera1.t);
    matrix_scale(3, 1, camera1.t, -1.0, camera1.t);
    Triangulate(p, q, camera0, camera1, error1, in_front1, angle1, true);

    double error2, angle2;
    bool in_front2;
    memcpy(camera1.R, R2, sizeof(double) * 9);
    // memcpy(camera1.t, t2, sizeof(double) * 3);
    // Triangulate(p, q, camera0, camera1, error2, in_front2, angle2, false);
    matrix_transpose_product(3, 3, 3, 1, R2, t2, camera1.t);
    matrix_scale(3, 1, camera1.t, -1.0, camera1.t);
    Triangulate(p, q, camera0, camera1, error2, in_front2, angle2, true);

    double error3, angle3;
    bool in_front3;
    memcpy(camera1.R, R3, sizeof(double) * 9);
    // memcpy(camera1.t, t3, sizeof(double) * 3);
    // Triangulate(p, q, camera0, camera1, error3, in_front3, angle3, false);
    matrix_transpose_product(3, 3, 3, 1, R3, t3, camera1.t);
    matrix_scale(3, 1, camera1.t, -1.0, camera1.t);
    Triangulate(p, q, camera0, camera1, error3, in_front3, angle3, true);

    double error4, angle4;
    bool in_front4;
    memcpy(camera1.R, R4, sizeof(double) * 9);
    // memcpy(camera1.t, t4, sizeof(double) * 3);
    // Triangulate(p, q, camera0, camera1, error4, in_front4, angle4, false);
    matrix_transpose_product(3, 3, 3, 1, R4, t4, camera1.t);
    matrix_scale(3, 1, camera1.t, -1.0, camera1.t);
    Triangulate(p, q, camera0, camera1, error4, in_front4, angle4, true);

    matrix_print(3, 3, R1);
    printf("\n");
    matrix_print(1, 3, t1);
    matrix_print(1, 3, n1);
    printf("error1: %0.3f   angle1: %0.3f  in-front1: %d\n",
           error1, RAD2DEG(angle1), in_front1 ? 1 : 0);

    matrix_print(3, 3, R2);
    printf("\n");
    matrix_print(1, 3, t2);
    matrix_print(1, 3, n2);
    printf("error2: %0.3f   angle2: %0.3f  in-front2: %d\n",
           error2, RAD2DEG(angle2), in_front2 ? 1 : 0);

    matrix_print(3, 3, R3);
    printf("\n");
    matrix_print(1, 3, t3);
    matrix_print(1, 3, n3);
    printf("error3: %0.3f   angle3: %0.3f  in-front3: %d\n",
           error3, RAD2DEG(angle3), in_front3 ? 1 : 0);

    matrix_print(3, 3, R4);
    printf("\n");
    matrix_print(1, 3, t4);
    matrix_print(1, 3, n4);
    printf("error4: %0.3f   angle4: %0.3f  in-front4: %d\n",
           error4, RAD2DEG(angle4), in_front4 ? 1 : 0);


    double max_plane_norm = 0.0;
    double error_best = 0.0;
    int best = -1, good1 = -1, good2 = -1;
    if (in_front1 && fabs(n1[2]) > max_plane_norm) {
        max_plane_norm = fabs(n1[2]);
        best = 1;
        error_best = error1;
    } 

    if (in_front1) {
        good1 = 1;
    }

    if (in_front2 && fabs(n2[2]) > max_plane_norm) {
        max_plane_norm = fabs(n2[2]);
        best = 2;
        error_best = error2;
    }

    if (in_front2) {
        if (good1 == -1)
            good1 = 2;
        else
            good2 = 2;        
    }
    

    if (in_front3 && fabs(n3[2]) > max_plane_norm) {
        max_plane_norm = fabs(n3[2]);
        best = 3;
        error_best = error3;
    }

    if (in_front3) {
        if (good1 == -1)
            good1 = 3;
        else
            good2 = 3;
        
    }

    if (in_front4 && fabs(n4[2]) > max_plane_norm) {
        max_plane_norm = fabs(n4[2]);
        best = 4;
        error_best = error4;
    }

    if (in_front4) {
        if (good1 == -1)
            good1 = 4;
        else
            good2 = 4;        
    }

    if (best == -1) {
        printf("[DecomposeHomography] No solution found\n");
        return false;
    }

    if (Ra != NULL && ta != NULL) {
        switch (best) {
        case 1:
            memcpy(Ra, R1, sizeof(double) * 9);
            memcpy(ta, t1, sizeof(double) * 3);
            break;
        case 2:
            memcpy(Ra, R2, sizeof(double) * 9);
            memcpy(ta, t2, sizeof(double) * 3);
            break;
        case 3:
            memcpy(Ra, R3, sizeof(double) * 9);
            memcpy(ta, t3, sizeof(double) * 3);
            break;
        case 4:
            memcpy(Ra, R4, sizeof(double) * 9);
            memcpy(ta, t4, sizeof(double) * 3);
            break;
        }

        if (Rb != NULL && tb != NULL) {
            int other = -1;
            if (good1 == best) {
                other = good2;
            } else {
                assert(good2 == best);
                other = good1;
            }

            if (other == -1) {
                printf("[DecomposeHomography] Error: other = -1\n");
                matrix_ident(3, Rb);
                tb[0] = tb[1] = tb[2];
            } else {
                switch (other) {
                case 1:
                    memcpy(Rb, R1, sizeof(double) * 9);
                    memcpy(tb, t1, sizeof(double) * 3);
                    break;
                case 2:
                    memcpy(Rb, R2, sizeof(double) * 9);
                    memcpy(tb, t2, sizeof(double) * 3);
                    break;
                case 3:
                    memcpy(Rb, R3, sizeof(double) * 9);
                    memcpy(tb, t3, sizeof(double) * 3);
                    break;
                case 4:
                    memcpy(Rb, R4, sizeof(double) * 9);
                    memcpy(tb, t4, sizeof(double) * 3);
                    break;
                }
            }
        }
    }
    
    if (error_best > 6.0) {
        printf("[DecomposeHomograph] Error [%0.3f] too high!\n", error_best);
        return false;
    }

    return true;
}

bool ComputeFundamentalMatrix(double f1, double f2, 
                              double *R2, double *t2, double *F)
{
    double t2_cross[9] = { 0.0, -t2[2], t2[1],
                           t2[2], 0.0, -t2[0],
                           -t2[1], t2[0], 0.0 };

    double E[9];
    matrix_product(3, 3, 3, 3, t2_cross, R2, E);

    /* Compute F matrix */
    double K1[9] = { f1, 0, 0,
                     0, f1, 0,
                     0, 0, 1 };
    double K2[9] = {f2, 0, 0,
                    0, f2, 0,
                    0, 0, 1 };
    
    double K1inv[9], K2inv[9];
    matrix_invert(3, K1, K1inv);
    matrix_invert(3, K2, K2inv);

    double tmp[9];
    matrix_transpose_product(3, 3, 3, 3, K2inv, E, tmp);
    matrix_product(3, 3, 3, 3, tmp, K1inv, F);

    return true;
}
