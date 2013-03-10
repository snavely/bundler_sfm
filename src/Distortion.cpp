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

/* Distortion.cpp */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Distortion.h"

#include "matrix.h"
#include "sfm.h"
#include "vector.h"

void InvertDistortion(int n_in, int n_out, double r0, double r1, 
                      double *k_in, double *k_out)
{
    const int NUM_SAMPLES = 20;
    
    int num_eqns = NUM_SAMPLES;
    int num_vars = n_out;
    
    double *A = new double[num_eqns * num_vars];
    double *b = new double[num_eqns];

    for (int i = 0; i < NUM_SAMPLES; i++) {
        double t = r0 + i * (r1 - r0) / (NUM_SAMPLES - 1);

        double p = 1.0;
        double a = 0.0;
        for (int j = 0; j < n_in; j++) {
            a += p * k_in[j];
            p = p * t;
        }

        double ap = 1.0;
        for (int j = 0; j < n_out; j++) {
            A[i * num_vars + j] = ap;
            ap = ap * a;
        }

        b[i] = t;
    }

    dgelsy_driver(A, b, k_out, num_eqns, num_vars, 1);

#if 0
    printf("[InvertDistortion] Output poly: ");
    matrix_print(1, n_out, k_out);

    for (int i = 0; i < NUM_SAMPLES; i++) {
        double t = r0 + i * (r1 - r0) / (NUM_SAMPLES - 1);

        double p = 1.0;
        double a = 0.0;
        for (int j = 0; j < n_in; j++) {
            a += p * k_in[j];
            p = p * t;
        }

        double ap = 1.0;
        double a2 = 0.0;
        for (int j = 0; j < n_out; j++) {
            a2 += ap * k_out[j];
            ap = ap * a;
        }

        printf("%0.3f -> %0.3f -> %0.3f\n", t, a, a2);
    }
#endif

    delete [] A;
    delete [] b;
}

v2_t UndistortNormalizedPoint(v2_t p, camera_params_t c) 
{
    double r = sqrt(Vx(p) * Vx(p) + Vy(p) * Vy(p));
    if (r == 0.0)
        return p;
    
    double t = 1.0;
    double a = 0.0;

    for (int i = 0; i < POLY_INVERSE_DEGREE; i++) {
        a += t * c.k_inv[i];
        t = t * r;
    }

    double factor = a / r;
    
    return v2_scale(factor, p);
}
