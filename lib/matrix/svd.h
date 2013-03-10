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

#ifndef __svd_h__
#define __svd_h__

#ifdef __cplusplus
extern "C" {
#endif

int svd(int m,int n,int withu,int withv,double eps,double tol,
        double *a, double *q, double *u, double *v, double *vt);

#ifdef __cplusplus
}
#endif

#endif /* __matrix_h__ */
