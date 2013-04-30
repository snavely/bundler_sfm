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

/*  svd.c -- Singular value decomposition. Translated to 'C' from the
 *           original Algol code in "Handbook for Automatic
 *           Computation,
 *           vol. II, Linear Algebra", Springer-Verlag.
 *
 *  (C) 2000, C. Bond. All rights reserved.
 *
 *  This is almost an exact translation from the original, except that
 *  an iteration counter is added to prevent stalls. This corresponds
 *  to similar changes in other translations.
 *
 *  Returns an error code = 0, if no errors and 'k' if a failure to
 *  converge at the 'kth' singular value.
 * 
 */
#ifdef __APPLE__
#include <malloc/malloc.h> /* for array allocation */
#else
#include <malloc.h> /* for array allocation */
#endif
#include <math.h>   /* for 'fabs'           */

#include "matrix.h"

int svd(int m,int n,int withu,int withv,double eps,double tol,
        double *a, double *q, double *u, double *v, double *vt)
{
    
    int i,j,k,l,l1,iter,retval;
    double c,f,g,h,s,x,y,z;
    double *e;
    
    e = (double *)calloc(n,sizeof(double));
    
    retval = 0;

    /* Copy 'a' to 'u' */    
    for (i=0;i<m;i++) {
        for (j=0;j<n;j++)
            u[i*n + j] = a[i*n + j];
    }
    
    /* Householder's reduction to bidiagonal form. */
    g = x = 0.0;
        
    for (i=0;i<n;i++) {
        
        e[i] = g;
        s = 0.0;
        l = i+1;
        
        for (j=i;j<m;j++)
            s += (u[j*n+i]*u[j*n+i]);
        
        if (s < tol)
            g = 0.0;
        
        else {
            
            f = u[i*n+i];
            
            g = (f < 0) ? sqrt(s) : -sqrt(s);
            
            h = f * g - s;
            
            u[i*n+i] = f - g;
            
            for (j=l;j<n;j++) {
                
                s = 0.0;
                
                for (k=i;k<m;k++)
                    s += (u[k*n+i] * u[k*n+j]);
                
                f = s / h;
                
                for (k=i;k<m;k++)
                    u[k*n+j] += (f * u[k*n+i]);
                
            }
            /* end j */
        }
        /* end s */
        q[i] = g;
        
        s = 0.0;
        
        for (j=l;j<n;j++)
            s += (u[i*n+j] * u[i*n+j]);
        
        if (s < tol)
            g = 0.0;
        
        else {
            
            f = u[i*n+i+1];
            
            g = (f < 0) ? sqrt(s) : -sqrt(s);
            
            h = f * g - s;
            
            u[i*n+i+1] = f - g;
            
            for (j=l;j<n;j++) 
                e[j] = u[i*n+j]/h;
            
            for (j=l;j<m;j++) {
                
                s = 0.0;
                
                for (k=l;k<n;k++) 
                    s += (u[j*n+k] * u[i*n+k]);
                
                for (k=l;k<n;k++)
                    u[j*n+k] += (s * e[k]);
                
            }
            /* end j */
        }
        /* end s */
        y = fabs(q[i]) + fabs(e[i]);
                                 
        if (y > x)
            x = y;
        
    }
    /* end i */

    /* accumulation of right-hand transformations */
    if (withv) {
        
        for (i=n-1;i>=0;i--) {
            
            if (g != 0.0) {
                
                h = u[i*n+i+1] * g;
                
                for (j=l;j<n;j++)
                    v[j*n+i] = u[i*n+j]/h;
                
                for (j=l;j<n;j++) {
                    
                    s = 0.0;
                    
                    for (k=l;k<n;k++) 
                        s += (u[i*n+k] * v[k*n+j]);
                    
                    for (k=l;k<n;k++)
                        v[k*n+j] += (s * v[k*n+i]);
                    
                    
                }
                /* end j */
            }
            /* end g */
            for (j=l;j<n;j++)
                v[i*n+j] = v[j*n+i] = 0.0;
            
            v[i*n+i] = 1.0;
            
            g = e[i];
            
            l = i;
            
        }
        /* end i */
        
    }
    /* end withv, parens added for clarity */

    /* accumulation of left-hand transformations */
    if (withu) {
        
        for (i=n;i<m;i++) {
            
            for (j=n;j<m;j++)
                u[i*n+j] = 0.0;
            
            u[i*n+i] = 1.0;
            
        }
        
    }
    
    if (withu) {
        
        for (i=n-1;i>=0;i--) {
            
            l = i + 1;
            
            g = q[i];
            
            for (j=l;j<m;j++)  /* upper limit was 'n' */
                u[i*n+j] = 0.0;
            
            if (g != 0.0) {
                
                h = u[i*n+i] * g;
                
                for (j=l;j<m;j++) {
                    /* upper limit was 'n' */
                    s = 0.0;
                    
                    for (k=l;k<m;k++)
                        s += (u[k*n+i] * u[k*n+j]);
                    
                    f = s / h;
                    
                    for (k=i;k<m;k++) 
                        u[k*n+j] += (f * u[k*n+i]);
                    
                }
                /* end j */
                for (j=i;j<m;j++) 
                    u[j*n+i] /= g;
                
            }
            /* end g */
            else {
                
                for (j=i;j<m;j++)
                    u[j*n+i] = 0.0;
                
            }
            
            u[i*n+i] += 1.0;
            
        }
        /* end i*/
    }
    /* end withu, parens added for clarity */

    /* diagonalization of the bidiagonal form */
    eps *= x;
    
    for (k=n-1;k>=0;k--) {
        
        iter = 0;
        
    test_f_splitting:
        for (l=k;l>=0;l--) {
            
            if (fabs(e[l]) <= eps) goto test_f_convergence;
            
            if (fabs(q[l-1]) <= eps) goto cancellation;
            
        }
        /* end l */

        /* cancellation of e[l] if l > 0 */
    cancellation:
        c = 0.0;
        
        s = 1.0;
        
        l1 = l - 1;
        
        for (i=l;i<=k;i++) {
            
            f = s * e[i];
            
            e[i] *= c;
            
            if (fabs(f) <= eps) goto test_f_convergence;
            
            g = q[i];
            
            h = q[i] = sqrt(f*f + g*g);
            
            c = g / h;
            
            s = -f / h;
            
            if (withu) {
                
                for (j=0;j<m;j++) {
                    
                    y = u[j*n+l1];
                    
                    z = u[j*n+i];
                    
                    u[j*n+l1] = y * c + z * s;
                    
                    u[j*n+i] = -y * s + z * c;
                    
                }
                /* end j */
            }
            /* end withu, parens added for clarity */
        }
        /* end i */
    test_f_convergence:
        z = q[k];
        
        if (l == k) goto convergence;
        

        /* shift from bottom 2x2 minor */
        iter++;
        
        if (iter > 30) {
            
            retval = k;
            
            break;
            
        }
        
        x = q[l];
        
        y = q[k-1];
        
        g = e[k-1];
        
        h = e[k];
        
        f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
        
        g = sqrt(f*f + 1.0);
        
        f = ((x-z)*(x+z) + h*(y/((f<0)?(f-g):(f+g))-h))/x;
        
        /* next QR transformation */
        c = s = 1.0;
        
        for (i=l+1;i<=k;i++) {
            
            g = e[i];
            
            y = q[i];
            
            h = s * g;
            
            g *= c;
            
            e[i-1] = z = sqrt(f*f+h*h);
            
            c = f / z;
            
            s = h / z;
            
            f = x * c + g * s;
            
            g = -x * s + g * c;
            
            h = y * s;
            
            y *= c;
            
            if (withv) {
                
                for (j=0;j<n;j++) {
                    
                    x = v[j*n+i-1];
                    
                    z = v[j*n+i];
                    
                    v[j*n+i-1] = x * c + z * s;
                    
                    v[j*n+i] = -x * s + z * c;
                    
                }
                /* end j */
            }
            /* end withv, parens added for clarity */
            q[i-1] = z = sqrt(f*f + h*h);
            
            c = f/z;
            
            s = h/z;
            
            f = c * g + s * y;
            
            x = -s * g + c * y;
            
            if (withu) {
                
                for (j=0;j<m;j++) {
                    
                    y = u[j*n+i-1];
                    
                    z = u[j*n+i];
                    
                    u[j*n+i-1] = y * c + z * s;
                    
                    u[j*n+i] = -y * s + z * c;
                    
                }
                /* end j */
            }
            /* end withu, parens added for clarity */
        }
        /* end i */
        e[l] = 0.0;
        
        e[k] = f;
        
        q[k] = x;
        
        goto test_f_splitting;
        
    convergence:
        if (z < 0.0) {
            
            /* q[k] is made non-negative */
            q[k] = - z;
            
            if (withv) {
                
                for (j=0;j<n;j++)
                    v[j*n+k] = -v[j*n+k];
                
            }
            /* end withv, parens added for clarity */
        }
        /* end z */
    }
    /* end k */
    
    free(e);

    matrix_transpose(m, n, v, vt);

    return retval;
    
}
