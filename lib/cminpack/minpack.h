/* minpack.h
 *
 * This is a C version of the minpack minimization package.
 * It has been derived from the fortran code using f2c.
 * See http://www.netlib.org/minpack/ for more information.
 * Minpack was developed by Jorge More', Burt Garbow, and 
 * Ken Hillstrom at Argonne National Laboratory
 * Converted to C by Manolis Lourakis, July 2002
 *
 */

#ifndef __MINPACK_H
#define __MINPACK_H

#define __CMINPACK_VERSION "1.0"

#include <f2c.h>

#ifdef __cplusplus
extern "C" {
#endif

/* automatically generated declarations using cproto */

/* chkder.c */
int chkder_(integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *xp, doublereal *fvecp, integer *mode, doublereal *err);

/* dogleg.c */
int dogleg_(integer *n, doublereal *r__, integer *lr, doublereal *diag, doublereal *qtb, doublereal *delta, doublereal *x, doublereal *wa1, doublereal *wa2);

/* dpmpar.c */
doublereal dpmpar_(integer *i__);

/* enorm.c */
doublereal enorm_(integer *n, doublereal *x);

/* fdjac1.c */
int fdjac1_(S_fp fcn, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, integer *iflag, integer *ml, integer *mu, doublereal *epsfcn, doublereal *wa1, doublereal *wa2);

/* fdjac2.c */
int fdjac2_(S_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, integer *iflag, doublereal *epsfcn, doublereal *wa);

/* hybrd1.c */
int hybrd1_(U_fp fcn, integer *n, doublereal *x, doublereal *fvec, doublereal *tol, integer *info, doublereal *wa, integer *lwa);

/* hybrd.c */
int hybrd_(S_fp fcn, integer *n, doublereal *x, doublereal *fvec, doublereal *xtol, integer *maxfev, integer *ml, integer *mu, doublereal *epsfcn, doublereal *diag, integer *mode, doublereal *factor, integer *nprint, integer *info, integer *nfev, doublereal *fjac, integer *ldfjac, doublereal *r__, integer *lr, doublereal *qtf, doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *wa4);
/* hybrj1.c */

int hybrj1_(U_fp fcn, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *tol, integer *info, doublereal *wa, integer *lwa);

/* hybrj.c */
int hybrj_(S_fp fcn, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *xtol, integer *maxfev, doublereal *diag, integer *mode, doublereal *factor, integer *nprint, integer *info, integer *nfev, integer *njev, doublereal *r__, integer *lr, doublereal *qtf, doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *wa4);

/* lmder1.c */
int lmder1_(U_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *tol, integer *info, integer *ipvt, doublereal *wa, integer *lwa);

/* lmder.c */
int lmder_(S_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *ftol, doublereal *xtol, doublereal *gtol, integer *maxfev, doublereal *diag, integer *mode, doublereal *factor, integer *nprint, integer *info, integer *nfev, integer *njev, integer *ipvt, doublereal *qtf, doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *wa4);

/* lmdif1.c */
int lmdif1_(U_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *tol, integer *info, integer *iwa, doublereal *wa, integer *lwa);

/* lmdif.c */
int lmdif_(S_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *ftol, doublereal *xtol, doublereal *gtol, integer *maxfev, doublereal *epsfcn, doublereal *diag, integer *mode, doublereal *factor, integer *nprint, integer *info, integer *nfev, doublereal *fjac, integer *ldfjac, integer *ipvt, doublereal *qtf, doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *wa4);

/* lmpar.c */
int lmpar_(integer *n, doublereal *r__, integer *ldr, integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *delta, doublereal *par, doublereal *x, doublereal *sdiag, doublereal *wa1, doublereal *wa2);

/* lmstr1.c */
int lmstr1_(U_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *tol, integer *info, integer *ipvt, doublereal *wa, integer *lwa);

/* lmstr.c */
int lmstr_(S_fp fcn, integer *m, integer *n, doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *ftol, doublereal *xtol, doublereal *gtol, integer *maxfev, doublereal *diag, integer *mode, doublereal *factor, integer *nprint, integer *info, integer *nfev, integer *njev, integer *ipvt, doublereal *qtf, doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *wa4);

/* qform.c */
int qform_(integer *m, integer *n, doublereal *q, integer *ldq, doublereal *wa);

/* qrfac.c */
int qrfac_(integer *m, integer *n, doublereal *a, integer *lda, logical *pivot, integer *ipvt, integer *lipvt, doublereal *rdiag, doublereal *acnorm, doublereal *wa);

/* qrsolv.c */
int qrsolv_(integer *n, doublereal *r__, integer *ldr, integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *x, doublereal *sdiag, doublereal *wa);

/* r1mpyq.c */
int r1mpyq_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *w);

/* r1updt.c */
int r1updt_(integer *m, integer *n, doublereal *s, integer *ls, doublereal *u, doublereal *v, doublereal *w, logical *sing);

/* rwupdt.c */
int rwupdt_(integer *n, doublereal *r__, integer *ldr, doublereal *w, doublereal *b, doublereal *alpha, doublereal *cos__, doublereal *sin__);

#ifdef __cplusplus
}; /* end of extern "C" */
#endif

#endif /* __MINPACK_H */
