/* chkder.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int chkder_(integer *m, integer *n, doublereal *x, 
	doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *xp, 
	doublereal *fvecp, integer *mode, doublereal *err)
{
    /* Initialized data */

    static doublereal factor = 100.;
    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal), d_lg10(doublereal *);

    /* Local variables */
    static integer i__, j;
    static doublereal eps, epsf, temp, epsmch;
    extern doublereal dpmpar_(integer *);
    static doublereal epslog;

/*     ********** */

/*     subroutine chkder */

/*     this subroutine checks the gradients of m nonlinear functions */
/*     in n variables, evaluated at a point x, for consistency with */
/*     the functions themselves. the user must call chkder twice, */
/*     first with mode = 1 and then with mode = 2. */

/*     mode = 1. on input, x must contain the point of evaluation. */
/*               on output, xp is set to a neighboring point. */

/*     mode = 2. on input, fvec must contain the functions and the */
/*                         rows of fjac must contain the gradients */
/*                         of the respective functions each evaluated */
/*                         at x, and fvecp must contain the functions */
/*                         evaluated at xp. */
/*               on output, err contains measures of correctness of */
/*                          the respective gradients. */

/*     the subroutine does not perform reliably if cancellation or */
/*     rounding errors cause a severe loss of significance in the */
/*     evaluation of a function. therefore, none of the components */
/*     of x should be unusually small (in particular, zero) or any */
/*     other value which may cause loss of significance. */

/*     the subroutine statement is */

/*       subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of functions. */

/*       n is a positive integer input variable set to the number */
/*         of variables. */

/*       x is an input array of length n. */

/*       fvec is an array of length m. on input when mode = 2, */
/*         fvec must contain the functions evaluated at x. */

/*       fjac is an m by n array. on input when mode = 2, */
/*         the rows of fjac must contain the gradients of */
/*         the respective functions evaluated at x. */

/*       ldfjac is a positive integer input parameter not less than m */
/*         which specifies the leading dimension of the array fjac. */

/*       xp is an array of length n. on output when mode = 1, */
/*         xp is set to a neighboring point of x. */

/*       fvecp is an array of length m. on input when mode = 2, */
/*         fvecp must contain the functions evaluated at xp. */

/*       mode is an integer input variable set to 1 on the first call */
/*         and 2 on the second. other values of mode are equivalent */
/*         to mode = 1. */

/*       err is an array of length m. on output when mode = 2, */
/*         err contains measures of correctness of the respective */
/*         gradients. if there is no severe loss of significance, */
/*         then if err(i) is 1.0 the i-th gradient is correct, */
/*         while if err(i) is 0.0 the i-th gradient is incorrect. */
/*         for values of err between 0.0 and 1.0, the categorization */
/*         is less certain. in general, a value of err(i) greater */
/*         than 0.5 indicates that the i-th gradient is probably */
/*         correct, while a value of err(i) less than 0.5 indicates */
/*         that the i-th gradient is probably incorrect. */

/*     subprograms called */

/*       minpack supplied ... dpmpar */

/*       fortran supplied ... dabs,dlog10,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --err;
    --fvecp;
    --fvec;
    --xp;
    --x;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = dpmpar_(&c__1);

    eps = sqrt(epsmch);

    if (*mode == 2) {
	goto L20;
    }

/*        mode = 1. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = eps * (d__1 = x[j], abs(d__1));
	if (temp == zero) {
	    temp = eps;
	}
	xp[j] = x[j] + temp;
/* L10: */
    }
    goto L70;
L20:

/*        mode = 2. */

    epsf = factor * epsmch;
    epslog = d_lg10(&eps);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	err[i__] = zero;
/* L30: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = (d__1 = x[j], abs(d__1));
	if (temp == zero) {
	    temp = one;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    err[i__] += temp * fjac[i__ + j * fjac_dim1];
/* L40: */
	}
/* L50: */
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = one;
	if (fvec[i__] != zero && fvecp[i__] != zero && (d__2 = fvecp[i__] - 
		fvec[i__], abs(d__2)) >= epsf * (d__1 = fvec[i__], abs(d__1)))
		 {
	    temp = eps * (d__3 = (fvecp[i__] - fvec[i__]) / eps - err[i__], 
		    abs(d__3)) / ((d__4 = fvec[i__], abs(d__4)) + (d__5 = 
		    fvecp[i__], abs(d__5)));
	}
	err[i__] = one;
	if (temp > epsmch && temp < eps) {
	    err[i__] = (d_lg10(&temp) - epslog) / epslog;
	}
	if (temp >= eps) {
	    err[i__] = zero;
	}
/* L60: */
    }
L70:

    return 0;

/*     last card of subroutine chkder. */

} /* chkder_ */

