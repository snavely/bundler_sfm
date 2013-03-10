/* dogleg.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dogleg_(integer *n, doublereal *r__, integer *lr, 
	doublereal *diag, doublereal *qtb, doublereal *delta, doublereal *x, 
	doublereal *wa1, doublereal *wa2)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, jj, jp1;
    static doublereal sum, temp, alpha, bnorm;
    extern doublereal enorm_(integer *, doublereal *);
    static doublereal gnorm, qnorm, epsmch;
    extern doublereal dpmpar_(integer *);
    static doublereal sgnorm;

/*     ********** */

/*     subroutine dogleg */

/*     given an m by n matrix a, an n by n nonsingular diagonal */
/*     matrix d, an m-vector b, and a positive number delta, the */
/*     problem is to determine the convex combination x of the */
/*     gauss-newton and scaled gradient directions that minimizes */
/*     (a*x - b) in the least squares sense, subject to the */
/*     restriction that the euclidean norm of d*x be at most delta. */

/*     this subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     qr factorization of a. that is, if a = q*r, where q has */
/*     orthogonal columns and r is an upper triangular matrix, */
/*     then dogleg expects the full upper triangle of r and */
/*     the first n components of (q transpose)*b. */

/*     the subroutine statement is */

/*       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2) */

/*     where */

/*       n is a positive integer input variable set to the order of r. */

/*       r is an input array of length lr which must contain the upper */
/*         triangular matrix r stored by rows. */

/*       lr is a positive integer input variable not less than */
/*         (n*(n+1))/2. */

/*       diag is an input array of length n which must contain the */
/*         diagonal elements of the matrix d. */

/*       qtb is an input array of length n which must contain the first */
/*         n elements of the vector (q transpose)*b. */

/*       delta is a positive input variable which specifies an upper */
/*         bound on the euclidean norm of d*x. */

/*       x is an output array of length n which contains the desired */
/*         convex combination of the gauss-newton direction and the */
/*         scaled gradient direction. */

/*       wa1 and wa2 are work arrays of length n. */

/*     subprograms called */

/*       minpack-supplied ... dpmpar,enorm */

/*       fortran-supplied ... dabs,dmax1,dmin1,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --x;
    --qtb;
    --diag;
    --r__;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = dpmpar_(&c__1);

/*     first, calculate the gauss-newton direction. */

    jj = *n * (*n + 1) / 2 + 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	j = *n - k + 1;
	jp1 = j + 1;
	jj -= k;
	l = jj + 1;
	sum = zero;
	if (*n < jp1) {
	    goto L20;
	}
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    sum += r__[l] * x[i__];
	    ++l;
/* L10: */
	}
L20:
	temp = r__[jj];
	if (temp != zero) {
	    goto L40;
	}
	l = j;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = temp, d__3 = (d__1 = r__[l], abs(d__1));
	    temp = max(d__2,d__3);
	    l = l + *n - i__;
/* L30: */
	}
	temp = epsmch * temp;
	if (temp == zero) {
	    temp = epsmch;
	}
L40:
	x[j] = (qtb[j] - sum) / temp;
/* L50: */
    }

/*     test whether the gauss-newton direction is acceptable. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = zero;
	wa2[j] = diag[j] * x[j];
/* L60: */
    }
    qnorm = enorm_(n, &wa2[1]);
    if (qnorm <= *delta) {
	goto L140;
    }

/*     the gauss-newton direction is not acceptable. */
/*     next, calculate the scaled gradient direction. */

    l = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = qtb[j];
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    wa1[i__] += r__[l] * temp;
	    ++l;
/* L70: */
	}
	wa1[j] /= diag[j];
/* L80: */
    }

/*     calculate the norm of the scaled gradient and test for */
/*     the special case in which the scaled gradient is zero. */

    gnorm = enorm_(n, &wa1[1]);
    sgnorm = zero;
    alpha = *delta / qnorm;
    if (gnorm == zero) {
	goto L120;
    }

/*     calculate the point along the scaled gradient */
/*     at which the quadratic is minimized. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = wa1[j] / gnorm / diag[j];
/* L90: */
    }
    l = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = zero;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum += r__[l] * wa1[i__];
	    ++l;
/* L100: */
	}
	wa2[j] = sum;
/* L110: */
    }
    temp = enorm_(n, &wa2[1]);
    sgnorm = gnorm / temp / temp;

/*     test whether the scaled gradient direction is acceptable. */

    alpha = zero;
    if (sgnorm >= *delta) {
	goto L120;
    }

/*     the scaled gradient direction is not acceptable. */
/*     finally, calculate the point along the dogleg */
/*     at which the quadratic is minimized. */

    bnorm = enorm_(n, &qtb[1]);
    temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / *delta);
/* Computing 2nd power */
    d__1 = sgnorm / *delta;
/* Computing 2nd power */
    d__2 = temp - *delta / qnorm;
/* Computing 2nd power */
    d__3 = *delta / qnorm;
/* Computing 2nd power */
    d__4 = sgnorm / *delta;
    temp = temp - *delta / qnorm * (d__1 * d__1) + sqrt(d__2 * d__2 + (one - 
	    d__3 * d__3) * (one - d__4 * d__4));
/* Computing 2nd power */
    d__1 = sgnorm / *delta;
    alpha = *delta / qnorm * (one - d__1 * d__1) / temp;
L120:

/*     form appropriate convex combination of the gauss-newton */
/*     direction and the scaled gradient direction. */

    temp = (one - alpha) * min(sgnorm,*delta);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = temp * wa1[j] + alpha * x[j];
/* L130: */
    }
L140:
    return 0;

/*     last card of subroutine dogleg. */

} /* dogleg_ */

