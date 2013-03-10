/* fdjac1.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int fdjac1_(S_fp fcn, integer *n, doublereal *x, doublereal *
	fvec, doublereal *fjac, integer *ldfjac, integer *iflag, integer *ml, 
	integer *mu, doublereal *epsfcn, doublereal *wa1, doublereal *wa2)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal h__;
    static integer i__, j, k;
    static doublereal eps, temp;
    static integer msum;
    static doublereal epsmch;
    extern doublereal dpmpar_(integer *);

/*     ********** */

/*     subroutine fdjac1 */

/*     this subroutine computes a forward-difference approximation */
/*     to the n by n jacobian matrix associated with a specified */
/*     problem of n functions in n variables. if the jacobian has */
/*     a banded form, then function evaluations are saved by only */
/*     approximating the nonzero terms. */

/*     the subroutine statement is */

/*       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn, */
/*                         wa1,wa2) */

/*     where */

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions. fcn must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         subroutine fcn(n,x,fvec,iflag) */
/*         integer n,iflag */
/*         double precision x(n),fvec(n) */
/*         ---------- */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/*         ---------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of fdjac1. */
/*         in this case set iflag to a negative integer. */

/*       n is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       x is an input array of length n. */

/*       fvec is an input array of length n which must contain the */
/*         functions evaluated at x. */

/*       fjac is an output n by n array which contains the */
/*         approximation to the jacobian matrix evaluated at x. */

/*       ldfjac is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array fjac. */

/*       iflag is an integer variable which can be used to terminate */
/*         the execution of fdjac1. see description of fcn. */

/*       ml is a nonnegative integer input variable which specifies */
/*         the number of subdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         ml to at least n - 1. */

/*       epsfcn is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. this */
/*         approximation assumes that the relative errors in the */
/*         functions are of the order of epsfcn. if epsfcn is less */
/*         than the machine precision, it is assumed that the relative */
/*         errors in the functions are of the order of the machine */
/*         precision. */

/*       mu is a nonnegative integer input variable which specifies */
/*         the number of superdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         mu to at least n - 1. */

/*       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at */
/*         least n, then the jacobian is considered dense, and wa2 is */
/*         not referenced. */

/*     subprograms called */

/*       minpack-supplied ... dpmpar */

/*       fortran-supplied ... dabs,dmax1,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --fvec;
    --x;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = dpmpar_(&c__1);

    eps = sqrt((max(*epsfcn,epsmch)));
    msum = *ml + *mu + 1;
    if (msum < *n) {
	goto L40;
    }

/*        computation of dense approximate jacobian. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = x[j];
	h__ = eps * abs(temp);
	if (h__ == zero) {
	    h__ = eps;
	}
	x[j] = temp + h__;
	(*fcn)(n, &x[1], &wa1[1], iflag);
	if (*iflag < 0) {
	    goto L30;
	}
	x[j] = temp;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fjac[i__ + j * fjac_dim1] = (wa1[i__] - fvec[i__]) / h__;
/* L10: */
	}
/* L20: */
    }
L30:
    goto L110;
L40:

/*        computation of banded approximate jacobian. */

    i__1 = msum;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n;
	i__3 = msum;
	for (j = k; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {
	    wa2[j] = x[j];
	    h__ = eps * (d__1 = wa2[j], abs(d__1));
	    if (h__ == zero) {
		h__ = eps;
	    }
	    x[j] = wa2[j] + h__;
/* L60: */
	}
	(*fcn)(n, &x[1], &wa1[1], iflag);
	if (*iflag < 0) {
	    goto L100;
	}
	i__3 = *n;
	i__2 = msum;
	for (j = k; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
	    x[j] = wa2[j];
	    h__ = eps * (d__1 = wa2[j], abs(d__1));
	    if (h__ == zero) {
		h__ = eps;
	    }
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		fjac[i__ + j * fjac_dim1] = zero;
		if (i__ >= j - *mu && i__ <= j + *ml) {
		    fjac[i__ + j * fjac_dim1] = (wa1[i__] - fvec[i__]) / h__;
		}
/* L70: */
	    }
/* L80: */
	}
/* L90: */
    }
L100:
L110:
    return 0;

/*     last card of subroutine fdjac1. */

} /* fdjac1_ */

