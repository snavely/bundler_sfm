/* qform.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <f2c.h>

/* Subroutine */ int qform_(integer *m, integer *n, doublereal *q, integer *
	ldq, doublereal *wa)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, jm1, np1;
    static doublereal sum, temp;
    static integer minmn;

/*     ********** */

/*     subroutine qform */

/*     this subroutine proceeds from the computed qr factorization of */
/*     an m by n matrix a to accumulate the m by m orthogonal matrix */
/*     q from its factored form. */

/*     the subroutine statement is */

/*       subroutine qform(m,n,q,ldq,wa) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of a and the order of q. */

/*       n is a positive integer input variable set to the number */
/*         of columns of a. */

/*       q is an m by m array. on input the full lower trapezoid in */
/*         the first min(m,n) columns of q contains the factored form. */
/*         on output q has been accumulated into a square matrix. */

/*       ldq is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array q. */

/*       wa is a work array of length m. */

/*     subprograms called */

/*       fortran-supplied ... min0 */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;

    /* Function Body */

/*     zero out upper triangle of q in the first min(m,n) columns. */

    minmn = min(*m,*n);
    if (minmn < 2) {
	goto L30;
    }
    i__1 = minmn;
    for (j = 2; j <= i__1; ++j) {
	jm1 = j - 1;
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + j * q_dim1] = zero;
/* L10: */
	}
/* L20: */
    }
L30:

/*     initialize remaining columns to those of the identity matrix. */

    np1 = *n + 1;
    if (*m < np1) {
	goto L60;
    }
    i__1 = *m;
    for (j = np1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + j * q_dim1] = zero;
/* L40: */
	}
	q[j + j * q_dim1] = one;
/* L50: */
    }
L60:

/*     accumulate q from its factored form. */

    i__1 = minmn;
    for (l = 1; l <= i__1; ++l) {
	k = minmn - l + 1;
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    wa[i__] = q[i__ + k * q_dim1];
	    q[i__ + k * q_dim1] = zero;
/* L70: */
	}
	q[k + k * q_dim1] = one;
	if (wa[k] == zero) {
	    goto L110;
	}
	i__2 = *m;
	for (j = k; j <= i__2; ++j) {
	    sum = zero;
	    i__3 = *m;
	    for (i__ = k; i__ <= i__3; ++i__) {
		sum += q[i__ + j * q_dim1] * wa[i__];
/* L80: */
	    }
	    temp = sum / wa[k];
	    i__3 = *m;
	    for (i__ = k; i__ <= i__3; ++i__) {
		q[i__ + j * q_dim1] -= temp * wa[i__];
/* L90: */
	    }
/* L100: */
	}
L110:
/* L120: */
	;
    }
    return 0;

/*     last card of subroutine qform. */

} /* qform_ */

