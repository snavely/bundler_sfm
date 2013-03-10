/* r1mpyq.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <f2c.h>

/* Subroutine */ int r1mpyq_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *v, doublereal *w)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nm1, nmj;
    static doublereal cos__, sin__, temp;

/*     ********** */

/*     subroutine r1mpyq */

/*     given an m by n matrix a, this subroutine computes a*q where */
/*     q is the product of 2*(n - 1) transformations */

/*           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) */

/*     and gv(i), gw(i) are givens rotations in the (i,n) plane which */
/*     eliminate elements in the i-th and n-th planes, respectively. */
/*     q itself is not given, rather the information to recover the */
/*     gv, gw rotations is supplied. */

/*     the subroutine statement is */

/*       subroutine r1mpyq(m,n,a,lda,v,w) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of a. */

/*       n is a positive integer input variable set to the number */
/*         of columns of a. */

/*       a is an m by n array. on input a must contain the matrix */
/*         to be postmultiplied by the orthogonal matrix q */
/*         described above. on output a*q has replaced a. */

/*       lda is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array a. */

/*       v is an input array of length n. v(i) must contain the */
/*         information necessary to recover the givens rotation gv(i) */
/*         described above. */

/*       w is an input array of length n. w(i) must contain the */
/*         information necessary to recover the givens rotation gw(i) */
/*         described above. */

/*     subroutines called */

/*       fortran-supplied ... dabs,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --w;
    --v;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */

/*     apply the first set of givens rotations to a. */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L50;
    }
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
	j = *n - nmj;
	if ((d__1 = v[j], abs(d__1)) > one) {
	    cos__ = one / v[j];
	}
	if ((d__1 = v[j], abs(d__1)) > one) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(one - d__2 * d__2);
	}
	if ((d__1 = v[j], abs(d__1)) <= one) {
	    sin__ = v[j];
	}
	if ((d__1 = v[j], abs(d__1)) <= one) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(one - d__2 * d__2);
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__ * a[i__ + j * a_dim1] - sin__ * a[i__ + *n * a_dim1];
	    a[i__ + *n * a_dim1] = sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + *n * a_dim1];
	    a[i__ + j * a_dim1] = temp;
/* L10: */
	}
/* L20: */
    }

/*     apply the second set of givens rotations to a. */

    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	if ((d__1 = w[j], abs(d__1)) > one) {
	    cos__ = one / w[j];
	}
	if ((d__1 = w[j], abs(d__1)) > one) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(one - d__2 * d__2);
	}
	if ((d__1 = w[j], abs(d__1)) <= one) {
	    sin__ = w[j];
	}
	if ((d__1 = w[j], abs(d__1)) <= one) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(one - d__2 * d__2);
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__ * a[i__ + j * a_dim1] + sin__ * a[i__ + *n * a_dim1];
	    a[i__ + *n * a_dim1] = -sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + *n * a_dim1];
	    a[i__ + j * a_dim1] = temp;
/* L30: */
	}
/* L40: */
    }
L50:
    return 0;

/*     last card of subroutine r1mpyq. */

} /* r1mpyq_ */

