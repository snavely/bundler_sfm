#include "f2c.h"
#include "blaswrap.h"

integer icamax_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    integer i__, ix;
    real smax;
    extern doublereal scabs1_(complex *);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    smax = scabs1_(&cx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (scabs1_(&cx[ix]) <= smax) {
	    goto L5;
	}
	ret_val = i__;
	smax = scabs1_(&cx[ix]);
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    smax = scabs1_(&cx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (scabs1_(&cx[i__]) <= smax) {
	    goto L30;
	}
	ret_val = i__;
	smax = scabs1_(&cx[i__]);
L30:
	;
    }
    return ret_val;
} /* icamax_ */
