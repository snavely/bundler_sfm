#include "f2c.h"
#include "blaswrap.h"

doublereal dcabs1_(doublecomplex *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. */
/*  Purpose */
/*  ======= */

/*  DCABS1 computes absolute value of a double complex number */

/*     .. Intrinsic Functions .. */

    ret_val = (d__1 = z__->r, abs(d__1)) + (d__2 = d_imag(z__), abs(d__2));
    return ret_val;
} /* dcabs1_ */
