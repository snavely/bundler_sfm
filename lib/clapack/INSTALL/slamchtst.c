#include "f2c.h"
#include "blaswrap.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__4 = 4;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    real t, rnd, eps, base, emin, prec, emax, rmin, rmax, sfmin;
    extern doublereal slamch_(char *);

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };



/*  -- LAPACK test routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    eps = slamch_("Epsilon");
    sfmin = slamch_("Safe minimum");
    base = slamch_("Base");
    prec = slamch_("Precision");
    t = slamch_("Number of digits in mantissa");
    rnd = slamch_("Rounding mode");
    emin = slamch_("Minimum exponent");
    rmin = slamch_("Underflow threshold");
    emax = slamch_("Largest exponent");
    rmax = slamch_("Overflow threshold");

    s_wsle(&io___11);
    do_lio(&c__9, &c__1, " Epsilon                      = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&eps, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, " Safe minimum                 = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&sfmin, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___13);
    do_lio(&c__9, &c__1, " Base                         = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&base, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___14);
    do_lio(&c__9, &c__1, " Precision                    = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&prec, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___15);
    do_lio(&c__9, &c__1, " Number of digits in mantissa = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&t, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___16);
    do_lio(&c__9, &c__1, " Rounding mode                = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&rnd, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___17);
    do_lio(&c__9, &c__1, " Minimum exponent             = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&emin, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, " Underflow threshold          = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, " Largest exponent             = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&emax, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, " Overflow threshold           = ", (ftnlen)32);
    do_lio(&c__4, &c__1, (char *)&rmax, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, " Reciprocal of safe minimum   = ", (ftnlen)32);
    r__1 = 1 / sfmin;
    do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsle();

    return 0;
} /* MAIN__ */

/* Main program alias */ int test2_ () { MAIN__ (); return 0; }
