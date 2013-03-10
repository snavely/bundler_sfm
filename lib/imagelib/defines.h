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

/* defines.h */
/* Contains standard definitions */

#ifndef __defines_h__
#define __defines_h__

#define CLAMP(x,mn,mx) (((x) < mn) ? mn : (((x) > mx) ? mx : (x)))

#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

#define MAX4(a,b,c,d) (MAX(MAX(a,b), MAX(c,d)))
#define MIN4(a,b,c,d) (MIN(MIN(a,b), MIN(c,d)))

#define DEG2RAD(d) ((d) * (M_PI / 180.0))
#define RAD2DEG(r) ((r) * (180.0 / M_PI))

#define SQ(x) ((x) * (x))
#define CB(x) ((x) * (x) * (x))

/* Return the index of entry (x, y) in a 2D array */
#define INDEX(x, y, w) ((y) * (w) + (x))

#define SGN(x) ((x) < 0 ? (-1) : (1))

#ifdef WIN32
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif
#endif

#endif /* __defines_h__ */
