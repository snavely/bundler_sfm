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

/* pad.h */
/* Contains definitions for padding integers in various ways */

#ifndef __pad_h__
#define __pad_h__

#define BIT_PAD_BYTE(x) ((((x) % 8) == 0) ? (x) : (x) + 8 - ((x) % 8))
#define BIT_PAD_WORD(x) ((((x) % 32) == 0) ? (x) : (x) + 32 - ((x) % 32))
#define BYTE_PAD_WORD(x) ((((x) % 4) == 0) ? (x) : (x) + 4 - ((x) % 4))

#endif /* __pad_h__ */
