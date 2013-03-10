/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
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

/* LoadJPEG.h */
/* Code for reading a jpeg file */

#ifndef __load_jpeg_h__
#define __load_jpeg_h__

img_t *LoadJPEG(const char *filename);
void GetJPEGDimensions(const char *filename, int &w, int &h);
void WriteJPEG(const img_t *img, const char *filename);

#endif /* __load_jpeg_h__ */
