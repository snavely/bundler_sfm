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

/* pgm.c */
/* Read PGM files */

#include "image.h"

static img_t *img_read_pgm(FILE *fp);

/* PGM files allow a comment starting with '#' to end-of-line.  
 * Skip white space including any comments. */
void skip_comments(FILE *fp)
{
    int ch;

    fscanf(fp," ");      /* Skip white space. */
    while ((ch = fgetc(fp)) == '#') {
      while ((ch = fgetc(fp)) != '\n'  &&  ch != EOF)
	;
      fscanf(fp," ");
    }
    ungetc(ch, fp);      /* Replace last character read. */
}

/* Read a PGM file */
img_t *img_read_pgm_file(const char *filename)
{
    FILE *f;
    img_t *img;

    /* The "b" option is for binary input, which is needed if this is
     * compiled under Windows.  It has no effect in Linux. */

    f = fopen(filename, "rb");

    if (f == NULL) {
	printf("Error: could not open file %s", filename);
	return NULL;
    }

    img = img_read_pgm(f);
    fclose(f);

    return img;
}

/* Read the dimensions of an image */
int img_read_pgm_dimensions(char *filename, int *w, int *h) 
{
    FILE *fp = fopen(filename, "rb");
    int char1, char2, c1, c2;

    if (fp == NULL) {
	printf("Error: could not open file %s", filename);
	return -1;
    }

    char1 = fgetc(fp);
    char2 = fgetc(fp);
    skip_comments(fp);
    c1 = fscanf(fp, "%d", w);
    skip_comments(fp);
    c2 = fscanf(fp, "%d", h);
    
    fclose(fp);

    return 0;
}


/* Read a PGM file from the given file pointer and return it as a
 * float Image structure with pixels in the range [0,1].  If the file
 * contains more than one image, then the images will be returned
 * linked by the "next" field of the Image data structure.  
 * See "man pgm" for details on PGM file format.  This handles only
 * the usual 8-bit "raw" PGM format.  Use xv or the PNM tools (such as
 * pnmdepth) to convert from other formats. */

static img_t *img_read_pgm(FILE *fp) {
    int char1, char2, w, h, max, c1, c2, c3, x, y;
    img_t *img;

    char1 = fgetc(fp);
    char2 = fgetc(fp);
    skip_comments(fp);
    c1 = fscanf(fp, "%d", &w);
    skip_comments(fp);
    c2 = fscanf(fp, "%d", &h);
    skip_comments(fp);
    c3 = fscanf(fp, "%d", &max);

    printf("[pgm.c] (%d, %d, %d)\n", w, h, max);

    if (char1 != 'P' || char2 != '5' || 
	c1 != 1 || c2 != 1 || c3 != 1 || max > 255) {
	printf("Input is not a standard raw 8-bit PGM file.\n"
	       "Use xv or pnmdepth to convert file to 8-bit PGM format.\n");
	return NULL;
    }
    
    fgetc(fp);  /* Discard exactly one byte after header. */

    /* Create floating point image with pixels in range [0,1]. */
    img = img_new(w, h);

    for (y = h - 1; y > 0; y--) {
	color_t *c = img->pixels + y * w;

	for (x = 0; x < w; x++) {
	    int v = (int) fgetc(fp);
	    // img_set_pixel(img, x, y, v, v, v);
	    c->r = c->g = c->b = v;
	    img_set_valid_pixel(img, x, y);

	    c++;
	}
    }

    return img;
}
