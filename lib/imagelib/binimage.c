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

/* binimage.c */
/* Routines dealing with binary images */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "binimage.h"
#include "bmp.h"
#include "pad.h"
#include "transform.h"
#include "util.h"

#define BIN_IMAGE_HEADER "BIMG"

#define NUMBITS(b) (BIT_PAD_WORD(BIT_PAD_WORD((b)->w) * ((b)->h)))
#define NUMBYTES(b) ((NUMBITS(b) >> 3) + ((NUMBITS(b) % 8 == 0) ? 0 : 1))

/* Turn on a bit in the image */
void bimg_setbit(bimg_t *b, int x, int y) {
    if (x < 0 || y < 0 || x >= (int) b->w || y >= (int) b->h) {
        return;
    } else {
        int bit_num = y * BIT_PAD_WORD(b->w) + x;
        int byte_num = (bit_num >> 3);
        int bit_idx = 7 - (bit_num % 8);
        b->bits[byte_num] |= (1 << bit_idx);
    }
}

/* Get the value of a bit in the image */
u_int8_t bimg_getbit(bimg_t *b, int x, int y) {
    if (x < 0 || y < 0 || x >= (int) b->w || y >= (int) b->h) {
        return 0;
    } else {
        int bit_num = y * BIT_PAD_WORD(b->w) + x;
        int byte_num = (bit_num >> 3);
        int bit_idx = 7 - (bit_num % 8);
        return (b->bits[byte_num] & (1 << bit_idx)) ? 1 : 0;
    }
}

bimg_t *new_bin_image(int w, int h) {
    bimg_t *b = malloc(sizeof(bimg_t));
    int num_bytes;

    b->w = w; b->h = h;
    num_bytes = NUMBYTES(b);
    b->bits = calloc(num_bytes, 1);

    return b;
}

void bimg_clear(bimg_t *b) 
{
    int bytes = NUMBYTES(b);
    int i;
    
    for (i = 0; i < bytes; i++) {
        b->bits[i] = 0;
    }
}

int write_bin_image(FILE *f, bimg_t *b) {
    /* Write the header */
    fwrite(BIN_IMAGE_HEADER, 1, 4, f);
    
    /* Write the width and height */
    fwrite(&b->w, 4, 1, f);
    fwrite(&b->h, 4, 1, f);
    
    /* Write the data */
    fwrite(b->bits, 1, NUMBYTES(b), f);

    return 0;
}

bimg_t *read_bin_image(FILE *f) {
    bimg_t *b = NULL;
    char hdr[5] = "";
    u_int32_t w, h, bytes;

    /* Read the header */
    fread(hdr, 1, 4, f);
    hdr[4] = 0;
    
    if (strcmp(hdr, BIN_IMAGE_HEADER) != 0) {
        printf("Error: file is not a valid binary image\n");
        return NULL;
    }

    /* Read the width and height */
    fread(&w, 4, 1, f);
    fread(&h, 4, 1, f);

    b = new_bin_image(w, h);

    /* Read in the data */
    bytes = NUMBYTES(b);
    fread(b->bits, 1, bytes, f);

    return b;
}

bimg_t *bmp2bmg(bmp_t *bmp) {
    bimg_t *b = NULL;
    unsigned int x, y, w, h;

    w = bmp->info_header.width;
    h = bmp->info_header.height;

    b = new_bin_image(w, h);

    /* Read the pixel data */
    for (y = 0; y < b->h; y++) {
        for (x = 0; x < b->w; x++) {
            color_t pixel = bmp->pixels[y * w + x];

            if (pixel.r == 0x0 && pixel.g == 0x0 && pixel.b == 0x0)
                bimg_setbit(b, x, y);
            else if (pixel.r != 0xff || pixel.g != 0xff || pixel.b != 0xff)
                printf("Weird color at (%d, %d): (%d,%d,%d)\n", x, y, pixel.r, pixel.g, pixel.b);
        }
    }

    return b;
}

bimg_t *img2bmg(img_t *img) {
    bimg_t *b = NULL;
    unsigned int x, y, w, h;

    w = img->w;
    h = img->h;

    b = new_bin_image(w, h);

    /* Read the pixel data */
    for (y = 0; y < b->h; y++) {
        for (x = 0; x < b->w; x++) {
            color_t pixel = img_get_pixel(img, x, y);

            if (pixel.r == 0x0 && pixel.g == 0x0 && pixel.b == 0x0)
                bimg_setbit(b, x, y);
            else if (pixel.r != 0xff || pixel.g != 0xff || pixel.b != 0xff)
                printf("Weird color at (%d, %d): (%d,%d,%d)\n", x, y, pixel.r, pixel.g, pixel.b);
        }
    }

    return b;
}

bmp_t *bmg2bmp(bimg_t *b) {
    bmp_t *bmp = malloc(sizeof(bmp_t));

    bmp_file_header_t file_header;
    bmp_info_header_t info_header;
    palette_t palette;

    color_t black = { 0x0, 0x0, 0x0, 0x0 };
    color_t white = { 0xff, 0xff, 0xff, 0xff };

    unsigned int x, y, count = 0;
    int pad_bytes, data_bytes;

    pad_bytes = BYTE_PAD_WORD(b->w) - b->w;
    data_bytes = NUMBYTES(b) + pad_bytes * b->h;

    /* Fill the file header */
    file_header.filesize = 14 + sizeof(bmp_info_header_t) + sizeof(color_t) * 2 + data_bytes;

    file_header.offset = 14 + sizeof(bmp_info_header_t) + sizeof(color_t) * 2;

    /* Fill the info header */
    info_header.header_size = sizeof(bmp_info_header_t);
    info_header.width = b->w;
    info_header.height = b->h;
    info_header.num_planes = 1;
    info_header.bits_per_pixel = 1;
    info_header.compression_type = 0;
    info_header.image_size = data_bytes;
    info_header.x_pixels_per_meter = 2834;
    info_header.y_pixels_per_meter = 2834;
    info_header.colors_used = 2;
    info_header.colors_important = 2;
  
    /* Fill the palette */
    palette.num_colors = 2;
    palette.colors = malloc(sizeof(color_t) * 2);
    palette.colors[0] = black;
    palette.colors[1] = white;

    /* Fill the bitmap struct */
    bmp->file_header = file_header;
    bmp->info_header = info_header;
    bmp->palette = palette;
  
    /* Fill in the pixel data */
    bmp->pixels = malloc(sizeof(color_t) * b->w * b->h);
    for (y = 0; y < b->h; y++) {
        for (x = 0; x < b->w; x++, count++) {
            bmp->pixels[count] = (bimg_getbit(b, x, y) ? black : white);
        }
    }

    return bmp;
}

img_t *bmg2img(bimg_t *b) {
    img_t *img = img_new(b->w, b->h);

    unsigned int x, y, count = 0;
    int pad_bytes, data_bytes;

    pad_bytes = BYTE_PAD_WORD(b->w) - b->w;
    data_bytes = NUMBYTES(b) + pad_bytes * b->h;

    for (y = 0; y < b->h; y++) {
        for (x = 0; x < b->w; x++, count++) {
            if (bimg_getbit(b, x, y))
                img_set_pixel(img, x, y, 0xff, 0xff, 0xff);
            else
                img_set_pixel(img, x, y, 0x0, 0x0, 0x0);
        }
    }

    return img;
}

void free_bin_image(bimg_t *b) {
    free(b->bits);
    free(b);
}

void print_bin_image(bimg_t *b) {
    int x, y;

    for (y = b->h - 1; y >= 0; y--) {
        for (x = 0; x < b->w; x++) {
            printf(" %c ", (bimg_getbit(b, x, y) ? '*' : ' '));    
        }
        printf("\n");
    }
}

bimg_t *transform_bin_image(bimg_t *b, trans2D_t *T) {
    /* Test each point in the plane for inclusion in the new
     * (transformed) binary image */
    bimg_t *Tb = new_bin_image(b->w, b->h);

    /* Compute the inverse of T */
    trans2D_t *Tinv = transform_invert(T);

    unsigned int x, y;
    double xnew, ynew;
    int xc, yc;

    for (y = 0; y < b->h; y++) {
        for (x = 0; x < b->w; x++) {
            // x -= b->w / 2;
            // y -= b->h / 2;
            transform_point(Tinv, (double) x, (double) y, &xnew, &ynew);
            // x += b->w / 2;
            // y += b->h / 2;
            // xnew += b->w / 2;
            // ynew += b->h / 2;

            /* Check if the point lies in the set */
            xc = iround(xnew);
            yc = iround(ynew);

            if (bimg_getbit(b, xc, yc))
                bimg_setbit(Tb, x, y);
        }
    }

    transform_free(Tinv);

    return Tb;
}

/* Return a list of all points that are turned on in the binary image*/
iv2_t *bin_image_point_list(bimg_t *b, int *num_pts) {
    int n = 0, idx;
    unsigned int x, y;
    iv2_t *pts = NULL;

    /* Count the number of points turned on */
    for (y = 0; y < b->h; y++) {
	for (x = 0; x < b->w; x++) {
	    if (bimg_getbit(b, x, y))
		n++;
	}
    }

    *num_pts = n;

    pts = (iv2_t *)malloc(sizeof(iv2_t) * n);

    idx = 0;
    for (y = 0; y < b->h; y++) {
	for (x = 0; x < b->w; x++) {
	    if (bimg_getbit(b, x, y)) {
		pts[idx++] = iv2_new(x, y);
	    }
	}
    }
    
    return pts;
}
