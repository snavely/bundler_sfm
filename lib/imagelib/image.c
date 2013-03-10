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

/* image.c */

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <float.h>

#include "bmp.h"
#include "fileio.h"
#include "defines.h"
#include "fileio.h"
#include "filter.h"
#include "image.h"
#include "pad.h"
// #include "random.h"
#include "transform.h"
#include "util.h"
#include "vector.h"

#define IMG_FILE_HEADER "IMGX"

#define NULL_PIXEL 0xcc

#define BITS2BYTES(b) (((b) >> 3) + (((b) % 8) == 0 ? 0 : 1))

/* Create a new image with the given width and height */
img_t *img_new(int w, int h) {
    img_t *img = malloc(sizeof(img_t));
    img->w = w; img->h = h;
    img->pixels = calloc(sizeof(color_t), w * h);
    img->origin = v2_new(0.0, 0.0);

    /* All pixels start out as invalid */
    img->pixel_mask = (u_int8_t *)calloc(BITS2BYTES(w * h), sizeof(u_int8_t));

    return img;
}

/* Create a new float image with the given width and height */
fimg_t *fimg_new(int w, int h) {
    fimg_t *img = malloc(sizeof(fimg_t));
    img->w = w; img->h = h;
    img->pixels = calloc(sizeof(float), w * h);
    img->origin = v2_new(0.0, 0.0);

    /* All pixels start out as invalid */
    img->pixel_mask = (u_int8_t *)calloc(BITS2BYTES(w * h), sizeof(u_int8_t));

    return img;
}

color_t img_get_pixel(img_t *img, int x, int y) {
    /* Find a pixel inside the image */
    if (x < 0) 
	x = 0;
    else if (x >= img->w)
	x = img->w - 1;
    
    if (y < 0)
	y = 0;
    else if (y >= img->h)
	y = img->h -1;

    return img->pixels[y * img->w + x];
}

float fimg_get_pixel(fimg_t *img, int x, int y) {
    /* Find a pixel inside the image */
    if (x < 0) 
	x = 0;
    else if (x >= img->w)
	x = img->w - 1;
    
    if (y < 0)
	y = 0;
    else if (y >= img->h)
	y = img->h -1;

    return img->pixels[y * img->w + x];
}

/* Set the alpha on a pixel */
void img_set_alpha(img_t *img, int x, int y, u_int8_t a) {
    color_t *p = img->pixels + y * img->w + x;
    p->extra = a;
}

/* Fill in the given array with the 2x2 square of pixels whose
 * bottom-left corner is (x,y) */
void img_get_pixel_square(img_t *img, int x, int y, color_t **sq) {
    int base = y * img->w + x;
    sq[0] = img->pixels + base;

    if (x < img->w - 1)
	sq[1] = img->pixels + base + 1;
    else 
	sq[1] = NULL;

    if (y < img->h - 1)
	sq[2] = img->pixels + base + img->w;
    else
	sq[2] = NULL;

    if (x < img->w - 1 && y < img->h - 1)
	sq[3] = img->pixels + base + img->w + 1;
    else
	sq[3] = NULL;
}

int img_pixel_in_range(img_t *img, double x, double y)
{
    if (x < 0.0 || x >= img->w - 1 || y < 0.0 || y >= img->h - 1)
	return 0;
    
    return 1;
}

fcolor_t img_average_2x2_square(img_t *img, int x, int y) 
{
    color_t *sq[4];
    fcolor_t result = fcolor_new(0.0, 0.0, 0.0);
    int i;

    img_get_pixel_square(img, x, y, sq);

    for (i = 0; i < 4; i++) {
	result.r += (float) sq[i]->r;
	result.g += (float) sq[i]->g;
	result.b += (float) sq[i]->b;
    }
    
    result.r *= 0.25;
    result.g *= 0.25;
    result.b *= 0.25;

    return result;
}

void img_set_valid_pixel(img_t *img, int x, int y) 
{
	int idx, byte, bit;
    if (x < 0 || y < 0 || x >= img->w || y >= img->h)
	printf("[img_set_valid_pixel] Error: pixel (%d, %d) "
	       "out of range (%d, %d)\n", x, y, img->w, img->h);

    idx = y * img->w + x;
    byte = (idx >> 3);
    bit = idx - (byte << 3);
    
    img->pixel_mask[byte] |= (1 << bit);
}

void fimg_set_valid_pixel(fimg_t *img, int x, int y) 
{
	int idx, byte, bit;

    if (x < 0 || y < 0 || x >= img->w || y >= img->h)
	printf("[img_set_valid_pixel] Error: pixel (%d, %d) "
	       "out of range (%d, %d)\n", x, y, img->w, img->h);

    idx = y * img->w + x;
    byte = (idx >> 3);
    bit = idx - (byte << 3);
    
    img->pixel_mask[byte] |= (1 << bit);
}


void img_set_all_valid_pixels(img_t *img) 
{
    int x, y;
    
    for (y = 0; y < img->h; y++)
	for (x = 0; x < img->w; x++)
	    img_set_valid_pixel(img, x, y);
}

void img_set_pixel(img_t *img, int x, int y, 
                   u_int8_t r, u_int8_t g, u_int8_t b) 
{
    color_t *p = img->pixels + y * img->w + x;

    p->r = r;
    p->g = g;
    p->b = b;
    p->extra = 0;

    /* Mark the pixel as valid */
    img_set_valid_pixel(img, x, y);
}

void fimg_set_pixel(fimg_t *img, int x, int y, float p) 
{
    img->pixels[y * img->w + x] = p;

    /* Mark the pixel as valid */
    fimg_set_valid_pixel(img, x, y);    
}

void img_set_RGBA(img_t *img, int x, int y, 
		  u_int8_t r, u_int8_t g, u_int8_t b, u_int8_t a)
{
    color_t *p = img->pixels + y * img->w + x;

    p->r = r;
    p->g = g;
    p->b = b;
    p->extra = a;

    /* Mark the pixel as valid */
    img_set_valid_pixel(img, x, y);    
}

int img_pixel_is_valid(img_t *img, int x, int y) 
{
    /* Check if pixel is inside the image */
    if (x < 0 || x >= img->w || y < 0 || y >= img->h) {
	return 0;
    } else {
	int idx = y * img->w + x;
	int byte = (idx >> 3);
	int bit = idx - (byte << 3);

	return ((img->pixel_mask[byte] & (1 << bit)) == 0) ? 0 : 1;
    }
}

/* Invalidate the given pixel */
int img_invalidate_pixel(img_t *img, int x, int y) {
    /* Check if pixel is inside the image */
    if (x < 0 || x >= img->w || y < 0 || y >= img->h) {
	return 0;
    } else {
	int idx = y * img->w + x;
	int byte = (idx >> 3);
	int bit = idx - (byte << 3);

	img->pixel_mask[byte] &= ~(1 << bit);
    }
    
    return 0;
}

#if 0
void img_nullify_pixel(img_t *img, int x, int y) {
    color_t *p = img->pixels + y * img->w + x;

    p->r = 0x0;
    p->g = 0x0;
    p->b = 0x0;
    p->extra = NULL_PIXEL;
}
#endif

/* Convert from a bitmap to an image */
img_t *bmp2img(bmp_t *bmp) {
    img_t *img = img_new(BMP_WIDTH(bmp), BMP_HEIGHT(bmp));
    memcpy(img->pixels, bmp->pixels, sizeof(color_t) * img->w * img->h);
    img_set_all_valid_pixels(img);

    return img;
}

/* Convert from an image to grayscale bitmap */
bmp_t *img2bmp_gray(img_t *img, int blue) {
    bmp_t *bmp = malloc(sizeof(bmp_t));
    int w = img->w, h = img->h;
    int data_bytes = h * BYTE_PAD_WORD(w);
    int i;

    /* Setup the file header */
    bmp->file_header.filesize = 14 + sizeof(bmp_info_header_t) + 1024 + data_bytes;
    bmp->file_header.offset = 14 + sizeof(bmp_info_header_t) + 1024;

    // printf("filesize is %d\n", 14 + sizeof(bmp_info_header_t) + data_bytes);

    /* Yay!  Make the palette */
    bmp->palette.num_colors = 256;
    bmp->palette.colors = (color_t *)malloc(sizeof(color_t) * 256);

    for (i = 0; i < 256; i++) {
	bmp->palette.colors[i].r = i;
	bmp->palette.colors[i].g = i;
	bmp->palette.colors[i].b = i;
	bmp->palette.colors[i].extra = 0x0;
    }

    if (blue) 
	bmp->palette.colors[0].b = 0xff;

    /* Setup the info header */
    bmp->info_header.header_size = sizeof(bmp_info_header_t);
    bmp->info_header.width = w;
    bmp->info_header.height = h;
    bmp->info_header.num_planes = 1;
    bmp->info_header.bits_per_pixel = 8;
    bmp->info_header.compression_type = 0;
    bmp->info_header.image_size = data_bytes;
    bmp->info_header.x_pixels_per_meter = 2834;
    bmp->info_header.y_pixels_per_meter = 2834;
    bmp->info_header.colors_used = 256;
    bmp->info_header.colors_important = 0;
    
    /* Fill in the pixel data */
    bmp->pixels = malloc(sizeof(color_t) * w * h);
    memcpy(bmp->pixels, img->pixels, sizeof(color_t) * w * h);

    if (blue) {
	/* Change black pixels to blue */
	for (i = 0; i < w * h; i++) {
	    if (bmp->pixels[i].r == 0x0 && 
		bmp->pixels[i].g == 0x0 && 
		bmp->pixels[i].b == 0x0)
		bmp->pixels[i].b = 0xff;
	}
    }

    return bmp;
}

/* Convert from an image to a bitmap */
bmp_t *img2bmp(img_t *img) {
    bmp_t *bmp = malloc(sizeof(bmp_t));
    int w = img->w, h = img->h;
    int data_bytes = h * BYTE_PAD_WORD(3 * w);

    /* Setup the file header */
    bmp->file_header.filesize = 14 + sizeof(bmp_info_header_t) + data_bytes;
    bmp->file_header.offset = 14 + sizeof(bmp_info_header_t);

    /* Setup the info header */
    bmp->info_header.header_size = sizeof(bmp_info_header_t);
    bmp->info_header.width = w;
    bmp->info_header.height = h;
    bmp->info_header.num_planes = 1;
    bmp->info_header.bits_per_pixel = 24;
    bmp->info_header.compression_type = 0;
    bmp->info_header.image_size = data_bytes;
    bmp->info_header.x_pixels_per_meter = 2834;
    bmp->info_header.y_pixels_per_meter = 2834;
    bmp->info_header.colors_used = 0;
    bmp->info_header.colors_important = 0;
    bmp->palette.num_colors = 0;
    bmp->palette.colors = NULL;
    
    /* Fill in the pixel data */
    bmp->pixels = malloc(sizeof(color_t) * w * h);
    memcpy(bmp->pixels, img->pixels, sizeof(color_t) * w * h);

    return bmp;
}

bmp_t *img2bmp32(img_t *img) {
    bmp_t *bmp = malloc(sizeof(bmp_t));
    int w = img->w, h = img->h;
    int data_bytes = h * BYTE_PAD_WORD(4 * w);

    /* Setup the file header */
    bmp->file_header.filesize = 14 + sizeof(bmp_info_header_t) + data_bytes;
    bmp->file_header.offset = 14 + sizeof(bmp_info_header_t);

    /* Setup the info header */
    bmp->info_header.header_size = sizeof(bmp_info_header_t);
    bmp->info_header.width = w;
    bmp->info_header.height = h;
    bmp->info_header.num_planes = 1;
    bmp->info_header.bits_per_pixel = 32;
    bmp->info_header.compression_type = 0;
    bmp->info_header.image_size = data_bytes;
    bmp->info_header.x_pixels_per_meter = 2834;
    bmp->info_header.y_pixels_per_meter = 2834;
    bmp->info_header.colors_used = 0;
    bmp->info_header.colors_important = 0;
    bmp->palette.num_colors = 0;
    bmp->palette.colors = NULL;
    
    /* Fill in the pixel data */
    bmp->pixels = malloc(sizeof(color_t) * w * h);
    memcpy(bmp->pixels, img->pixels, sizeof(color_t) * w * h);

    return bmp;
}

/* Make a deep copy of the given image */
img_t *img_copy(img_t *img) {
    img_t *copy = img_new(img->w, img->h);
    copy->origin = img->origin;
    memcpy(copy->pixels, img->pixels, sizeof(color_t) * img->w * img->h);
    memcpy(copy->pixel_mask, img->pixel_mask, BITS2BYTES(img->w * img->h));

    return copy;
}
/* Add white noise with the given amplitude to the image.
 * The amplitude may be between 0.0 and 1.0 */
void img_add_white_noise(img_t *img, double a) {
    int x, y;
    double mult = a * 256.0;
    
    for (y = 0; y < img->h; y++) {
        for (x = 0; x < img->w; x++) {
            int n = iround((rand_unit() * 2.0 - 1.0) * mult);
            color_t c = img->pixels[y * img->w + x];

            c.r = CLAMP((int)c.r + n, 0, 255);
            c.g = CLAMP((int)c.g + n, 0, 255);
            c.b = CLAMP((int)c.b + n, 0, 255);
	    c.extra = 0;

            img->pixels[y * img->w + x] = c;
        }
    }
}

/* Convert the image to grayscale (returning a new image) */
img_t *img_convert_grayscale(img_t *img) {
    int x, y;
    img_t *img_gs = img_new(img->w, img->h);
    
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    double intensity = (c.r + c.g + c.b) / 3.0;
	    int gcol = iround(intensity);
	    img_set_pixel(img_gs, x, y, gcol, gcol, gcol);
	}
    }
    
    return img_gs;
}

img_t *img_sub_image(img_t *img, int x, int y, int w, int h) {
    img_t *subimg;
    int xi, yi;

    /* Check that the parameters specify a valid subimage */
    if (x < 0 || y < 0 || x + w > img->w || y + h > img->h)
        return NULL;
    
    subimg = img_new(w, h);
    
    for (yi = y; yi < y + h; yi++) {
        memcpy(subimg->pixels + (yi - y) * w,
               img->pixels + yi * img->w + x,
               sizeof(color_t) * w);
    }

    for (yi = y; yi < y + h; yi++) {
	for (xi = x; xi < x + w; xi++) {
	    if (img_pixel_is_valid(img, xi, yi))
		img_set_valid_pixel(subimg, xi - x, yi - y);
	}
    }

    return subimg;
}


/* Pad an image with empty pixels */
img_t *img_pad_top(img_t *img, int pixels) {
    img_t *out = img_new(img->w, img->h + pixels);
    int x, y;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    if (img_pixel_is_valid(img, x, y)) {
		color_t c = img_get_pixel(img, x, y);
		img_set_pixel(out, x, y, c.r, c.g, c.b);
	    }
	}
    }
    
    return out;
}

img_t *img_pad_bottom(img_t *img, int pixels) {
    img_t *out = img_new(img->w, img->h + pixels);
    int x, y;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    if (img_pixel_is_valid(img, x, y)) {
		color_t c = img_get_pixel(img, x, y);
		img_set_pixel(out, x, y + pixels, c.r, c.g, c.b);
	    }
	}
    }
    
    return out;
}

img_t *img_pad_left(img_t *img, int pixels) {
    img_t *out = img_new(img->w + pixels, img->h);
    int x, y;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    if (img_pixel_is_valid(img, x, y)) {
		color_t c = img_get_pixel(img, x, y);
		img_set_pixel(out, x + pixels, y, c.r, c.g, c.b);
	    }
	}
    }
    
    return out;
}

img_t *img_pad_right(img_t *img, int pixels) {
    img_t *out = img_new(img->w + pixels, img->h);
    int x, y;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    if (img_pixel_is_valid(img, x, y)) {
		color_t c = img_get_pixel(img, x, y);
		img_set_pixel(out, x, y, c.r, c.g, c.b);
	    }
	}
    }
    
    return out;
}

img_t *img_merge(img_t *i1, img_t *i2) {     
    int w = i1->w + i2->w;
    int h = MAX(i1->h, i2->h);
    img_t *img = img_new(w, h);
    int y;

    /* Fill in the pixels */
    for (y = 0; y < h; y++) {
        if (y < i1->h) {
            memcpy(img->pixels + y * w, 
                   i1->pixels + y * i1->w, 
                   sizeof(color_t) * i1->w);
	} else {
            memset(img->pixels + y * w, 1, sizeof(color_t) * i1->w);
	}

        if (y < i2->h) {
            memcpy(img->pixels + y * w + i1->w, 
                   i2->pixels + y * i2->w, 
                   sizeof(color_t) * i2->w);
	} else {
            memset(img->pixels + y * w + i1->w, 1, sizeof(color_t) * i2->w);
	}
    }

    return img;
}

/* Fuse the list of images into a single image.  
 * num = number of arguments */
img_t *img_fuse(int num, int blend, ...) {
    va_list ap;

    img_t **img_list = malloc(sizeof(img_t *) * num);
    int i;
    img_t *img_out;

    int w_new, h_new;
    v2_t min = v2_new(DBL_MAX, DBL_MAX);
    v2_t max = v2_new(-DBL_MAX, -DBL_MAX);
    v2_t origin;

    /* Gather arguments and find the bounds of the new image */
    va_start(ap, blend);
    for (i = 0; i < num; i++) {
	v2_t ll, ur;

	img_t *curr = img_list[i] = va_arg(ap, img_t *);

	ll = curr->origin;
	ur = v2_new(Vx(curr->origin) + curr->w, Vy(curr->origin) + curr->h);

	min = v2_minimum(min, ll);
	max = v2_maximum(max, ur);
    }
    va_end(ap);

    w_new = iround(ceil(Vx(max) - Vx(min) + 1));
    h_new = iround(ceil(Vy(max) - Vy(min) + 1));
    origin = min;

    img_out = img_new(w_new, h_new);
    img_out->origin = origin;
    
    for (i = 0; i < num; i++) {
	img_t *alpha = img_copy(img_list[i]);
	int x, y;

	if (blend == 1) {
	    for (y = 0; y < img_list[i]->h; y++) {
		for (x = 0; x < img_list[i]->w; x++) {
		    color_t pix = img_get_pixel(img_list[i], x, y);
		    pix.r /= num;
		    pix.g /= num;
		    pix.b /= num;

		    img_set_pixel(alpha, x, y, pix.r, pix.g, pix.b);
		}
	    }
	}

	if (blend == 0)
	    img_paste(img_out, alpha);
	else
	    img_blend(img_out, alpha);

	img_free(alpha);
    }

    free(img_list);

    return img_out;
}

/* Paste the source image into the destination image (in place) */
void img_paste(img_t *dest, img_t *src) {
    int x, y;

    for (y = 0; y < src->h; y++) {
	for (x = 0; x < src->w; x++) {
	    int xdest, ydest;
	    color_t col, col_dest;

	    if (!img_pixel_is_valid(src, x, y))
		continue;

	    /* Compute the position of this pixel in the destination */
	    xdest = iround(x + Vx(src->origin) - Vx(dest->origin));
	    ydest = iround(y + Vy(src->origin) - Vy(dest->origin));

	    if (xdest < 0 || xdest > dest->w - 1 || 
		ydest < 0 || ydest > dest->h - 1)
		continue;

	    col = img_get_pixel(src, x, y);
	    col_dest = img_get_pixel(dest, xdest, ydest);

	    if (col_dest.r == 0x0 && col_dest.g == 0x0 && col_dest.b == 0x0) {
		img_set_RGBA(dest, xdest, ydest, 
			     col.r, col.g, col.b, col.extra);
	    } else {
		img_set_RGBA(dest, xdest, ydest, 
			     (col.r + col_dest.r) / 2,
			     (col.g + col_dest.g) / 2,
			     (col.b + col_dest.b) / 2,
			     (col.extra + col_dest.extra) / 2);
	    }
	}
    }
}


/* Add the second image to the first */
void img_blend(img_t *dest, img_t *src) {
    int x, y;

    for (y = 0; y < src->h; y++) {
	for (x = 0; x < src->w; x++) {
	    int xdest, ydest;
	    color_t cols, cold;

	    if (!img_pixel_is_valid(src, x, y))
		continue;

	    /* Compute the position of this pixel in the destination */
	    xdest = iround(x + Vx(src->origin) - Vx(dest->origin));
	    ydest = iround(y + Vy(src->origin) - Vy(dest->origin));

	    if (xdest < 0 || xdest > dest->w - 1 || ydest < 0 || ydest > dest->h - 1)
		continue;

	    cols = img_get_pixel(src, x, y);
	    cold = img_get_pixel(dest, xdest, ydest);

	    img_set_pixel(dest, xdest, ydest, 
			  CLAMP(cols.r + cold.r, 0, 255), 
			  CLAMP(cols.g + cold.g, 0, 255), 
			  CLAMP(cols.b + cold.b, 0, 255));
	}
    }
}

void img_write(FILE *f, img_t *img) {
    /* Write the header */
    fwrite(IMG_FILE_HEADER, 1, 4, f);

    /* Write the width, height */
    write_short(&img->w, f);
    write_short(&img->h, f);
    
    /* Write the pixel data */
    fwrite(img->pixels, sizeof(color_t), img->w * img->h, f);
}

img_t *img_read(FILE *f) {
    char hdr[5];
    img_t *img = malloc(sizeof(img_t));
    
    /* Read the file header */
    fread(hdr, 1, 4, f);
    hdr[4] = 0;
    
    if (strcmp(hdr, IMG_FILE_HEADER) != 0) {
        printf("Invalid image file\n");
        return NULL;
    }
    
    /* Read the width, height */
    read_short(&img->w, f);
    read_short(&img->h, f);
    
    /* Read the pixel data */
    img->pixels = malloc(sizeof(color_t) * img->w * img->h);
    fread(img->pixels, sizeof(color_t), img->w * img->h, f);

    return img;    
}

img_t *img_read_bmp_file(char *fname) {
    img_t *img;
    bmp_t *bmp;
    FILE *f = open_file(fname, "rb");

    bmp = read_bmp(f);
    if (bmp == NULL) {
	printf("[img_read_bmp_file] Error reading bitmap %s.\n", fname);
	fclose(f);
	return NULL;
    }

    img = bmp2img(bmp);
    if (img == NULL) {
	printf("[img_read_bmp_file] Error in image conversion.\n");
    }
    

    free_bmp(bmp);
    
    fclose(f);

    return img;
}

void img_write_bmp_file(img_t *img, char *fname) {
    bmp_t *bmp;
    FILE *f = open_file(fname, "wb");
    
    bmp = img2bmp(img);

    if (bmp == NULL) {
	printf("[img_write_bmp_file] Error in image conversion.\n");
	fclose(f);
	return;
    }

    write_bmp(f, bmp);
    free_bmp(bmp);
    fclose(f);
}

void img_write_bmp_file_gray(img_t *img, char *fname, int blue) {
    bmp_t *bmp;
    FILE *f = open_file(fname, "w");
    
    bmp = img2bmp_gray(img, blue);

    if (bmp == NULL) {
	printf("[img_write_bmp_file] Error in image conversion.\n");
	fclose(f);
	return;
    }

    write_bmp(f, bmp);
    free_bmp(bmp);
    fclose(f);
}

void img_free(img_t *img) {
    free(img->pixels);
    free(img->pixel_mask);
    free(img);
}

void fimg_free(fimg_t *img) {
    free(img->pixels);
    free(img->pixel_mask);
    free(img);
}

/* Equalize an image, returning the result as a new image */
img_t *img_equalize(img_t *img) 
{
    double sum = 0.0;
    double mean, variance;
    int x, y;
    img_t *img_eq;
    double range;

    /* Compute the mean image intensity */
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    double I = color_intensity(c);

	    sum += I;
	}
    }
    
    mean = sum / (img->w * img->h);

    /* Compute the variance */
    sum = 0.0;
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    double I = color_intensity(c);

	    sum += (I - mean) * (I - mean);
	}
    }

    variance = sum / (img->w * img->h);
    range = 2.5 * sqrt(variance);

    img_eq = img_new(img->w, img->h);

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    double I = color_intensity(c);
	    double I_new = 128.0 + 128.0 * (I - mean) / range;
	    int I_i;

	    I_new = CLAMP(I_new, 0.0, 255.0);
	    I_i = iround(I_new);

	    img_set_pixel(img_eq, x, y, I_i, I_i, I_i);
	}
    }

    return img_eq;
}

/* Normalize an image, returning the result as a new image */
img_t *img_normalize(img_t *img) 
{
    double Imin = DBL_MAX, Imax = 0.0;
    int x, y;
    img_t *img_norm;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    double I = color_intensity(c);
	    
	    if (I < Imin)
		Imin = I;
	    if (I > Imax)
		Imax = I;
	}
    }
    
    img_norm = img_new(img->w, img->h);

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    double I = color_intensity(c);
	    double Inorm = 255.0 * (I - Imin) / (Imax - Imin);
	    int I_i = iround(Inorm);
	    
	    img_set_pixel(img_norm, x, y, I_i, I_i, I_i);
	}
    }

    return img_norm;
}

/* Compute the mean color of the neighborhood of radius `rad'
 * surrounding the pixel xc, yc */
void img_mean_color(img_t *img, int xc, int yc, int rad, double *mean) {
    int r = 0, g = 0, b = 0;
    int x, y;
    // int diam = 2 * rad - 1;
    int n = 0;

    int xmin = MAX(0, xc - (rad-1));
    int xmax = MIN(img->w - 1, xc + (rad-1));

    int ymin = MAX(0, yc - (rad-1));
    int ymax = MIN(img->h - 1, yc + (rad-1));
    

    for (y = ymin; y <= ymax; y++) {
	for (x = xmin; x <= xmax; x++) {
	    color_t c = img_get_pixel(img, x, y);

	    r += c.r;
	    g += c.g;
	    b += c.b;

	    n++;
	}
    }

    mean[0] = ((double) r) / ((double) n);
    mean[1] = ((double) g) / ((double) n);
    mean[2] = ((double) b) / ((double) n);
}

double img_variance(img_t *img, int xc, int yc, int rad) {
    int x, y;
    // int diam = 2 * rad - 1;
    double rvar = 0.0, gvar = 0.0, bvar = 0.0;
    double mean[3];
    int n = 0;

    int xmin = MAX(0, xc - (rad-1));
    int xmax = MIN(img->w - 1, xc + (rad-1));

    int ymin = MAX(0, yc - (rad-1));
    int ymax = MIN(img->h - 1, yc + (rad-1));

    img_mean_color(img, xc, yc, rad, mean);

    for (y = ymin; y <= ymax; y++) {
	for (x = xmin; x <= xmax; x++) {
	    color_t c = img_get_pixel(img, x, y);

	    rvar += (c.r - mean[0]) * (c.r - mean[0]);
	    gvar += (c.g - mean[1]) * (c.g - mean[1]);
	    bvar += (c.b - mean[2]) * (c.b - mean[2]);

	    n++;
	}
    }

    rvar /= ((double) n);
    gvar /= ((double) n);
    bvar /= ((double) n);

    return sqrt((rvar + gvar + bvar) / 3.0);
}


/* Compute the mean grayscale color of the neighborhood of radius
 * `rad' surrounding the pixel xc, yc */
double img_mean_grayscale(img_t *img, int xc, int yc, int rad) {
    int x, y;
    int out = 0;
    int diam = 2 * rad - 1;

    for (y = yc - (rad-1); y <= yc + (rad-1); y++) {
	for (x = xc - (rad-1); x <= xc + (rad-1); x++) {
	    color_t c = img_get_pixel(img, x, y);
	    out += c.r;
	}
    }
    
    return ((double) out) / ((double) (diam * diam));
}

/* Compute the variance of a neighborhood of radius `rad' surrounding
 * the pixel xc, yc */
double img_variance_grayscale(img_t *img, int xc, int yc, int rad) {
    int x, y;
    int diam = 2 * rad - 1;
    double var = 0;
    double mean = img_mean_grayscale(img, xc, yc, rad);

    for (y = yc - (rad-1); y <= yc + (rad-1); y++) {
	for (x = xc - (rad-1); x <= xc + (rad-1); x++) {
	    color_t c = img_get_pixel(img, x, y);

	    var += (c.r - mean) * (c.r - mean);
	}
    }

    var /= (diam * diam);

    return var;
}

/* Find the maximum and minimum variances of all neighborhoods of size
 * rad in the image */
void img_find_min_max_variance(img_t *img, int rad, double *min, double *max) {
    int x, y;
    double vmax = 0.0, vmin = DBL_MAX;
    
    for (y = rad; y < img->h - rad; y++) {
	for (x = rad; x < img->w - rad; x++) {
	    double var = img_variance_grayscale(img, x, y, rad);
	    if (var > vmax)
		vmax = var;
	    if (var < vmin)
		vmin = var;
	}
    }
    
    *max = vmax;
    *min = vmin;
}

/* Create gradient maps for the given image */
void img_create_gradient_maps(img_t *img, img_t **x_out, img_t **y_out) {
    img_t *gx = img_new(img->w, img->h);
    img_t *gy = img_new(img->w, img->h);
    
    int x, y;
    int xmax = 0, ymax = 0;
    double xratio, yratio;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    if (x == 0 || x == img->w - 1) {
		/* Boundary conditions */
		img_set_pixel(gx, x, y, 0, 0, 0);
	    } else {
		/* Compute the gradient */
		color_t p1 = img_get_pixel(img, x - 1, y);
		color_t p2 = img_get_pixel(img, x + 1, y);
		double i1 = (p1.r + p1.g + p1.b) / 3.0;
		double i2 = (p2.r + p2.g + p2.b) / 3.0;
		int diff = iround(fabs(i2 - i1));

		if (diff > xmax)
		    xmax = diff;

		img_set_pixel(gx, x, y, diff, diff, diff);
	    }

	    if (y == 0 || y == img->h - 1) {
		img_set_pixel(gy, x, y, 0, 0, 0);
	    } else {
		color_t p1 = img_get_pixel(img, x, y - 1);
		color_t p2 = img_get_pixel(img, x, y + 1);
		double i1 = (p1.r + p1.g + p1.b) / 3.0;
		double i2 = (p2.r + p2.g + p2.b) / 3.0;
		int diff = iround(fabs(i2 - i1));

		if (diff > ymax)
		    ymax = diff;

		img_set_pixel(gy, x, y, diff, diff, diff);
	    }
	}
    }

    xratio = 256.0 / xmax;
    yratio = 256.0 / ymax;

    /* Normalize */
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {	 
	    color_t cx = img_get_pixel(gx, x, y);
	    color_t cy = img_get_pixel(gy, x, y);
   
	    int xc = iround(xratio * cx.r);
	    int yc = iround(yratio * cy.r);

	    img_set_pixel(gx, x, y, xc, xc, xc);
	    img_set_pixel(gy, x, y, yc, yc, yc);
	}
    }
    
    *x_out = gx;
    *y_out = gy;
}

/* Compute the magnitude of the gradient of the image at the point (x,y) */
double img_gradient(img_t *img, int x, int y) {
    color_t p1, p2;
    double i1, i2, dx, dy;

    /* Compute partial derivative in x */
    p1 = img_get_pixel(img, x - 1, y);
    p2 = img_get_pixel(img, x + 1, y);
    
    i1 = (double) (p1.r + p1.g + p1.b) / 3.0;
    i2 = (double) (p2.r + p2.g + p2.b) / 3.0;

    dx = 0.5 * (i2 - i1);
    
    /* Compute partial derivative in y */
    p1 = img_get_pixel(img, x, y - 1);
    p2 = img_get_pixel(img, x, y + 1);
    
    i1 = (double) (p1.r + p1.g + p1.b) / 3.0;
    i2 = (double) (p2.r + p2.g + p2.b) / 3.0;

    dy = 0.5 * (i2 - i1);

    return sqrt(dx * dx + dy * dy);
}

void img_draw_pt(img_t *img, int x, int y, int size, u_int8_t r, u_int8_t g, u_int8_t b)
{
    int xp, yp;

    int ymin = iround(MAX(0, y - size / 2 - Vy(img->origin)));
    int ymax = iround(MIN(y + size / 2 - Vy(img->origin), img->h - 1));
    int xmin = iround(MAX(0, x - size / 2 - Vx(img->origin)));
    int xmax = iround(MIN(x + size / 2 - Vx(img->origin), img->w - 1));

    for (yp = ymin; yp <= ymax; yp++) {
	for (xp = xmin; xp <= xmax; xp++) {
	    img_set_pixel(img, xp, yp, r, g, b);
	}
    }
}

/* Draw a line from (x0, y0) to (x1, y1) */
void img_draw_line(img_t *img, int x0, int y0, int x1, int y1, u_int8_t r, u_int8_t g, u_int8_t b)
{
    /* Compute the slope */
    int dx = x1 - x0;
    int dy = y1 - y0;
    int tmp;
    double curr, slope;
    int x, y;

    if (abs(dy) > abs(dx)) {
	if (y0 > y1) {
	    tmp = y0;
	    y0 = y1;
	    y1 = tmp;

	    tmp = x0;
	    x0 = x1;
	    x1 = tmp;
	}

	slope = ((double) dx) / ((double) dy);

	for (y = y0, curr = x0; y <= y1; y++, curr += slope) {
	    x = iround(curr);

	    if (x >= 0 && x < img->w && y >= 0 && y < img->h)
		img_draw_pt(img, x, y, 1, r, g, b);
	}
    } else {
	if (x0 > x1) {
	    tmp = x0;
	    x0 = x1;
	    x1 = tmp;

	    tmp = y0;
	    y0 = y1;
	    y1 = tmp;
	}

	slope = ((double) dy) / ((double) dx);

	for (x = x0, curr = y0; x <= x1; x++, curr += slope) {
	    y = iround(curr);

	    if (x >= 0 && x < img->w && y >= 0 && y < img->h)
		img_draw_pt(img, x, y, 1, r, g, b);
	}	
    }
}

/* Combine a new color with the existing color of pixel (x,y) in the
 * image, weighted to give equal weight to all parts already
 * contributing to the pixel color */
void pixel_combine(img_t *img, int x, int y, double r, double g, double b, int num_parts) {    
    if (num_parts == 1) {
	int ri = iround(r);
	int gi = iround(g);
	int bi = iround(b);
	img_set_pixel(img, x, y, ri, gi, bi);
    } else {
	double ratio0 = ((double) num_parts - 1) / ((double) num_parts);
	double ratio1 = 1.0 / ((double) num_parts);
	color_t old_pixel = img_get_pixel(img, x, y);

	int ri = iround(ratio0 * ((double) old_pixel.r) + ratio1 * r);
	int gi = iround(ratio0 * ((double) old_pixel.g) + ratio1 * g);
	int bi = iround(ratio0 * ((double) old_pixel.b) + ratio1 * b);
	
	img_set_pixel(img, x, y, ri, gi, bi);
    }
}

/* Check if the image of the neighborhood under the given transform
 * belongs to the original image.  This is true if all four corners of
 * the neighborhood are inside the original image */
int nhood_image_within_image(int w, int h, v2_t origin,
			     trans2D_t *T, 
			     double x, double y, int nhood_radius) {
    double Tx, Ty;

    double corners[4][2] = 
	{ { x - nhood_radius, y - nhood_radius },
	  { x + nhood_radius, y - nhood_radius },
	  { x - nhood_radius, y + nhood_radius },
	  { x + nhood_radius, y + nhood_radius } };

    int cidx;

    for (cidx = 0; cidx < 4; cidx++) {
	transform_point(T,
			corners[cidx][0], corners[cidx][1],
			&Tx, &Ty);

	if (Tx < Vx(origin) || Ty < Vy(origin) || Tx > Vx(origin) + w - 1 || Ty > Vy(origin) + h - 1)
	    return 0;
    }
    
    return 1;
}

/* Upsize the image by putting a border of pixels on the right and
 * upper sides */
img_t *img_upsize(img_t *img, int w_new, int h_new) 
{
    img_t *img_out;

    if (w_new == img->w && h_new == img->h) {
	return img_copy(img);
    } else if (w_new < img->w || h_new < img->h) {
	printf("[img_upsize] Error: New dimensions must be no smaller than the old dimensions\n");
	return NULL;
    }

    img_out = img_new(w_new, h_new);
    img_out->origin = img->origin;

    /* Paste the old image into the new image */
    img_paste(img_out, img);

    return img_out;
}

/* Upsize the image so that is a square and width and height are a
 * power of two */
img_t *img_upsize_square_power_of_two(img_t *img)
{
    int w_new, h_new, size;

    w_new = least_larger_power_of_two(img->w);
    h_new = least_larger_power_of_two(img->h);

    size = MAX(w_new, h_new);

    return img_upsize(img, size, size);
}

/* Upsize the image so that width and height are a power of two */
img_t *img_upsize_power_of_two(img_t *img)
{
    int w_new, h_new;

    w_new = least_larger_power_of_two(img->w);
    h_new = least_larger_power_of_two(img->h);

    return img_upsize(img, w_new, h_new);
}

/* Returns true is the given subregion of an image consists of valid
 * pixels */
int img_region_is_valid(img_t *img, int xmin, int xmax, int ymin, int ymax) 
{
    int xp, yp;

    if (xmin < 0 || xmax >= img->w || ymin < 0 || ymax >= img->h)
	return 0;
    
    for (yp = ymin; yp <= ymax; yp++) {
	for (xp = xmin; xp <= xmax; xp++) {
	    if (!img_pixel_is_valid(img, xp, yp))
		return 0;
	}
    }

    return 1;
}

/* Compute the width, height, and origin of the image under the given
 * transform */
void img_compute_Tstats(img_t *img, trans2D_t *T, int *w_out, int *h_out, v2_t *origin_out) {
    int w = img->w, h = img->h, w_new, h_new;
    v2_t crs[4];
    v2_t min = v2_new(DBL_MAX, DBL_MAX);
    v2_t max = v2_new(-DBL_MAX, -DBL_MAX);
    v2_t origin;
    int i;

    crs[0] = v2_new(0, 0);
    crs[1] = v2_new(0, h - 1);
    crs[2] = v2_new(w - 1, 0);
    crs[3] = v2_new(w - 1, h - 1);

    for (i = 0; i < 4; i++) {
	crs[i] = v2_add(crs[i], img->origin);
	crs[i] = transform_vector(T, crs[i]);
	min = v2_minimum(min, crs[i]);
	max = v2_maximum(max, crs[i]);
    }

    w_new = iround(ceil(Vx(max) - Vx(min) + 1));
    h_new = iround(ceil(Vy(max) - Vy(min) + 1));
    origin = min;

    *w_out = w_new;
    *h_out = h_new;
    *origin_out = origin;
}

/* Light up invalid pixels */
void img_light_invalid_pixels(img_t *img) 
{
    int x, y;
    
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    if (!img_pixel_is_valid(img, x, y)) {
		img_set_pixel(img, x, y, 0xff, 0x0, 0xff);
	    }
	}
    }
}

/* Turn off all pixels of a given color */
void img_set_invalid_pixels(img_t *img, unsigned char r, unsigned char g, unsigned char b) {
    int x, y;
    
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);

	    if (c.r == r && c.g == g && c.b == b) {
		img_invalidate_pixel(img, x, y);
	    }
	}
    }
}

/* Convert an image to a 16-bit raw format */
void img_write_16bit_raw(img_t *img, char *filename) {    
    int w = img->w, h = img->h;
    int x, y;
    
    FILE *f = open_file(filename, "w");

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    long int c24 = 
		(((long) c.b) << 16) | (((long) c.g) << 8) | (((long) c.r));
	    unsigned short int c16 = (unsigned short int) (c24 >> 8);
	    
	    fwrite(&c16, sizeof(unsigned short), 1, f);
	}
    }

    fclose(f);
}

/* Shrink-wrap an image */
img_t *img_shrink_wrap(img_t *img) {
    img_t *img_out;
    int w, h, w_new = 0, h_new = 0;
    int x, y;

    w = img->w;
    h = img->h;
    
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    if (img_pixel_is_valid(img, x, y)) {
		w_new = MAX(w_new, x+1);
		h_new = MAX(h_new, y+1);
	    }
	}
    }
    
    img_out = img_new(w_new, h_new);

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    if (img_pixel_is_valid(img, x, y)) {
		color_t c = img_get_pixel(img, x, y);
		img_set_pixel(img_out, x, y, c.r, c.g, c.b);
	    }
	}
    }

    img_out->origin = img->origin;

    return img_out;
}

/* Copy the image data into a buffer */
void img_buffer(img_t *img, unsigned char *buffer) 
{
    // int x, y;
    color_t *c;
    int npixels = img->w * img->h;
    color_t *cend = img->pixels + npixels;
    unsigned char *b = buffer;

    for (c = img->pixels; c != cend; c++, b += 3) {
        b[0] = c->r;
        b[1] = c->g;
        b[2] = c->b;
    }

#if 0    
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    int base = 3 * (y * img->w + x);
	    buffer[base + 0] = c.r;
	    buffer[base + 1] = c.g;
	    buffer[base + 2] = c.b;
	}
    }
#endif
}

/* Copy the image data into a buffer, with alpha */
void img_buffer_alpha(img_t *img, unsigned char *buffer, unsigned char alpha) 
{
    int x, y;
    
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    int base = 4 * (y * img->w + x);
	    buffer[base + 0] = c.r;
	    buffer[base + 1] = c.g;
	    buffer[base + 2] = c.b;

	    if (img_pixel_is_valid(img, x, y))
		buffer[base + 3] = alpha;
	    else
		buffer[base + 3] = 0;
	}
    }
}

/* Copy the (vertically flipped) image data into a buffer */
void img_buffer_flip(img_t *img, unsigned char *buffer)
{
    int x, y;
    
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    int base = 3 * ((img->h - y - 1) * img->w + x);
	    buffer[base + 0] = c.r;
	    buffer[base + 1] = c.g;
	    buffer[base + 2] = c.b;
	}
    }
}

/* Convert the given image coordinates to NDC (normalized device
 * coordinates */
void img_to_NDC(img_t *img, double x, double y, double *x_out, double *y_out)
{
    double size = MAX(img->w, img->h);
    double size2 = 0.5 * size;

    *x_out = (x - 0.5 * img->w) / size2;
    *y_out = (y - 0.5 * img->h) / size2;
}

/* Convert the given image coordinates from NDC (normalized device
 * coordinates */
void img_from_NDC(img_t *img, double x, double y, 
		  double *x_out, double *y_out)
{
    double size = MAX(img->w, img->h);
    double size2 = 0.5 * size;

    *x_out = x * size2 + 0.5 * img->w;
    *y_out = y * size2 + 0.5 * img->h;
}


/* Convert an image to a locally normalized form */
img_t *img_normalize_local(img_t *img, double sigma) 
{
    /* Compute the filter */
    int size;
    double *filter = compute_gaussian_filter(sigma, 2.0, &size);
    int rad = size / 2;
    int x, y;

    img_t *img_out = img_new(img->w, img->h);

    /* Normalize the image */
    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    /* Compute the mean and variance in a region */
	    double total_weight = 0.0;
	    double sum = 0.0;
	    double mean, variance;
	    int dx, dy;
	    double d;
	    int dint;

	    for (dy = -rad; dy <= rad; dy++) {
		for (dx = -rad; dx <= rad; dx++) {
		    int xi = x + dx;
		    int yi = y + dy;
		    double weight;

		    /* Check bounds */
		    if (xi < 0 || xi >= img->w || yi < 0 || yi >= img->h)
			continue;

		    weight = filter[dx + rad] * filter[dy + rad];
		    sum += weight * 
			color_intensity(img_get_pixel(img, xi, yi));

		    total_weight += weight;
		}
	    }

	    mean = sum / total_weight;
	    
	    /* Compute variance */
	    total_weight = 0.0;
	    sum = 0.0;

	    for (dy = -rad; dy <= rad; dy++) {
		for (dx = -rad; dx <= rad; dx++) {
		    int xi = x + dx;
		    int yi = y + dy;

		    double weight, diff;

		    /* Check bounds */
		    if (xi < 0 || xi >= img->w || yi < 0 || yi >= img->h)
			continue;
		
		    weight = filter[dx + rad] * filter[dy + rad];
		    diff = mean - color_intensity(img_get_pixel(img, xi, yi));
		    sum += weight * diff * diff;
		    
		    total_weight += weight;
		}
	    }

	    variance = sum / total_weight;
	    if (variance < 1.0e-4)
		variance = 1.0e-4;

	    d = (color_intensity(img_get_pixel(img, x, y)) - mean) / 
		sqrt(variance);

	    d = d * 64.0 + 127.5;
	    d = CLAMP(d, 0.0, 255.0);
	    dint = iround(d);

	    img_set_pixel(img_out, x, y, dint, dint, dint); 
	    // img.Pixel(x, y, 0) - mean) / sqrt(variance);
	}
    }

    free(filter);

    return img_out;
}

img_t *fimg2img(fimg_t *img) {
    img_t *img_out = img_new(img->w, img->h);

    int x, y;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    float p = fimg_get_pixel(img, x, y);
	    int v = iround(CLAMP(p, 0.0, 255.0));
	    img_set_pixel(img_out, x, y, v, v, v);
	}
    }

    return img_out;
}
