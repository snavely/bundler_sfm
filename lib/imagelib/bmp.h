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

#ifndef __bmp_h__
#define __bmp_h__

#include <stdio.h>
#include <sys/types.h>

#include "color.h"

/*
 * File format for .bmp files 
 * BM
 * 4 - filesize
 * 4 - reserved=0
 * 4 - offset_of_SOD
 * ----
 * 4 - headersize
 * 4 - width
 * 4 - height
 * 2 - numplanes
 * 2 - bitsperpixel
 * 4 - compressiontype
 * 4 - imagesize
 * 4 - Xpixelspermeter
 * 4 - Ypixelspermeter
 * 4 - colorsused
 * 4 - colorsimportant
 * ----
 * {
 * 	1 - red
 * 	1 - green
 * 	1 - blue
 * 	1 - extra
 * }x(1 >> bitsperpixel)
 * 
 * data
 */

typedef struct {
    u_int32_t filesize;
    u_int32_t offset;
} bmp_file_header_t;

typedef struct {
    u_int32_t header_size;
    u_int32_t width;
    u_int32_t height;
    u_int16_t num_planes;
    u_int16_t bits_per_pixel;
    u_int32_t compression_type;
    u_int32_t image_size;
    u_int32_t x_pixels_per_meter;
    u_int32_t y_pixels_per_meter;
    u_int32_t colors_used;
    u_int32_t colors_important;
} bmp_info_header_t;

typedef struct {
    int num_colors;
    color_t *colors;
} palette_t;

typedef struct {
    bmp_file_header_t file_header;
    bmp_info_header_t info_header;
    palette_t palette;
    color_t *pixels;
} bmp_t;

/* Read a bitmap and return it */
bmp_t *read_bmp(FILE *f);
bmp_t *read_bmp_file(char *f);

/* Get the dimensions (w, h) of the given bmp file */
int bmp_file_get_dimensions(char *filename, int *w, int *h);

/* Write a bitmap to the given file */
int write_bmp(FILE *f, bmp_t *bmp);
int write_bmp_file(char *f, bmp_t *bmp);

int bmp_get_width(bmp_t *bmp);
int bmp_get_height(bmp_t *bmp);

/* Create a sub-bitmap from the given bitmap */
bmp_t *sub_bitmap(bmp_t *bmp, int x, int y, int w, int h);

/* Paste two bitmaps together horizontally */
bmp_t *merge_bmps(bmp_t *b1, bmp_t *b2);

/* Free a bitmap structure */
void free_bmp(bmp_t *b);

#define BMP_WIDTH(b) ((b)->info_header.width)
#define BMP_HEIGHT(b) ((b)->info_header.height)
#define BMP_GET_PIXEL(b, x, y) ((b)->pixels[(y) * BMP_WIDTH(b) + (x)].r)

#endif /* __bmp_h__ */
