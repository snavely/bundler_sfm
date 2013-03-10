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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "binvolume.h"
#include "bmp.h"
#include "defines.h"
#include "image.h"
#include "pad.h"

#define BIN_VOLUME_HEADER "BVOL"

// #define NUMBITS(b) (BIT_PAD_WORD((b)->w) * BIT_PAD_WORD((b)->h) * BIT_PAD_WORD((b)->d))
#define NUMBITS(b) (BIT_PAD_WORD((b)->d * BIT_PAD_WORD((b)->h * BIT_PAD_WORD((b)->w))))
#define NUMBYTES(b) ((NUMBITS(b) >> 3) + ((NUMBITS(b) % 8 == 0) ? 0 : 1))

#define BVOL_PLANE_BYTES(b) (BIT_PAD_WORD((b)->h * BIT_PAD_WORD((b)->w)) / 8)


/* Turn on a bit in the volume */
void bvol_setbit(bvol_t *b, int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 || x >= (int) b->w || y >= (int) b->h || z >= (int) b->d) {
        return;
    } else {
        if (b->type == BVOL_BITS) {
            int bit_num = z * BIT_PAD_WORD(b->h * BIT_PAD_WORD(b->w)) + 
                y * BIT_PAD_WORD(b->w) + x;
            int byte_num = (bit_num >> 3);
            int bit_idx = 7 - (bit_num % 8);
            b->data.bits[byte_num] |= (1 << bit_idx);
        } else {
            b->data.heights[y * b->w + x] = z;
        }
    }
}

/* Get the value of a bit in the volume */
u_int8_t bvol_getbit(bvol_t *b, int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 || x >= (int) b->w || y >= (int) b->h || z >= (int) b->d) {
        return 0;
    } else {
        if (b->type == BVOL_BITS) {
            int bit_num = z * BIT_PAD_WORD(b->h * BIT_PAD_WORD(b->w)) + 
                y * BIT_PAD_WORD(b->w) + x;
            int byte_num = (bit_num >> 3);
            int bit_idx = 7 - (bit_num % 8);
            return (b->data.bits[byte_num] & (1 << bit_idx)) ? 1 : 0;
        } else {
            return (b->data.heights[y * b->w + x] == z) ? 1 : 0;
        }
    }
}

/* Returns the number of bits set in the binary volume */
u_int32_t bvol_num_points(bvol_t *b) {
    if (b->type == BVOL_BITS) {
        int x, y, z;
        u_int32_t num_bits = 0;
        
        for (z = 0; z < (int) b->d; z++)
            for (y = 0; y < (int) b->h; y++)
                for (x = 0; x < (int) b->w; x++)
                    if (bvol_getbit(b, x, y, z))
                        num_bits++;

        return num_bits;
    } else {
        return b->w * b->h;
    }
}

/* Create and return an empty binary volume */
bvol_t *new_bin_volume(int w, int h, int d, bvol_type_t type) {
    bvol_t *b = malloc(sizeof(bvol_t));
    b->w = w; b->h = h; b->d = d;
    b->type = type;

    if (type == BVOL_BITS)
        b->data.bits = calloc(NUMBYTES(b), 1);
    else {
        if (d != 256) 
            printf("Error: binary volume with heights should have depth 256\n");
        d = 256;
        b->data.heights = calloc(b->w * b->h, 1);
    }

    return b;
}

/* Create a copy of a binary volume */
bvol_t *bvol_copy(bvol_t *b) {
    bvol_t *bnew = new_bin_volume(b->w, b->h, b->d, b->type);
    
    if (bnew->type == BVOL_BITS)
        memcpy(bnew->data.bits, b->data.bits, NUMBYTES(bnew));
    else
        memcpy(bnew->data.heights, b->data.heights, bnew->w * bnew->h);
    
    return bnew;
}

double bvol_gradient(bvol_t *b, int x, int y, int z) {
    double dx, dy;
    double max_grad = 2 * b->d * b->d;

    if (x == 0)
        dx = BVOL_HEIGHT(b, x + 1, y) - BVOL_HEIGHT(b, x, y);
    else if (x == b->w - 1)
        dx = BVOL_HEIGHT(b, x, y) - BVOL_HEIGHT(b, x - 1, y);
    else 
        dx = BVOL_HEIGHT(b, x + 1, y) - BVOL_HEIGHT(b, x - 1, y);

    if (y == 0)
        dy = BVOL_HEIGHT(b, x, y + 1) - BVOL_HEIGHT(b, x, y);
    else if (y == b->h - 1)
        dy = BVOL_HEIGHT(b, x, y) - BVOL_HEIGHT(b, x, y - 1);
    else
        dy = BVOL_HEIGHT(b, x, y + 1) - BVOL_HEIGHT(b, x, y - 1);
    
    return sqrt(dx * dx + dy * dy) / max_grad;
}

bvol_t *bmp2bvm(bmp_t *bmp) {
    int w = BMP_WIDTH(bmp), h = BMP_HEIGHT(bmp), d = 256;
    int x, y;
    bvol_t *b = new_bin_volume(w, h, d, BVOL_HEIGHTS);

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            bvol_setbit(b, x, y, BMP_GET_PIXEL(bmp, x, y));
        }
    }

    return b;
}

bvol_t *img2bvm(img_t *img) {
    bvol_t *bvol = new_bin_volume(img->w, img->h, 256, BVOL_HEIGHTS);
    int x, y;
    
    for (y = 0; y < img->h; y++) {
        for (x = 0; x < img->w; x++) {
	    int z = img->pixels[y * img->w + x].r;
	    bvol_setbit(bvol, x, y, z);
        }
    }

    return bvol;
}

img_t *bvm2img(bvol_t *b) {
    img_t *img = NULL;
    int x, y, z;
    int set;

    if (b->d != 256) {
        printf("Cannot create image\n");
        return NULL;
    }
    
    img = img_new(b->w, b->h);
    
    for (y = 0; y < (int) b->h; y++) {
        for (x = 0; x < (int) b->w; x++) {
            set = 0;
            for (z = 0; z < (int) b->d; z++) {
                if (bvol_getbit(b, x, y, z)) {
                    img_set_pixel(img, x, y, z, z, z);
                    set = 1;
                }
            }

            if (!set) {
                img_set_pixel(img, x, y, 0, 0, 0);
            }
        }
    }
    
    return img;
}

bmp_t *bvm2bmp(bvol_t *b) {
    bmp_t *bmp = malloc(sizeof(bmp_t));

    bmp_file_header_t file_header;
    bmp_info_header_t info_header;

    int x, y, z = 0;
    int data_bytes;
    int set;

    if (b->d != 256) {
        printf("Cannot create bitmap\n");
        free(bmp);
        return NULL;
    }

    data_bytes = BYTE_PAD_WORD(3 * b->w) * b->h;

    /* Fill the file header */
    file_header.filesize = 14 + sizeof(bmp_info_header_t) + data_bytes;
    file_header.offset = 14 + sizeof(bmp_info_header_t);

    /* Fill the info header */
    info_header.header_size = sizeof(bmp_info_header_t);
    info_header.width = b->w;
    info_header.height = b->h;
    info_header.num_planes = 1;
    info_header.bits_per_pixel = 24;
    info_header.compression_type = 0;
    info_header.image_size = data_bytes;
    info_header.x_pixels_per_meter = 2834;
    info_header.y_pixels_per_meter = 2834;
    info_header.colors_used = 0x0;
    info_header.colors_important = 0x0;
  
    /* Fill the bitmap struct */
    bmp->file_header = file_header;
    bmp->info_header = info_header;
  
    /* Fill in the pixel data */
    bmp->pixels = malloc(sizeof(color_t) * b->w * b->h);
    for (y = 0; y < (int) b->h; y++) {
        for (x = 0; x < (int) b->w; x++) {
            set = 0;
            for (z = 0; z < (int) b->d; z++) {
                if (bvol_getbit(b, x, y, z)) {
                    bmp->pixels[y * b->w + x].r = z;
                    bmp->pixels[y * b->w + x].g = z;
                    bmp->pixels[y * b->w + x].b = z;
                    set = 1;
                }
            }

            if (!set) {
                bmp->pixels[y * b->w + x].r = 0xff;
                bmp->pixels[y * b->w + x].g = 0xff;
                bmp->pixels[y * b->w + x].b = 0xff;
            }
        }
    }

    return bmp;
}

void write_bin_volume(FILE *f, bvol_t *b) {
    u_int32_t type;

    /* Write the header */
    fwrite(BIN_VOLUME_HEADER, 4, 1, f);
    
    /* Write the width, height, depth */
    fwrite(&b->w, 4, 1, f);
    fwrite(&b->h, 4, 1, f);
    fwrite(&b->d, 4, 1, f);
    
    /* Write the type */
    type = b->type;
    fwrite(&type, 4, 1, f);

    /* Write the data */
    if (type == BVOL_BITS) {
        int bytes = NUMBYTES(b);
        fwrite(b->data.bits, 1, bytes, f);
    } else {
        fwrite(b->data.heights, 1, b->w * b->h, f);
    }
}

bvol_t *read_bin_volume(FILE *f) {
    bvol_t *b = NULL;
    char hdr[5] = "";
    u_int32_t w, h, d, type, bytes;

    /* Read the header */
    fread(hdr, 1, 4, f);
    hdr[4] = 0;
    
    if (strcmp(hdr, BIN_VOLUME_HEADER) != 0) {
        printf("Error: file is not a valid binary volume\n");
        fclose(f);
        return NULL;
    }
    
    /* Read the width, height, and depth */
    fread(&w, 4, 1, f);
    fread(&h, 4, 1, f);
    fread(&d, 4, 1, f);

    /* Read the type */
    fread(&type, 4, 1, f);
    b = new_bin_volume(w, h, d, type);
    
    /* Read in the data */
    if (type == BVOL_BITS) {
        bytes = NUMBYTES(b);
        fread(b->data.bits, 1, bytes, f);
    } else {
        fread(b->data.heights, 1, w * h, f);
    }

    return b;
}

/* Slice a plane out of the given binary volume, storing the result in 
 * bimg */
bimg_t *slice_bin_volume(bvol_t *bvol, int z) {
    bimg_t *bimg = new_bin_image(bvol->w, bvol->h);

    if (bvol->type == BVOL_BITS) {
        memcpy(bimg->bits, 
               bvol->data.bits + z * BVOL_PLANE_BYTES(bvol), 
               BVOL_PLANE_BYTES(bvol));
    } else {
        int x, y;
        for (y = 0; y < (int) bvol->h; y++) {
            for (x = 0; x < (int) bvol->w; x++) {
                if (bvol_getbit(bvol, x, y, z))
                    bimg_setbit(bimg, x, y);
            }
        }
    }

    return bimg;
}

void free_bin_volume(bvol_t *b) {
    if (b->type == BVOL_BITS)
        free(b->data.bits);
    else
        free(b->data.heights);
    free(b);
}

