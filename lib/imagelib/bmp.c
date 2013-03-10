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

/* bmp.c */
/* Routines for reading/writing bmp files */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bmp.h"
#include "defines.h"
#include "fileio.h"
#include "pad.h"

int read_bmp_file_header(FILE *f, bmp_file_header_t *hdr) {
    u_int8_t sig[3];
    u_int32_t word;

    /* Read the signature */
    fread(sig, 1, 2, f);
    sig[2] = 0;
    if (sig[0] != 'B' || sig[1] != 'M')
        return 1;

    /* Read the filesize */
    read_word(&hdr->filesize, f);

    // printf("filesize = %d\n", hdr->filesize);

    /* Read the reserved word */
    read_word(&word, f);

    /* Read the offset of the start of data */
    read_word(&hdr->offset, f);

    return 0;
}

int read_bmp_info_header(FILE *f, bmp_info_header_t *hdr) {
    read_word(&hdr->header_size, f);
    read_word(&hdr->width, f);
    read_word(&hdr->height, f);
    read_short(&hdr->num_planes, f);
    read_short(&hdr->bits_per_pixel, f);
    read_word(&hdr->compression_type, f);
    read_word(&hdr->image_size, f);
    read_word(&hdr->x_pixels_per_meter, f);
    read_word(&hdr->y_pixels_per_meter, f);
    read_word(&hdr->colors_used, f);
    read_word(&hdr->colors_important, f);

    return 0;
}

/* Get the dimensions (w, h) of the given bmp file */
int bmp_file_get_dimensions(char *filename, int *w, int *h) 
{
    FILE *f = fopen(filename, "rb");
    bmp_t b;
    
    if (f == NULL)
	return -1;
    
    /* Read the file header */
    if (read_bmp_file_header(f, &(b.file_header)) != 0) {
	fclose(f);
        return -1;
    }

    /* Read the info header */
    if (read_bmp_info_header(f, &(b.info_header)) != 0) {
	fclose(f);
        return -1;
    }    

    *w = b.info_header.width;
    *h = b.info_header.height;

    fclose(f);

    return 0;
}


int read_palette(FILE *f, u_int16_t bpp, palette_t *p) {
    int i, num_colors;

    p->num_colors = num_colors = (1 << bpp);
    p->colors = malloc(sizeof(color_t) * num_colors);

    for (i = 0; i < num_colors; i++) {
        read_byte(&p->colors[i].b, f);
        read_byte(&p->colors[i].g, f);
        read_byte(&p->colors[i].r, f);
        read_byte(&p->colors[i].extra, f);

#if 0
	if (i == 0) {
	    printf("color[0] = (%d,%d,%d)\n", p->colors[i].r, p->colors[i].g, p->colors[i].b);
	}
#endif
    }
    
    return 0;
}

bmp_t *read_bmp_file(char *f) {
    FILE *fd = open_file(f, "rb");
    bmp_t *out = read_bmp(fd);
    fclose(fd);
    
    return out;
    
}

bmp_t *read_bmp(FILE *f) {
    bmp_t *b = NULL;
    int i, num_pixels, w, h, x, y;

    b = malloc(sizeof(bmp_t));

    /* Read the file header */
    if (read_bmp_file_header(f, &b->file_header) != 0) {
        free(b);
        return NULL;
    }

    /* Read the info header */
    if (read_bmp_info_header(f, &b->info_header) != 0) {
        free(b);
        return NULL;
    }

    num_pixels = b->info_header.width * b->info_header.height;
    b->pixels = malloc(sizeof(color_t) * num_pixels);
    w = b->info_header.width;
    h = b->info_header.height;

    // printf("bpp: %d\n", b->info_header.bits_per_pixel);

    if (b->info_header.bits_per_pixel <= 8) {
        int bpp = b->info_header.bits_per_pixel;
        u_int8_t mask[9] = { 0x0, 0x1, 0x3, 0x0, 0xf, 0x0, 0x0, 0x0, 0xff };

	// printf("image is indexed...\n");

        /* Read the palette */
        if (read_palette(f, b->info_header.bits_per_pixel, &b->palette) != 0) {
            free(b);
            return NULL;
        }

        for (y = 0; y < h; y++) {
            u_int8_t extra[4];
            u_int32_t bytes_read = 0;

            for (x = 0; x < w; ) {
                /* Read a byte */
                u_int8_t byte, col;
                read_byte(&byte, f);
                bytes_read++;
                
                for (i=8 / bpp - 1; i >= 0; i--) {
                    /* Read the ith set of bpp bits */
                    col = (byte >> (i * bpp)) & mask[bpp];
                    b->pixels[y * w + x] = 
                        b->palette.colors[col];
                    x++;

                    if (x >= w)
                        break;
                }
            }

            /* Read the leftover bytes */
            if ((bytes_read % 4) != 0)
                fread(extra, 1, 4 - (bytes_read % 4), f);
        }
    } else {
	if (b->info_header.bits_per_pixel == 24) {
	    /* Skip to the data section */
	    int curr_byte = 6 + sizeof(bmp_file_header_t) + sizeof(bmp_info_header_t);
	    int start_byte = b->file_header.offset;
	    u_int8_t dummy;
		int line_bytes;
		char *line;

	    // printf("curr: %d\nstart: %d\n", curr_byte, start_byte);

	    while (curr_byte < start_byte) {
		read_byte(&dummy, f);
		curr_byte++;
	    }

#if 0
	    for (y = 0; y < h; y++) {
		u_int8_t extra[4];
		u_int32_t bytes_read = 0;

		for (x = 0; x < w; x ++) {
		    read_byte(&b->pixels[y * w + x].b, f);
		    read_byte(&b->pixels[y * w + x].g, f);
		    read_byte(&b->pixels[y * w + x].r, f);
		    b->pixels[y * w + x].extra = 0xff;
		}

		bytes_read = 3 * w;
            
		if ((bytes_read % 4) != 0)
		    fread(extra, 1, 4 - (bytes_read % 4), f);    
	    }
#else
	    line_bytes = 3 * w;

	    if ((line_bytes % 4) != 0)
		line_bytes += (4 - (line_bytes % 4));

	    line = (char *) malloc(sizeof(char) * line_bytes);

	    for (y = 0; y < h; y++) {
		int idx = 0;

		fread(line, 1, line_bytes, f);

		for (x = 0; x < w; x ++) {
		    b->pixels[y * w + x].b = line[idx++];
		    b->pixels[y * w + x].g = line[idx++];
		    b->pixels[y * w + x].r = line[idx++];
		    b->pixels[y * w + x].extra = 0xff;
		}
	    }

	    free(line);
#endif
	} else if (b->info_header.bits_per_pixel == 32) {
	    int curr_offset = 54;
	    int start_of_data = b->file_header.offset;
		char *ignore;

	    if (curr_offset != start_of_data) {
		printf("32-bit bmp, skipping %d bytes\n", 
		       (start_of_data - curr_offset));

		ignore = malloc(start_of_data - curr_offset);

		/* Skip 68 bytes */
		fread(ignore, 1, start_of_data - curr_offset, f);
		free(ignore);
	    }

	    for (y = 0; y < h; y++) {
		u_int8_t extra[4];
		u_int32_t bytes_read = 0;

		for (x = 0; x < w; x ++) {
		    read_byte(&b->pixels[y * w + x].b, f);
		    read_byte(&b->pixels[y * w + x].g, f);
		    read_byte(&b->pixels[y * w + x].r, f);
		    read_byte(&b->pixels[y * w + x].extra, f);
		}

		bytes_read = 4 * w;
            
		if ((bytes_read % 4) != 0)
		    fread(extra, 1, 4 - (bytes_read % 4), f);    
	    }
	}
    }

    return b;
}

int write_bmp_file_header(FILE *f, bmp_file_header_t *file_header) {
    u_int32_t word = 0;

    /* Write the signature */
    fwrite("BM", 1, 2, f);

    /* Write the file size */
    write_word(&file_header->filesize, f);
  
    /* Write the reserved word */
    write_word(&word, f);
  
    /* Write the offset */
    write_word(&file_header->offset, f);
  
    return 0;
}

int write_bmp_info_header(FILE *f, bmp_info_header_t *hdr) {
    write_word(&hdr->header_size, f);
    write_word(&hdr->width, f);
    write_word(&hdr->height, f);
    write_short(&hdr->num_planes, f);
    write_short(&hdr->bits_per_pixel, f);
    write_word(&hdr->compression_type, f);
    write_word(&hdr->image_size, f);
    write_word(&hdr->x_pixels_per_meter, f);
    write_word(&hdr->y_pixels_per_meter, f);
    write_word(&hdr->colors_used, f);
    write_word(&hdr->colors_important, f);

    return 0;
}

int write_palette(FILE *f, u_int16_t bpp, palette_t *p) {
    int i, num_colors;

    num_colors = (1 << bpp);

    for (i = 0; i < num_colors; i++) {
        write_byte(&p->colors[i].b, f);
        write_byte(&p->colors[i].g, f);
        write_byte(&p->colors[i].r, f);
        write_byte(&p->colors[i].extra, f);
    }
    
    return 0;
}

static int colors_are_equal(color_t *c1, color_t *c2) {
    return (c1->r == c2->r && c1->g == c2->g && c1->b == c2->b);
}

int write_bmp_file(char *f, bmp_t *bmp) {
    FILE *fd = open_file(f, "wb");
    int ret = write_bmp(fd, bmp);
    fclose(fd);

    return ret;
}

int write_bmp(FILE *f, bmp_t *bmp) {
    int bpp = bmp->info_header.bits_per_pixel;
    u_int32_t num_colors = (1 << bpp);
    int w = bmp->info_header.width;
    int h = bmp->info_header.height;
    int x, y;

    if (write_bmp_file_header(f, &bmp->file_header) != 0)
        return 1;
    
    if (write_bmp_info_header(f, &bmp->info_header) != 0)
        return 1;
  
    if (bmp->info_header.bits_per_pixel <= 8) {
        if (write_palette(f, bmp->info_header.bits_per_pixel, &bmp->palette) != 0)
            return 1;

        for (y = 0; y < h; y++) {
            int bytes_written = 0;

            for (x = 0; x < w; ) {
                u_int8_t byte = 0;
                int b;

                for (b = 8 / bpp - 1; b >= 0; b--) {
                    color_t *pixel = &bmp->pixels[y * w + x];
                    int col_idx = -1;
                    unsigned int i;

                    /* Match the pixel color to a palette color */
                    for (i = 0; i < num_colors; i++) {
                        if (colors_are_equal(pixel, bmp->palette.colors + i)) {
                            col_idx = i;
                            break;
                        }
                    }

                    if (col_idx == -1) {
                        printf("Error: pixel color was not in palette\n");
                        return 1;
                    }

                    byte |= (col_idx << (b * bpp));
                    x++;

                    if (x >= w)
                        break;
                }

                bytes_written++;
		
                /* Write the byte */
                write_byte(&byte, f);
            }

            /* Write extra padding bytes if needed */
            if ((bytes_written % 4) != 0) {
                u_int8_t extra[4] = { 0, 0, 0, 0 };
                fwrite(extra, 1, 4 - (bytes_written % 4), f);
            }
        }
    } else if (bmp->info_header.bits_per_pixel == 24) {
        for (y = 0; y < h; y++) {
            int bytes_written;

            for (x = 0; x < w; x++) {
                write_byte(&bmp->pixels[y * w + x].b, f);
                write_byte(&bmp->pixels[y * w + x].g, f);
                write_byte(&bmp->pixels[y * w + x].r, f);
            }

            bytes_written = 3 * w;
            
            /* Write extra padding bytes if needed */
            if ((bytes_written % 4) != 0) {
                u_int8_t extra[4] = { 0, 0, 0, 0 };
                fwrite(extra, 1, 4 - (bytes_written % 4), f);
            }
        }
    } else if (bmp->info_header.bits_per_pixel == 32) {
	printf("[write_bmp] Writing 32-bit bmp\n");

        for (y = 0; y < h; y++) {
            for (x = 0; x < w; x++) {
                write_byte(&bmp->pixels[y * w + x].b, f);
                write_byte(&bmp->pixels[y * w + x].g, f);
                write_byte(&bmp->pixels[y * w + x].r, f);
		write_byte(&bmp->pixels[y * w + x].extra, f);
            }
        }	
    } else {
	printf("[write_bmp] Error: invalid bpp\n");
    }
    

    return 0;
}

int bmp_get_width(bmp_t *bmp) {
    return BMP_WIDTH(bmp);
}

int bmp_get_height(bmp_t *bmp) {
    return BMP_HEIGHT(bmp);
}

bmp_t *sub_bitmap(bmp_t *bmp, int x, int y, int w, int h) {
    bmp_t *sb = malloc(sizeof(bmp_t));
    int bw = BMP_WIDTH(bmp), bh = BMP_HEIGHT(bmp);
    int data_bytes;
    int yi;

    /* Check that the parameters are legal, i.e. the sub-bitmap 
     * lies within the parent */

    if (x < 0 || y < 0 || x + w > bw || y + h > bh) {
        printf("Illegal sub-bitmap specified\n");
        return NULL;
    }

    /* Calculate the data size */
    if (bmp->info_header.bits_per_pixel <= 8)
        data_bytes = h * BYTE_PAD_WORD(BIT_PAD_BYTE(w) / bmp->info_header.bits_per_pixel);
    else
        data_bytes = h * BYTE_PAD_WORD(w * bmp->info_header.bits_per_pixel / 8);

    /* Adjust the file header */
    sb->file_header.filesize = 14 + sizeof(bmp_info_header_t) + 
        sizeof(color_t) * bmp->palette.num_colors + data_bytes;
    sb->file_header.offset = bmp->file_header.offset;

    /* Adjust the info header */
    sb->info_header.header_size = bmp->info_header.header_size;
    sb->info_header.width = w;
    sb->info_header.height = h;
    sb->info_header.num_planes = bmp->info_header.num_planes;
    sb->info_header.bits_per_pixel = bmp->info_header.bits_per_pixel;
    sb->info_header.compression_type = bmp->info_header.compression_type;
    sb->info_header.image_size = data_bytes;
    sb->info_header.x_pixels_per_meter = bmp->info_header.x_pixels_per_meter;
    sb->info_header.y_pixels_per_meter = bmp->info_header.y_pixels_per_meter;
    sb->info_header.colors_used = bmp->info_header.colors_used;
    sb->info_header.colors_important = bmp->info_header.colors_important;

    /* Copy the palette */
    sb->palette.num_colors = bmp->palette.num_colors;
    sb->palette.colors = malloc(sizeof(color_t) * sb->palette.num_colors);
    memcpy(sb->palette.colors, bmp->palette.colors, 
           sizeof(color_t) * sb->palette.num_colors);

    /* Fill in the pixel data */
    sb->pixels = malloc(sizeof(color_t) * w * h);

    for (yi = y; yi < y + h; yi++)
        memcpy(sb->pixels + (yi - y) * w, bmp->pixels + yi * bw + x, 
               sizeof(color_t) * w);

    return sb;
}

/* Paste two bitmaps together horizontally */
bmp_t *merge_bmps(bmp_t *b1, bmp_t *b2) {
    /* Make a 24-bit bitmap for simplicity */
    bmp_t *bmerge = malloc(sizeof(bmp_t));
    int w1 = BMP_WIDTH(b1), w2 = BMP_WIDTH(b2);
    int h1 = BMP_HEIGHT(b1), h2 = BMP_HEIGHT(b2);
    int w = w1 + w2, h = MAX(h1, h2);

    int data_bytes = h * BYTE_PAD_WORD(3 * w);
    int y;

    /* Setup the file header */
    bmerge->file_header.filesize = 14 + sizeof(bmp_info_header_t) + data_bytes;
    bmerge->file_header.offset = 14 + sizeof(bmp_info_header_t);

    /* Setup the info header */
    bmerge->info_header.header_size = sizeof(bmp_info_header_t);
    bmerge->info_header.width = w;
    bmerge->info_header.height = h;
    bmerge->info_header.num_planes = 1;
    bmerge->info_header.bits_per_pixel = 24;
    bmerge->info_header.compression_type = 0;
    bmerge->info_header.image_size = data_bytes;
    bmerge->info_header.x_pixels_per_meter = b1->info_header.x_pixels_per_meter;
    bmerge->info_header.y_pixels_per_meter = b1->info_header.y_pixels_per_meter;
    bmerge->info_header.colors_used = 0;
    // b1->info_header.colors_used; // + b2->info_header.colors_used;
    bmerge->info_header.colors_important = 0;
    // b1->info_header.colors_important; // + b2->info_header.colors_important;

    bmerge->palette.num_colors = 0;
    bmerge->palette.colors = NULL;
    
    /* Fill in the pixel data */
    bmerge->pixels = malloc(sizeof(color_t) * w * h);
    for (y = 0; y < h; y++) {
        if (y >= h1)
            memset(bmerge->pixels + y * w, 0xff, sizeof(color_t) * w1);
        else
            memcpy(bmerge->pixels + y * w, b1->pixels + y * w1, sizeof(color_t) * w1);

        if (y >= h2)
            memset(bmerge->pixels + y * w + w1, 0xff, sizeof(color_t) * w2);
        else
            memcpy(bmerge->pixels + y * w + w1, b2->pixels + y * w2, sizeof(color_t) * w2);
    }
    
    return bmerge;
}

void free_bmp(bmp_t *bmp) {     
    if (bmp->info_header.bits_per_pixel <= 8)
        free(bmp->palette.colors);
    free(bmp->pixels);
    free(bmp);
}
