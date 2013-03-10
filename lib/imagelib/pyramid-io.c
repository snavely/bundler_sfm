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

/* pyramid-io.c */
/* Routines for reading / writing / rendering pyramids */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dmap-io.h"
#include "fileio.h"
#include "image.h"
#include "pyramid.h"
#include "pyramid-io.h"

#define PYRAMID_FILE_HEADER "DPYR"

img_dist_pyr_t *img_dist_pyr_new(int num_levels, int w, int h) {
    img_dist_pyr_t *dpyr = (img_dist_pyr_t *)malloc(sizeof(img_dist_pyr_t));
    int i;
    
    dpyr->num_levels = num_levels;
    dpyr->w = w;
    dpyr->h = h;
    dpyr->em = EM_UNKNOWN;
    dpyr->zweight = 0.0;
    dpyr->nhood_radius = 0;

    dpyr->dmaps = (img_dmap_t *)malloc(sizeof(img_dmap_t) * num_levels);

    for (i = 0; i < num_levels; i++) {
	dpyr->dmaps[i].dists = (double *)malloc(sizeof(double) * (w >> i) * (h >> i));
	dpyr->dmaps[i].nns = (v2_t *)malloc(sizeof(v2_t) * (w >> i) * (h >> i));
    }

    return dpyr;
}

/* Render all images in the pyramid into a single image */
img_t *img_pyr_render(img_pyr_t *pyr) {
    img_t *img = NULL;
    int i;

    if (pyr->num_levels == 1)
	return img_copy(pyr->imgs[0]);

    img_light_invalid_pixels(pyr->imgs[0]);
    img_light_invalid_pixels(pyr->imgs[1]);

    img = img_merge(pyr->imgs[0], pyr->imgs[1]);
    
    for (i = 2; i < pyr->num_levels; i++) {
	img_t *tmp = img;

	img_light_invalid_pixels(pyr->imgs[i]);
	img = img_merge(img, pyr->imgs[i]);
	img_free(tmp);
    }

    return img;
}

/* Free a gaussian pyramid */
void img_pyr_free(img_pyr_t *pyr) {
    int i;
    
    for (i = 0; i < pyr->num_levels; i++)
	img_free(pyr->imgs[i]);
    
    free(pyr->imgs);
    free(pyr);
}

/* Convert a distance map to a distance pyramid */
img_dist_pyr_t *dmap2dpyr(img_dmap_t *dmap) 
{
    img_dist_pyr_t *dpyr = img_dist_pyr_new(1, dmap->w, dmap->h);    

    dpyr->dmaps[0].w = dmap->w;
    dpyr->dmaps[0].h = dmap->h;

    memcpy(dpyr->dmaps[0].dists, dmap->dists, sizeof(double) * dmap->w * dmap->h);
    memcpy(dpyr->dmaps[0].nns, dmap->nns, sizeof(v2_t) * dmap->w * dmap->h);
    dpyr->dmaps[0].uppers = NULL;

    return dpyr;
}

/* Write a distance pyramid to a file */
void img_write_distance_pyramid(FILE *f, img_dist_pyr_t *pyr) {
    int count;
    u_int16_t i, em;

    /* Write the identifier */
    write_word((u_int32_t *)PYRAMID_FILE_HEADER, f);

    /* Write the width and height */
    write_short(&pyr->w, f);
    write_short(&pyr->h, f);

    /* Write the number of levels */
    write_short(&pyr->num_levels, f);

    /* Write information about the metric */
    em = (unsigned short)pyr->em;
    write_short(&em, f);
    write_short(&pyr->nhood_radius, f);
    write_double(&pyr->zweight, f);

    for (i = 0; i < pyr->num_levels; i++) {
	// int num_entries = (pyr->w >> i) * (pyr->h >> i);
	int num_entries = pyr->dmaps[i].w * pyr->dmaps[i].h;
	short int has_uppers;

	/* Write the width and height of this level */
	write_short(&(pyr->dmaps[i].w), f);
	write_short(&(pyr->dmaps[i].h), f);

	if (pyr->dmaps[i].uppers == NULL)
	    has_uppers = 0;
	else
	    has_uppers = 1;

	write_short(&has_uppers, f);

	/* Write the distances */
	for (count = 0; count < num_entries; count++)
	    write_double(&(pyr->dmaps[i].dists[count]), f);

	/* Write the nearest neighbors */
	for (count = 0; count < num_entries; count++) {
	    write_double(&(Vx(pyr->dmaps[i].nns[count])), f);
	    write_double(&(Vy(pyr->dmaps[i].nns[count])), f);
	}

	if (has_uppers) {
	    for (count = 0; count < num_entries; count++) {
		write_short(&(Vx(pyr->dmaps[i].uppers[count])), f);
		write_short(&(Vy(pyr->dmaps[i].uppers[count])), f);
	    }
	}
    }
}

/* Write a distance pyramid to a file */
void img_write_distance_pyramid_file(img_dist_pyr_t *pyr, char *fname) {
    FILE *f = open_file(fname, "w");
    img_write_distance_pyramid(f, pyr);
    fclose(f);
}

/* Read a distance pyramid from a file */
img_dist_pyr_t *img_read_distance_pyramid_file(char *fname) {
    FILE *f = fopen(fname, "r");
    img_dist_pyr_t *pyr;

    if (f == NULL) {
	printf("[img_read_distance_pyramid_file] Error reading file %s\n", fname);
	return NULL;
    }
    
    pyr = img_read_distance_pyramid(f);

    fclose(f);
    
    return pyr;
}

/* Read a distance pyramid from a file */
img_dist_pyr_t *img_read_distance_pyramid(FILE *f) {
    int count;
    u_int16_t i, w, h, num_levels, em;
    char id[5];
    img_dist_pyr_t *pyr;

    /* Read the identifier */
    read_word((u_int32_t *)id, f);
    id[4] = 0;

    if (strcmp(id, PYRAMID_FILE_HEADER) != 0) {
	printf("[img_read_distance_pyramid] Invalid distance pyramid file\n");
	return NULL;
    }

    pyr = (img_dist_pyr_t *)malloc(sizeof(img_dist_pyr_t));

    /* Read the width and height */
    read_short(&w, f);
    read_short(&h, f);

    pyr->w = w;
    pyr->h = h;

    /* Read the number of levels */
    read_short(&num_levels, f);
    pyr->num_levels = num_levels;

    pyr->dmaps = (img_dmap_t *)malloc(sizeof(img_dmap_t) * num_levels);

    /* Read information about the metric */
    read_short(&em, f);
    read_short(&pyr->nhood_radius, f);
    read_double(&pyr->zweight, f);

    pyr->em = em;

    for (i = 0; i < num_levels; i++) {
	u_int16_t w_map, h_map;
	
	// int num_entries = (w >> i) * (h >> i);
	int num_entries;
	short has_uppers;

	read_short(&w_map, f);
	read_short(&h_map, f);

	pyr->dmaps[i].w = w_map; // (w >> i);
	pyr->dmaps[i].h = h_map; // (h >> i);

	num_entries = pyr->dmaps[i].w * pyr->dmaps[i].h;

	/* Initialize the maps */
	pyr->dmaps[i].dists = (double *)malloc(sizeof(double) * num_entries);
	pyr->dmaps[i].nns = (v2_t *)malloc(sizeof(v2_t) * num_entries);
	pyr->dmaps[i].uppers = NULL;

	read_short(&has_uppers, f);
	if (has_uppers)
	    pyr->dmaps[i].uppers = (iv2_t *)malloc(sizeof(iv2_t) * num_entries);

	/* Read the distances */
	for (count = 0; count < num_entries; count++) {
	    read_double(&(pyr->dmaps[i].dists[count]), f);
	}
    
	for (count = 0; count < num_entries; count++) {
	    read_double(&Vx(pyr->dmaps[i].nns[count]), f);
	    read_double(&Vy(pyr->dmaps[i].nns[count]), f);
	}

	if (has_uppers) {
	    for (count = 0; count < num_entries; count++) {
		read_short(&Vx(pyr->dmaps[i].uppers[count]), f);
		read_short(&Vy(pyr->dmaps[i].uppers[count]), f);
	    }
	}
    }

    return pyr;
}

/* Free a distance pyramid */
void img_free_distance_pyramid(img_dist_pyr_t *pyr) {
    int i;
    
    for (i = 0; i < pyr->num_levels; i++) {
	free(pyr->dmaps[i].dists);
	free(pyr->dmaps[i].nns);

	if (pyr->dmaps[i].uppers != NULL)
	    free(pyr->dmaps[i].uppers);
    }

    free(pyr->dmaps);

    // free(pyr->dists);
    // free(pyr->nns);

    free(pyr);
}

/* Render the distance pyramid as an image */
img_t *img_distance_pyramid_render(img_dist_pyr_t *pyr, int mode) {
    img_t *img = NULL, *dmap0 = NULL, *dmap1 = NULL;
    int i;

    if (pyr->num_levels == 1) {
	if (mode == 0)
	    return img_dmap_render(&(pyr->dmaps[0]));
	if (mode == 1) 
	    return img_dmap_render_disparity(&(pyr->dmaps[0]));
	if (mode == 2) 
	    return img_dmap_render_flow(&(pyr->dmaps[0]));
    }
    
    if (mode == 0) {
	dmap0 = img_dmap_render(&(pyr->dmaps[0]));
	dmap1 = img_dmap_render(&(pyr->dmaps[1]));
    } else if (mode == 1) {
	dmap0 = img_dmap_render_disparity(&(pyr->dmaps[0]));
	dmap1 = img_dmap_render_disparity(&(pyr->dmaps[1]));
    } else if (mode == 2) {
	dmap0 = img_dmap_render_flow(&(pyr->dmaps[0]));
	dmap1 = img_dmap_render_flow(&(pyr->dmaps[1]));
    }
    
    img = img_merge(dmap0, dmap1);

    img_free(dmap0);
    img_free(dmap1);
    
    for (i = 2; i < pyr->num_levels; i++) {
	img_t *tmp = img;
	img_t *dmap = NULL;

	if (mode == 0) {
	    dmap = img_dmap_render(&(pyr->dmaps[i]));
	} else if (mode == 1) {
	    dmap = img_dmap_render_disparity(&(pyr->dmaps[i]));
	} else if (mode == 2) {
	    dmap = img_dmap_render_flow(&(pyr->dmaps[i]));
	}
	
	img = img_merge(img, dmap);
	img_free(dmap);
	img_free(tmp);
    }

    return img;
}
