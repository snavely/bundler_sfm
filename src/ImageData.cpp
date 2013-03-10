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

/* ImageData.cpp */
/* Simple image storage and operations */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ImageData.h"
#include "SifterUtil.h"
#include "Register.h"

#ifndef __BUNDLER__
#include "UtilGL.h"
#endif /* __BUNDLER__ */

#include "canny.h"
#include "defines.h"
#include "filter.h"
#include "fit.h"
#include "matrix.h"
#include "pgm.h"
#include "resample.h"
#include "transform.h"
#include "util.h"

#include "LoadJPEG.h"

#define SUBSAMPLE_LEVEL 1

#if 1
#define FOCAL_LENGTH 490
#define NUM_ROTATIONS 1
#define DEFAULT_WIDTH 1536
#define DEFAULT_HEIGHT 1024
#else
#define FOCAL_LENGTH 700 // 661.62760416666674
#define NUM_ROTATIONS 0
#define DEFAULT_WIDTH 2074
#define DEFAULT_HEIGHT 1382
#endif

extern double global_scale;

int CompareImageDates(ImageDate *date1, ImageDate *date2)
{
    /* Check that both dates are known */
    if (!date1->m_known || !date2->m_known)
        return 0;

    if (date1->m_year < date2->m_year)
        return -1;
    else if (date2->m_year < date1->m_year)
        return -1;

    if (date1->m_month < date2->m_month)
        return -1;
    else if (date2->m_month < date1->m_month)
        return 1;

    if (date1->m_day < date2->m_day)
        return -1;
    else if (date2->m_day < date1->m_day)
        return 1;    

    return 0;
}

static char *month_strings[] = 
    { "January",
      "February",
      "March",
      "April",
      "May",
      "June",
      "July",
      "August",
      "September",
      "October",
      "November",
      "December" };

void ImageDate::GetMonthString(char *buf)
{
    if (!m_known || m_month < 1 || m_month > 12)
        strcpy(buf, "Unknown");
    else
        strcpy(buf, month_strings[m_month - 1]);
}

void ImageDate::GetString(char *buf)
{
    if (!m_known) {
        strcpy(buf, "Unknown");
        return;
    }

    if (m_month == 0) {
        sprintf(buf, "%d", m_year);
    } else {
        char month_str[256];
        GetMonthString(month_str);

        sprintf(buf, "%s %d, %d, %02d:%02d", 
            month_str, m_day, m_year, m_hour, m_minute);
    }
}

double ImageDate::GetDateDouble() 
{
    if (!m_known)
        return DBL_MAX;

    /* Return the date in the form of a double */
    return (double) m_second + 
        60.0 * (m_minute + 
        60.0 * (m_hour + 
        24.0 * (m_day + 
        31.0 * (m_month + 
        12.0 * (m_year - 1800)))));	
}

double ImageDate::GetDateOnlyDouble()
{
    if (!m_known)
        return DBL_MAX;

    /* Return the date in the form of a double */
    return (m_day + 31.0 * (m_month + 12.0 * (m_year - 1800)));
}

double ImageDate::GetTimeDouble() 
{
    if (!m_known)
        return DBL_MAX;

    /* Return the date in the form of a double */
    return (double) m_second + 60.0 * (m_minute + 60.0 * (m_hour));
}

void ImageNote::Read(FILE *f, double width, double height)
{
    char buf[256];
    fgets(buf, 255, f);

    if (buf[strlen(buf) - 1] == '\n')
        buf[strlen(buf) - 1] = 0;

    if (buf[strlen(buf) - 1] == '\r')
        buf[strlen(buf) - 1] = 0;

    for (int i = 0; i < (int) strlen(buf); i++) {
        if (buf[i] == '*') buf[i] = '\n';
    }

    m_text = strdup(buf);

    double x, y, w, h;
    fscanf(f, "%lf %lf %lf %lf\n", &x, &y, &w, &h);

    x -= 0.5 * width;
    y = height - y - 1 - 0.5 * height;

    m_bbox = BoundingBox(x, y - h, x + w, y);
}

/* Initialize the image data given a string description */
void ImageData::InitFromString(char *buf, char *path, bool fisheye_by_default) 
{
    /* Eat the newline */
    if (buf[strlen(buf)-1] == '\n')
        buf[strlen(buf)-1] = 0;

    if (buf[strlen(buf)-1] == '\r')
        buf[strlen(buf)-1] = 0;

    /* Split the buffer into tokens */
    std::string str(buf);
    std::vector<std::string> toks;
    Tokenize(str, toks, " ");

#if 0
    while (t.HasMoreTokens()) {
        tok = t.GetNextToken();
        toks.push_back(tok);
    }
#endif

    int num_toks = (int) toks.size();

    bool fisheye = fisheye_by_default;
    if (num_toks >= 2) {
        fisheye = (atoi(toks[1].c_str()) == 1);
    }

    bool has_init_focal = false;
    double init_focal = 0.0;
    if (num_toks >= 3) {
        has_init_focal = true;
        init_focal = atof(toks[2].c_str());
    }

    ImageData data;

    // printf("[ImageData::InitFromFile] Adding image %s\n", toks[0].c_str());

    if (path == NULL || strcmp(path, ".") == 0 || toks[0].c_str()[0] == '/') {
        m_name = strdup(toks[0].c_str());
    } else {
        char tmp_name[512];
        sprintf(tmp_name, "%s/%s", path, toks[0].c_str());
        m_name = strdup(tmp_name);
    }

    m_img = NULL;
    m_thumb = NULL;
    m_thumb8 = NULL;


    // m_wximage = NULL;
    m_image_loaded = false;
    m_keys_loaded = false;
    m_keys_scale_rot_loaded = false;

    m_fisheye = fisheye;
    m_has_init_focal = has_init_focal;
    m_init_focal = init_focal;
    m_camera.m_adjusted = false;
    m_texture_index = -1;

    /* Extract out the user name */
    char base_name[256];
    GetBaseName(base_name);

    char *name_end = base_name;
    // char user[256];

    while (*name_end != '_' && *name_end != 0)
        name_end++;

    if (*name_end == 0) {
        strcpy(m_user_name, "Unknown");
        strcpy(m_flickr_index, "Unknown");
    } else {
        strncpy(m_user_name, base_name, name_end - base_name);
        m_user_name[name_end - base_name] = 0;

        strcpy(m_flickr_index, name_end + 1);

        // for (int i = 0; m_user_name[i] != 0; i++) 
        //    m_user_name[i] = (m_user_name[i] - 'a' + 10) % 26 + 'a';
    }

    /* Try to find a keypoint file */
    char key_buf[256];
    strcpy(key_buf, m_name);
    key_buf[strlen(m_name) - 3] = 'k';
    key_buf[strlen(m_name) - 2] = 'e';
    key_buf[strlen(m_name) - 1] = 'y';

    m_key_name = strdup(key_buf);
}

void ImageData::LoadImage() {
    /* Check if there is a jpg file with the same basename */
    char jpeg_buf[256];
    strcpy(jpeg_buf, m_name);
    jpeg_buf[strlen(m_name) - 3] = 'j';
    jpeg_buf[strlen(m_name) - 2] = 'p';
    jpeg_buf[strlen(m_name) - 1] = 'g';

    /* Check if there is a bmp file with the same basename */
    char bmp_buf[256];
    strcpy(bmp_buf, m_name);
    bmp_buf[strlen(m_name) - 3] = 'b';
    bmp_buf[strlen(m_name) - 2] = 'm';
    bmp_buf[strlen(m_name) - 1] = 'p';

    img_t *img;

    if (FileExists(jpeg_buf)) {
        // printf("[ImageData::LoadImage] Reading JPEG...\n");
        img = LoadJPEG(jpeg_buf);

        // printf("[ImageData::LoadImage] Dimensions: %d x %d\n", 
        //        img->w, img->h);
    } else if (FileExists(bmp_buf)) {
        img = img_read_bmp_file(bmp_buf);
    } else {
        img = img_read_pgm_file(m_name);
    }

    m_img = img;
    m_image_loaded = true;
}

void ImageData::UnloadImage() 
{
    if (!m_image_loaded) {
        printf("[ImageData::UnloadImage] Image hasn't been loaded!\n");
        return;
    }

    img_free(m_img);

    m_image_loaded = false;
}


void ImageData::LoadBackTexImage()
{

	/* Check if there is a texture image */
	char texture_jpeg_buf[256];
	strcpy(texture_jpeg_buf, m_name);
	texture_jpeg_buf[strlen(m_name) - 3] = 't';
	texture_jpeg_buf[strlen(m_name) - 2] = 'e';
	texture_jpeg_buf[strlen(m_name) - 1] = 'x';
	texture_jpeg_buf[strlen(m_name) + 0] = 't';
	texture_jpeg_buf[strlen(m_name) + 1] = 'u';
	texture_jpeg_buf[strlen(m_name) + 2] = 'r';
	texture_jpeg_buf[strlen(m_name) + 3] = 'e';
	texture_jpeg_buf[strlen(m_name) + 4] = '3';
	texture_jpeg_buf[strlen(m_name) + 5] = '.';
	texture_jpeg_buf[strlen(m_name) + 6] = 'j';
	texture_jpeg_buf[strlen(m_name) + 7] = 'p';
	texture_jpeg_buf[strlen(m_name) + 8] = 'g';
	texture_jpeg_buf[strlen(m_name) + 9] = 0;
	/* Check if inpainted texture is here */
	char jpeg_buf[256];
	strcpy(jpeg_buf, m_name);
	jpeg_buf[strlen(m_name) - 3] = 'j';
	jpeg_buf[strlen(m_name) - 2] = 'p';
	jpeg_buf[strlen(m_name) - 1] = 'g';

	/* Otheriwse use the original image*/
	char bmp_buf[256];
	strcpy(bmp_buf, texture_jpeg_buf);
	bmp_buf[strlen(texture_jpeg_buf) - 5] = '.';
	bmp_buf[strlen(texture_jpeg_buf) - 4] = 'j';
	bmp_buf[strlen(texture_jpeg_buf) - 3] = 'p';
	bmp_buf[strlen(texture_jpeg_buf) - 2] = 'g';
	bmp_buf[strlen(texture_jpeg_buf) - 1] =0;

	img_t *img;

	if (FileExists(texture_jpeg_buf)) {  
            // printf("[ImageData::LoadImage] Reading JPEG...\n");
            img = LoadJPEG(texture_jpeg_buf);

            // printf("[ImageData::LoadImage] Dimensions: %d x %d\n", 
            //        img->w, img->h);
	} else if (FileExists(bmp_buf)) {
            img = LoadJPEG(bmp_buf);
	}
	else if (FileExists(jpeg_buf)) {
            img = LoadJPEG(jpeg_buf);
	} 
	else {
            img = img_read_pgm_file(m_name);
	}

	m_back_texture_img = img;
	m_texture_image_loaded = true;
}



void ImageData::LoadTexImage()
{
    /* Check if there is a texture image */
    char texture_jpeg_buf[256];
    strcpy(texture_jpeg_buf, m_name);
    texture_jpeg_buf[strlen(m_name) - 3] = 't';
    texture_jpeg_buf[strlen(m_name) - 2] = 'e';
    texture_jpeg_buf[strlen(m_name) - 1] = 'x';
    texture_jpeg_buf[strlen(m_name) + 0] = 't';
    texture_jpeg_buf[strlen(m_name) + 1] = 'u';
    texture_jpeg_buf[strlen(m_name) + 2] = 'r';
    texture_jpeg_buf[strlen(m_name) + 3] = 'e';
	//texture_jpeg_buf[strlen(m_name) + 4] = '1';
	texture_jpeg_buf[strlen(m_name) + 4] = '.';
	texture_jpeg_buf[strlen(m_name) + 5] = 'j';
	texture_jpeg_buf[strlen(m_name) + 6] = 'p';
	texture_jpeg_buf[strlen(m_name) + 7] = 'g';
	texture_jpeg_buf[strlen(m_name) + 8] = 0;

    /* Check if there is a jpg file with the same basename */
    char jpeg_buf[256];
    strcpy(jpeg_buf, m_name);
    jpeg_buf[strlen(m_name) - 3] = 'j';
    jpeg_buf[strlen(m_name) - 2] = 'p';
    jpeg_buf[strlen(m_name) - 1] = 'g';

    /* Check if there is a bmp file with the same basename */
    char bmp_buf[256];
    strcpy(bmp_buf, m_name);
    bmp_buf[strlen(m_name) - 3] = 'b';
    bmp_buf[strlen(m_name) - 2] = 'm';
    bmp_buf[strlen(m_name) - 1] = 'p';

    img_t *img;

    if (FileExists(texture_jpeg_buf)) {
	// printf("[ImageData::LoadImage] Reading JPEG...\n");
	img = LoadJPEG(texture_jpeg_buf);

	// printf("[ImageData::LoadImage] Dimensions: %d x %d\n", 
        //        img->w, img->h);
    } else if (FileExists(jpeg_buf)) {
	img = LoadJPEG(jpeg_buf);
    } else if (FileExists(bmp_buf)) {
	img = img_read_bmp_file(bmp_buf);
    } else {
	img = img_read_pgm_file(m_name);
    }
    
    m_texture_img = img;
    m_texture_image_loaded = true;
}


void ImageData::UnloadBackTexImage()
{
	//if (!m_texture_image_loaded) {
	//	printf("[ImageData::UnloadTextureImage] Image hasn't been loaded!\n");
	//	return;
	//}

	img_free(m_back_texture_img);

	m_texture_image_loaded = false;    
}


void ImageData::UnloadTexImage()
{
    if (!m_texture_image_loaded) {
        printf("[ImageData::UnloadTextureImage] Image hasn't been loaded!\n");
        return;
    }

    img_free(m_texture_img);

    m_texture_image_loaded = false;    
}

bool ImageData::TexImageExists()
{
    /* Check if there is a texture image */
    char texture_jpeg_buf[256];
    strcpy(texture_jpeg_buf, m_name);
    texture_jpeg_buf[strlen(m_name) - 3] = 't';
    texture_jpeg_buf[strlen(m_name) - 2] = 'e';
    texture_jpeg_buf[strlen(m_name) - 1] = 'x';
    texture_jpeg_buf[strlen(m_name) + 0] = 't';
    texture_jpeg_buf[strlen(m_name) + 1] = 'u';
    texture_jpeg_buf[strlen(m_name) + 2] = 'r';
    texture_jpeg_buf[strlen(m_name) + 3] = 'e';

    texture_jpeg_buf[strlen(m_name) + 4] = '.';
    texture_jpeg_buf[strlen(m_name) + 5] = 'j';
    texture_jpeg_buf[strlen(m_name) + 6] = 'p';
    texture_jpeg_buf[strlen(m_name) + 7] = 'g';
    texture_jpeg_buf[strlen(m_name) + 8] = 0;
    
    if (FileExists(texture_jpeg_buf)) {
        return true;
    } else {
        return false;
    }
}



void ImageData::CheckLoadFloatingThumb()
{
#ifndef __BUNDLER__
    if (m_thumb != NULL)
	return;

    LoadFloatingThumbnail();
    LoadTextureImage(m_thumb, m_thumb_texture_index, m_thumb_bounds);
#endif /* __BUNDLER__ */
}

void ImageData::LoadFloatingThumbnail() {
    printf("Loading floating thumb\n");

    /* Look for a cached image */
    char basename[256];
    strcpy(basename, m_name);
    basename[strlen(basename) - 4] = 0;

    char thumb_jpg_buf[256];
    sprintf(thumb_jpg_buf, "%s.thumb.float.jpg", basename);

    char thumb_bmp_buf[256];
    sprintf(thumb_bmp_buf, "%s.thumb.float.bmp", basename);

    m_thumb = NULL;    

    if (FileExists(thumb_jpg_buf)) {
	m_thumb = LoadJPEG(thumb_jpg_buf);
    }

    if (FileExists(thumb_bmp_buf)) {
	m_thumb = img_read_bmp_file(thumb_bmp_buf);
    }

    if (m_thumb != NULL) {
	/* Shrink the thumbnail */
	if (m_thumb->w >= 512 || m_thumb->h >= 512) {
	    int scale = 1;

	    while (m_thumb->w / scale >= 512 || m_thumb->h / scale >= 512)
		scale *= 2;
	    
	    img_t *thumb_new = img_scale_fast(m_thumb, scale);
	    img_free(m_thumb);
	    m_thumb = thumb_new;
	}

	return ;
    }

    img_t *img = UndistortImage(0.0, 0.0);
    img_t *thumb = img_scale_fast(img, 4);

    m_thumb = thumb;

    img_free(img);

    /* Cache the image */
    img_write_bmp_file(m_thumb, thumb_bmp_buf);
}

void ImageData::LoadThumb256() {
    /* Look for a cached image */
    char basename[256];
    strcpy(basename, m_name);
    basename[strlen(basename) - 4] = 0;

    char thumb_jpg_buf[256];
    sprintf(thumb_jpg_buf, "%s.thumb256.jpg", basename);

    char thumb_bmp_buf[256];
    sprintf(thumb_bmp_buf, "%s.thumb256.bmp", basename);

    m_thumb256 = NULL;

    if (FileExists(thumb_jpg_buf)) {
	m_thumb256 = LoadJPEG(thumb_jpg_buf);
    }

    if (FileExists(thumb_bmp_buf)) {
	m_thumb256 = img_read_bmp_file(thumb_bmp_buf);
    }

    if (m_thumb256 != NULL) {
	return;
    }

    // img_t *img = UndistortImage(0.0, 0.0);
    LoadImage();

    double scale;
    img_t *thumb256 = RescaleImage(m_img, 256, scale);

    m_thumb256 = thumb256;

    UnloadImage();

    /* Cache the image */
    img_write_bmp_file(m_thumb256, thumb_bmp_buf);    
}

void ImageData::UnloadFloatingThumbnail() {
    if (m_thumb != NULL)
	img_free(m_thumb);

    m_thumb = NULL;
    m_thumb_texture_index = -1;
}

void ImageData::CheckLoadFixedThumb()
{
#ifndef __BUNDLER__
    if (m_thumb_fixed != NULL)
	return;

    LoadFixedThumbnail(128, 128, m_rotation);
    LoadTextureImage(m_thumb_fixed, m_thumb_fixed_texture_index, 
		     m_thumb_fixed_bounds);
#endif /* __BUNDLER__ */
}

void ImageData::LoadFixedThumbnail(int w_max, int h_max, int rotation) {
    printf("Loading fixed thumb\n");

    /* Look for a cached image */
    char basename[256];
    strcpy(basename, m_name);
    basename[strlen(basename) - 4] = 0;

    char thumb_jpg_buf[256];
    sprintf(thumb_jpg_buf, "%s.thumb.jpg", basename);

    char thumb_bmp_buf[256];
    sprintf(thumb_bmp_buf, "%s.thumb.bmp", basename);

    if (FileExists(thumb_jpg_buf)) {
	m_thumb_fixed = LoadJPEG(thumb_jpg_buf);
	return;
    }

    if (FileExists(thumb_bmp_buf)) {
	m_thumb_fixed = img_read_bmp_file(thumb_bmp_buf);
	return;
    }

    int w = GetWidth();
    int h = GetHeight();

    double w_ratio = (double) w / (double) w_max;
    double h_ratio = (double) h / (double) h_max;

    double ratio;
    if (w_ratio > h_ratio) {
	ratio = w_ratio;
    } else {
	ratio = h_ratio;
    }

    img_t *img = UndistortImage(0.0, 0.0, rotation);

    printf("Blurring, sigma is %0.3f\n", 0.35 * ratio);
    img_t *blur = img_smooth(img, 0.35 * ratio, 0);
    trans2D_t *T = new_scaling_transform(1.0 / ratio, 1.0 / ratio);
    img_t *scaled = img_resample_bbox(blur, T);

    img_free(blur);
    transform_free(T);

    img_t *thumb = img_new(w_max, h_max);
    
    int x_start = iround(0.5 * (w_max - scaled->w));
    int y_start = iround(0.5 * (h_max - scaled->h));
    int x_end = x_start + scaled->w;
    int y_end = y_start + scaled->h;

    for (int y = 0; y < w_max; y++) {
	for (int x = 0; x < h_max; x++) {
	    if (x < x_start || y < y_start) {
		img_set_pixel(thumb, x, y, 0x0, 0x0, 0x0);
	    } else if (x >= x_end || y >= y_end) {
		img_set_pixel(thumb, x, y, 0x0, 0x0, 0x0);
	    } else {
		color_t c = img_get_pixel(scaled, x - x_start, y - y_start);
		img_set_pixel(thumb, x, y, c.r, c.g, c.b);
	    }
	}
    }

    m_thumb_fixed = thumb;

    img_free(img);

    /* Cache the image */
    img_write_bmp_file(m_thumb_fixed, thumb_bmp_buf);
}

void ImageData::UnloadFixedThumbnail() {
    if (m_thumb_fixed != NULL)
	img_free(m_thumb_fixed);

    m_thumb_fixed = NULL;
    m_thumb_fixed_texture_index = -1;
}

int ImageData::GetNumKeys()
{
#ifndef __DEMO__
    if (m_keys_loaded) {
        return (int) m_keys.size();
    } else {
        if (!m_cached_keys)
            CacheNumKeys();
        
        return m_num_keys;
    }
#else
	return 0;
#endif
}

void ImageData::LoadOrExtractKeys(char *sift_binary, bool undistort) 
{
    if (m_keys_loaded)
	return;   /* Already loaded the keys */

#if 0
    /* Try to find a keypoint file */
    char key_buf[256];
    strcpy(key_buf, m_name);
    key_buf[strlen(m_name) - 3] = 'k';
    key_buf[strlen(m_name) - 2] = 'e';
    key_buf[strlen(m_name) - 1] = 'y';
#endif
    char gzKeyName[512];
    sprintf(gzKeyName, "%s.gz", m_key_name);

    if (FileExists(m_key_name) || FileExists(gzKeyName)) {
	LoadKeys(true, undistort);
    } else {
	ExtractFeatures(sift_binary, undistort);
    }
}

void ImageData::ExtractFeatures(char *sift_binary, bool undistort) 
{
#ifndef __DEMO__
    /* Find the extension */
    const char *in = m_name;
    char *ext = strrchr((char *) in, '.');
    char out[256];
	
    strncpy(out, in, ext - in);
    out[ext - in] = 0;
    strcat(out, ".key");

    char cmd[2048];
    sprintf(cmd, "%s < %s > %s", sift_binary, in, out);

#if 0
    if (log) {
	log->AppendText("  ");
	log->AppendText("Running command '");
	log->AppendText(cmd);
	log->AppendText("'\n");
    }
#endif
    
    system(cmd);

    // m_key_names.push_back(wxString(out));
    m_key_name = strdup(out);

    /* Read back the keypoints */
    std::vector<Keypoint> kps = ReadKeyFile(out);

    /* Flip y-axis to make things easier */
    for (int k = 0; k < (int) kps.size(); k++) {
	kps[k].m_y = GetHeight() - kps[k].m_y - 1.0;
    }

    /* Now make the image center the origin */
    for (int k = 0; k < (int) kps.size(); k++) {
	kps[k].m_x -= 0.5 * GetWidth();
	kps[k].m_y -= 0.5 * GetHeight();
    }

    // m_keys.push_back(kps);
    m_keys = kps;

    if (undistort)
        UndistortKeys();    
#else
    printf("Feature extraction unavailable in the demo version\n");
#endif
}

void ImageData::LoadDescriptors(bool undistort) {
    if (!m_keys_desc_loaded) {
        LoadKeys(true, undistort);
        return;
    }
            
#if 0
    /* Keys are loaded, check if descriptors are loaded */
    int num_keys = (int) m_keys.size();
    if (num_keys == 0)
        return;

    if (m_keys[0].m_d != NULL)
        return; /* Descriptors appear to be loaded */

    /* Keys are loaded but not descriptors */
    std::vector<Keypoint> kps = ReadKeyFile(m_key_name, undistort);
    assert(kps.size() == m_keys.size());
    
    for (int i = 0; i < num_keys; i++) {
        m_keys[i].m_d = kps[i].m_d;
    }
#endif
}

void ImageData::LoadKeys(bool descriptor, bool undistort) {
#ifndef __DEMO__
    if (m_keys_loaded && !descriptor)
	return;   /* Already loaded the keys */

    if (m_keys_desc_loaded && descriptor)
        return;   /* Already loaded keys with descriptors */

    /* Try to find a keypoint file */
    if (!descriptor) {
        std::vector<Keypoint> kps = ReadKeyFile(m_key_name);

        /* Flip y-axis to make things easier */
        for (int k = 0; k < (int) kps.size(); k++) {
            kps[k].m_y = GetHeight() - kps[k].m_y - 1.0;
        }
        
        /* Now make the image center the origin */
        for (int k = 0; k < (int) kps.size(); k++) {
            kps[k].m_x -= 0.5 * (GetWidth() - 1);
            kps[k].m_y -= 0.5 * (GetHeight() - 1);
        }
        
        m_keys = kps;
    
        m_keys_loaded = true;
    } else {
        std::vector<KeypointWithDesc> kps = 
            ReadKeyFileWithDesc(m_key_name, true);

        /* Flip y-axis to make things easier */
        for (int k = 0; k < (int) kps.size(); k++) {
            kps[k].m_y = GetHeight() - kps[k].m_y - 1.0;
        }
        
        /* Now make the image center the origin */
        for (int k = 0; k < (int) kps.size(); k++) {
            kps[k].m_x -= 0.5 * GetWidth();
            kps[k].m_y -= 0.5 * GetHeight();
        }
        
        m_keys_desc = kps;
    
        m_keys_desc_loaded = true;            
    }
    
    if (undistort)
        UndistortKeys();
#else
    printf("Cannot load keys in demo version\n");
#endif
}

void ImageData::UnloadKeys() {
    if (m_keys_desc_loaded) {
        int num_keys = (int) m_keys_desc.size();

        for (int i = 0; i < num_keys; i++) {
            if (m_keys_desc[i].m_d)
                delete [] m_keys_desc[i].m_d;

            m_keys_desc[i].m_d = NULL;
        }

        m_keys_desc.clear();
        m_keys_desc_loaded = false;
    } else if (m_keys_loaded) {
        m_keys.clear();
        m_keys_loaded = false;
    }
}

void ImageData::LoadKeysWithScaleRot(bool descriptor, bool undistort) 
{
#ifndef __DEMO__
    if (m_keys_scale_rot_loaded)
	return;   /* Already loaded the keys */

    /* Try to find a keypoint file */    
    std::vector<KeypointWithScaleRot> kps = 
        ReadKeyFileWithScaleRot(m_key_name, descriptor);

    /* Flip y-axis to make things easier */
    for (int k = 0; k < (int) kps.size(); k++) {
	kps[k].m_y = GetHeight() - kps[k].m_y - 1.0;
    }

    /* Now make the image center the origin */
    for (int k = 0; k < (int) kps.size(); k++) {
	kps[k].m_x -= 0.5 * GetWidth();
	kps[k].m_y -= 0.5 * GetHeight();
    }

    m_keys_scale_rot = kps;
    m_keys_scale_rot_loaded = true;

    if (undistort)
        UndistortKeys();
#else
    printf("Cannot load keys in demo version\n");
#endif
}

void ImageData::UnloadKeysWithScaleRot() {
    // int num_keys = (int) m_keys_scale_rot.size();

#if 0
    for (int i = 0; i < num_keys; i++) {
	if (m_keys_scale_rot[i].m_d)
	    delete [] m_keys_scale_rot[i].m_d;
    }
#endif

    m_keys_scale_rot.clear();
    m_keys_scale_rot_loaded = false;
}

void ImageData::CacheNumKeys()
{
#ifndef __DEMO__
    m_num_keys = GetNumberOfKeys(m_key_name);
#else
    m_num_keys = 0;
#endif

    m_cached_keys = true;
}

/* Save the dimensions of the image */
void ImageData::CacheDimensions() 
{
    int w = 1504, h = 1000;

    FILE *f = fopen(m_name, "r");

    if (f != NULL) {
	fclose(f);

	  int len = strlen(m_name);

        if (strcmp(m_name + len - 3, "pgm") == 0)
            img_read_pgm_dimensions(m_name, &w, &h);
        else if (strcmp(m_name + len - 3, "bmp") == 0)
	    bmp_file_get_dimensions(m_name, &w, &h);
        else if (strcmp(m_name + len - 3, "jpg") == 0)
            GetJPEGDimensions(m_name, w, h);
    } else {
	/* Create the bmp file */
	char bmp_file[256];
	strcpy(bmp_file, m_name);

	bmp_file[strlen(m_name) - 3] = 'b';
	bmp_file[strlen(m_name) - 2] = 'm';
	bmp_file[strlen(m_name) - 1] = 'p';

	if (FileExists(bmp_file)) {
	    bmp_file_get_dimensions(bmp_file, &w, &h);
	} else {
	    char jpeg_file[256];
	    strcpy(jpeg_file, m_name);
	    jpeg_file[strlen(m_name) - 3] = 'j';
	    jpeg_file[strlen(m_name) - 2] = 'p';
	    jpeg_file[strlen(m_name) - 1] = 'g';
	    
            if (FileExists(jpeg_file)) {
                GetJPEGDimensions(jpeg_file, w, h);
            } else {
                printf("[ImageData::CacheDimensions] Fatal error: "
                       "couldn't read image %s\n", m_name);
            }
	}
    }

    m_width = w;
    m_height = h;
    m_cached_dimensions = true;
}

int ImageData::GetWidth() { 
    if (m_image_loaded) return iround(global_scale * m_img->w);
    else { 
	if (m_cached_dimensions)
	    return iround(global_scale * m_width);

	CacheDimensions();
	return iround(global_scale * m_width);
    }

    // return DEFAULT_WIDTH; /* hack */
}
int ImageData::GetHeight() { 
    if (m_image_loaded) return iround(global_scale * m_img->h); 
    else { 
	if (m_cached_dimensions)
	    return iround(global_scale * m_height);

	CacheDimensions();
	return iround(global_scale * m_height);
    }
    // return DEFAULT_HEIGHT; /* hack */
}

int ImageData::GetArea() {
    return GetWidth() * GetHeight();
}

void ImageData::GetBaseName(char *buf) {
    int len = strlen(m_name);
    
    int i;
    for (i = len - 1; i >= 0 && m_name[i] != '/'; i--) { }

    if (i == -1)
	strcpy(buf, m_name);
    else
	strcpy(buf, m_name + i + 1);
    
    buf[strlen(buf) - 4] = 0;
}

void ImageData::GetName(char *buf)
{
    if (m_real_name != NULL) {
	strcpy(buf, m_real_name);
    } else {
	GetBaseName(buf);
	// for (int i = 0; buf[i] != '_' && buf[i] != 0; i++)
	//     buf[i] = (buf[i] - 'a' + 10) % 26 + 'a';
    }
}

void ImageData::GetCompressedTextureFilename(char *buf) const
{
    strcpy(buf, m_name);
    buf[strlen(m_name) - 3] = 'c';
    buf[strlen(m_name) - 2] = 't';
    buf[strlen(m_name) - 1] = 'x';
}

bool ImageData::CompressedTextureExists()
{
    char buf[256];
    GetCompressedTextureFilename(buf);

    if (FileExists(buf))
        return true;
    else
        return false;
}

/* Return the bounding box of the image */
BoundingBox ImageData::GetBoundingBox()
{
    int w = GetWidth(), h = GetHeight();
    return BoundingBox(-0.5 * w, -0.5 * h, 0.5 * w, 0.5 * h);
}

/* Get the texture coordinates for a point */
v2_t ImageData::GetPointTexCoords(v2_t point)
{
    int w = GetWidth();
    int h = GetHeight();

    return v2_new((Vx(point) + 0.5 * w) / (double) w,
		  (Vy(point) + 0.5 * h) / (double) h);
}

/* Get the texture coordinates of a set of points */
std::vector<v2_t> ImageData::GetPointTexCoords(const std::vector<v2_t> &pts)
{
    int num_points = (int) pts.size();
    std::vector<v2_t> coords;

    for (int i = 0; i < num_points; i++) {
	coords.push_back(GetPointTexCoords(pts[i]));
    }

    return coords;
}



/* Returns true if the given pixel is in range */
bool ImageData::PixelInRange(double x, double y) {
    int w = GetWidth();
    int h = GetHeight();

    if (x < 0.0 || x >= w - 1 || y < 0.0 || y >= h - 1)
        return false;

    return true;
}


/* Convert to/from normalized device coordinates */
void ImageData::ToNDC(double x, double y, double &x_n, double &y_n) {
    double size= MAX(GetWidth(), GetHeight());
    double size2   = 0.5 * size;
    double size2i  = 1.0 / size2;

    x_n = x * size2i;
    y_n = y * size2i;
}

void ImageData::FromNDC(double x_n, double y_n, double &x, double &y) {
    double size= MAX(GetWidth(), GetHeight());
    double size2   = 0.5 * size;

    x = x_n * size2;
    y = y_n * size2;
}

void ImageData::DistortPoint(double x, double y, 
                             double &x_out, double &y_out) const
{
    double I[9];
    GetRotationFromSpherical(-0.5 * M_PI, 0.5 * M_PI, I);
    
    DistortPoint(x, y, I, x_out, y_out);
}

void ImageData::DistortPoint(double x, double y, double *R,
			     double &x_out, double &y_out) const
{
    if (!m_fisheye) {
	// DistortPointRD(x, y, x_out, y_out);
	x_out = x;
	y_out = y;
	return;
    }

    double xn = x; // - m_fCx;
    double yn = y; // - m_fCy;

    double ray[3] = { xn, yn, -m_fFocal }, ray_rot[3];
    matrix_product(3, 3, 3, 1, R, ray, ray_rot);

    if (ray_rot[2] <= 0.0) {
	xn = -DBL_MAX;
	yn = -DBL_MAX;
	return;
    } else {
	xn = ray_rot[0] * m_fFocal / ray_rot[2];
	yn = ray_rot[1] * m_fFocal / ray_rot[2];
    }
    
    double r = sqrt(xn * xn + yn * yn);
    double angle = RAD2DEG(atan(r / m_fFocal));
    double rnew = m_fRad * angle / (0.5 * m_fAngle);
    
    x_out = xn * (rnew / r) + m_fCx;
    y_out = yn * (rnew / r) + m_fCy;
}

void ImageData::UndistortPoint(double x, double y, 
			       double &x_out, double &y_out) const
{
    if (!m_fisheye) {
	// UndistortPointRD(x, y, x_out, y_out);

	x_out = x;
	y_out = y;

	return;
    }

    double xn = x - m_fCx;
    double yn = y - m_fCy;
    
    double r = sqrt(xn * xn + yn * yn);
    double angle = 0.5 * m_fAngle * (r / m_fRad);
    double rnew = m_fFocal * tan(DEG2RAD(angle));
    
    x_out = xn * (rnew / r); // + m_fCx;
    y_out = yn * (rnew / r); // + m_fCy;
}

/* Create a pinhole view for the keys */
void ImageData::UndistortKeys() {
    if (!m_fisheye)
	return;

    int num_keys = (int) m_keys.size();
    if (num_keys == 0)
	return;
    
    for (int i = 0; i < num_keys; i++) {
	double x = m_keys[i].m_x;
	double y = m_keys[i].m_y;
	double x_new, y_new;

	UndistortPoint(x, y, x_new, y_new);
	
	m_keys[i].m_x = x_new;
	m_keys[i].m_y = y_new;
    }
}

std::vector<Keypoint> ImageData::UndistortKeysCopy() {
    if (!m_fisheye)
	return m_keys;

    int num_keys = (int) m_keys.size();
    if (num_keys == 0)
	return m_keys;
    
    std::vector<Keypoint> keys_new;
    keys_new.resize(num_keys);

    for (int i = 0; i < num_keys; i++) {
	double x = m_keys[i].m_x;
	double y = m_keys[i].m_y;
	double x_new, y_new;

	UndistortPoint(x, y, x_new, y_new);

	keys_new[i].m_x = x_new;
	keys_new[i].m_y = y_new;        
    }

    return keys_new;
}

/* Undistort the image to a different size than the original */
img_t *ImageData::UndistortImageResize(int w_new, int h_new)
{
    bool unload = false;
    
    /* Look for a cached image */
    char basename[256];
    strcpy(basename, m_name);
    basename[strlen(basename) - 4] = 0;

    char undistort_jpg_buf[256];
    sprintf(undistort_jpg_buf, "%s.undistort.jpg", basename);

    char undistort_bmp_buf[256];
    sprintf(undistort_bmp_buf, "%s.undistort.bmp", basename);

    if (FileExists(undistort_jpg_buf)) {
	return LoadJPEG(undistort_jpg_buf);
    }

    if (FileExists(undistort_bmp_buf)) {
	return img_read_bmp_file(undistort_bmp_buf);
    }

    if (!m_image_loaded) {
	unload = true;
	LoadImage();
    }
    
    if (!m_fisheye) {
	img_t *out = img_copy(m_img);  /* No need to undistort */

	if (unload)
	    UnloadImage();

	return out;
    }


    double Rsph[9];
    GetRotationFromSpherical(-0.5 * M_PI, 0.5 * M_PI, Rsph);

    /* Rotate 90 degrees around the z-axis */
#if 0
    double Rz[9] = { 0, -1, 0,
		     1,  0, 0,
		     0,  0, 1  };
#else
    double Rz[9] = { 1, 0, 0,
		     0, 1, 0,
		     0, 0, 1  };
#endif

    double R[9];
    matrix_product(3, 3, 3, 3, Rsph, Rz, R);

    int width = GetWidth();
    int height = GetHeight();

    img_t *img_out = img_new(w_new, h_new);

    for (int y = 0; y < h_new; y++) {
	for (int x = 0; x < w_new; x++) {
	    double xn = x - 0.5 * w_new;
	    double yn = y - 0.5 * h_new;
	    double x_new, y_new;
	    
	    DistortPoint(xn, yn, R, x_new, y_new);

	    x_new += 0.5 * width;
	    y_new += 0.5 * height;

	    if (x_new < 0 || x_new >= width || y_new < 0 || y_new >= height)
		continue;

	    fcolor_t c = pixel_lerp(m_img, x_new, y_new);

	    img_set_pixel(img_out, x, y, 
			  iround(c.r), iround(c.g), iround(c.b));
	}
    }

    img_write_bmp_file(img_out, undistort_bmp_buf);

    if (unload)
	UnloadImage();

    return img_out;
}


img_t *ImageData::UndistortImage(double theta, double phi, int num_rot) {
    bool unload = false;
    
    if (!m_image_loaded) {
	printf("[ImageData::UndistortImage] Loading image...\n");
	unload = true;
	LoadImage();
    }

#if 1   
    if (!m_fisheye) {
	img_t *out = img_copy(m_img);  /* No need to undistort */

	if (unload)
	    UnloadImage();

	return out;
    }
#endif

    double Rsph[9];
    GetRotationFromSpherical(theta - 0.5 * M_PI, phi + 0.5 * M_PI, Rsph);

    /* Rotate 90 degrees around the z-axis */
    double R90[9] = { 0, -1, 0,
		      1,  0, 0,
		      0,  0, 1  };

    double Rz[9];
    matrix_power(3, R90, num_rot, Rz);

    double R[9];
    matrix_product(3, 3, 3, 3, Rsph, Rz, R);

    int width = GetWidth();
    int height = GetHeight();

    img_t *img_out = 
	img_new(width / SUBSAMPLE_LEVEL, height / SUBSAMPLE_LEVEL);

    for (int y = 0; y < height / SUBSAMPLE_LEVEL; y++) {
	for (int x = 0; x < width / SUBSAMPLE_LEVEL; x++) {
	    double xn = x - 0.5 * width;
	    double yn = y - 0.5 * height;
	    double x_new, y_new;
	    
	    DistortPoint(xn, yn, R, x_new, y_new);

	    x_new += 0.5 * width;
	    y_new += 0.5 * height;

	    if (x_new < 0 || x_new >= width || y_new < 0 || y_new >= height)
		continue;

	    fcolor_t c = pixel_lerp(m_img, x_new, y_new);

	    img_set_pixel(img_out, x, y, 
			  iround(c.r), iround(c.g), iround(c.b));
	}
    }
    
    if (unload) {
	printf("[ImageData::UndistortImage] Unloading image\n");
	UnloadImage();
    }

    return img_out;
}

static void ProjectOnLine(double *p, double *line, double *proj) 
{
    double vec[2] = { -line[2] * line[0], -line[2] * line[1] };
    double p_diff[2] = { p[0] - vec[0], p[1] - vec[1] };
    double dot = p_diff[0] * line[0] + p_diff[1] * line[1];
    double par[2] = { dot * line[0], dot * line[1] };
    double perp[2] = { p_diff[0] - par[0], p_diff[1] - par[1] };
    proj[0] = perp[0] + vec[0];
    proj[1] = perp[1] + vec[1];
}

/* Find line segments in the image */
void ImageData::DetectLineSegments(double sigma, 
				   double threshold1, double threshold2,
				   double min_line_segment_size) 
{
    bool unload = false;

    if (!m_image_loaded) {
        printf("[ImageData::DetectLineSegments] Loading image...\n");
        LoadImage();
        unload = true;
    }

    int w = m_img->w;
    int h = m_img->h;

#if 0
    const int num_angles = 3;
    double angles[num_angles] = { 0.0, -0.125 * M_PI, 0.125 * M_PI };
#else
    const int num_angles = 1;
    double angles[num_angles] = { 0.0 /* -0.125 * M_PI */ };
#endif

    for (int angle_i = 0; angle_i < num_angles; angle_i++) {
        double angle = angles[angle_i];

        // img_out = img_canny_edge_detect(img, 0.8, 0.0350, 0.005);
#if 0
        /* Detect lines on the fisheye images */
        img_t *img_out = img_canny_edge_detect(m_img, sigma, 
            threshold1, threshold2);

        int num_links;
        edge_link_t *links;
        img_link_edges(img_out, &num_links, &links);

        /* Undistort the edges given the fisheye parameters */
        if (m_fisheye) {
            for (int i = 0; i < num_links; i++) {
                edge_link_t *link = links[i].next;

                while (link != NULL) {
                    double x_new, y_new;

                    UndistortPoint(link->x - 0.5 * w, link->y - 0.5 * h, 
                        x_new, y_new);

                    link->x = x_new;
                    link->y = y_new;

                    link = link->next;
                }
            }
        }
#else
        img_t *img_undistort = UndistortImage(0.0, angle);
        img_t *img_out = 
            img_canny_edge_detect(img_undistort, sigma, 
            threshold1, threshold2);

        int num_links;
        edge_link_t *links;
        img_link_edges(img_out, &num_links, &links);

        /* Make the image center the origin */
        for (int i = 0; i < num_links; i++) {
            edge_link_t *link = links[i].next;

            while (link != NULL) {
                double x_new, y_new;

                x_new = link->x - 0.5 * w;
                y_new = link->y - 0.5 * h;

                link->x = x_new;
                link->y = y_new;

                link = link->next;
            }
        }

#if 0
        char buf[256];
        sprintf(buf, "canny%03d.bmp", angle_i);
        img_write_bmp_file(img_out, buf);
#endif

        img_free(img_undistort);
#endif

        int num_lines;
        edge_link_t *lines;
        img_break_edges(num_links, &links, &num_lines, &lines, 2.5 /* 4.0 */ );

        int num_lines_long;
        edge_link_t *lines_long;
        // edge_remove_small(num_lines, &lines, &num_lines_long, &lines_long, 25);
        edge_remove_small(num_lines, &lines, 
            &num_lines_long, &lines_long, 
            iround(min_line_segment_size));

        /* Copy the links into our own data structure */
        for (int i = 0; i < num_lines_long; i++) {

            edge_link_t *link = lines_long[i].next;
            int length = 0;
            while (link != NULL) {
                link = link->next;
                length++;
            }

            v2_t *pts = new v2_t[length];	
            link = lines_long[i].next;
            int count = 0;
            while (link != NULL) {
                pts[count] = v2_new(link->x, link->y);
                link = link->next;
                count++;
            }

            double params[3];
            // double error = 
            fit_2D_line_orthogonal_regression(length, pts, params);

            // printf("error = %0.3f\n", error);

            /* Project the endpoints on the line */
            v2_t *projs = new v2_t[length];
            for (int j = 0; j < length; j++) {
                double pt_arr[2] = { Vx(pts[j]), Vy(pts[j]) };
                double proj_arr[2];
                ProjectOnLine(pt_arr, params, proj_arr);
                projs[j] = v2_new(proj_arr[0], proj_arr[1]);
            }

            /* Compute the mean of the projections */
            v2_t mean = v2_new(0.0, 0.0);
            for (int j = 0; j < length; j++) {
                mean = v2_add(mean, projs[j]);
            }
            mean = v2_scale(1.0 / length, mean);

            /* Find the extrema of the projected points */
            double norm_max = 0.0;
            int idx_max = -1;
            for (int j = 0; j < length; j++) {
                v2_t p = v2_sub(projs[j], mean);
                double norm = v2_norm(p);
                if (norm > norm_max) {
                    norm_max = norm;
                    idx_max = j;
                }
            }

            int extrema1 = idx_max;
            norm_max = 0.0;
            idx_max = -1;
            for (int j = 0; j < length; j++) {
                v2_t p = v2_sub(projs[j], projs[extrema1]);
                double norm = v2_norm(p);
                if (norm > norm_max) {
                    norm_max = norm;
                    idx_max = j;
                }
            }

            int extrema2 = idx_max;

            /* Project the endpoints on the main image plane */
            if (m_fisheye) {
                double R[9];
                GetRotationFromSpherical(-0.5 * M_PI, angle + 0.5 * M_PI, R);

                double p1[3] = { Vx(projs[extrema1]), 
                    Vy(projs[extrema1]), -m_fFocal };
                double p2[3] = { Vx(projs[extrema2]), 
                    Vy(projs[extrema2]), -m_fFocal };

                double Rp1[3], Rp2[3];
                matrix_product(3, 3, 3, 1, R, p1, Rp1);
                matrix_product(3, 3, 3, 1, R, p2, Rp2);

                matrix_scale(3, 1, Rp1, m_fFocal / Rp1[2], Rp1);
                matrix_scale(3, 1, Rp2, m_fFocal / Rp2[2], Rp2);

                LineSegment2D line;
                line.m_p1[0] = Rp1[0]; // Vx(projs[extrema1]);
                line.m_p1[1] = Rp1[1]; // Vy(projs[extrema1]);

                line.m_p2[0] = Rp2[0]; // Vx(projs[extrema2]);
                line.m_p2[1] = Rp2[1]; // Vy(projs[extrema2]);

                m_line_segments.push_back(line);	
            } else {
                LineSegment2D line;
                line.m_p1[0] = Vx(projs[extrema1]);
                line.m_p1[1] = Vy(projs[extrema1]);

                line.m_p2[0] = Vx(projs[extrema2]);
                line.m_p2[1] = Vy(projs[extrema2]);                

                m_line_segments.push_back(line);	
            }

            delete [] pts;
            delete [] projs;
        }

        /* Cleanup */
        printf("[ImageData::DetectLineSegments] Detected %d lines\n", 
            num_lines_long);
        printf("[ImageData::DetectLineSegments] Freeing %d lines\n", 
            num_lines);

        for (int i = 0; i < num_lines; i++) {
            edge_free(&lines[i]);
        }

        free(links);
        free(lines);
        free(lines_long);

        img_free(img_out);
    }

    if (unload) {
        printf("[ImageData::DetectLineSegments] Unloading image\n");
        UnloadImage();
    }
}

/* Render line segments to an image */
img_t *ImageData::RenderLineSegments() 
{
    img_t *segment_img = UndistortImage();

    int w = GetWidth();
    int h = GetHeight();

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            color_t p = img_get_pixel(segment_img, x, y);
            double r = 0.75 * p.r + 0.25 * 0xff;
            double g = 0.75 * p.g + 0.25 * 0xff;
            double b = 0.75 * p.b + 0.25 * 0xff;

            int r_int = CLAMP(iround(r), 0, 255);
            int g_int = CLAMP(iround(g), 0, 255);
            int b_int = CLAMP(iround(b), 0, 255);

            img_set_pixel(segment_img, x, y, r_int, g_int, b_int);
        }
    }

    int num_segments = (int) m_line_segments.size();

    for (int i = 0; i < num_segments; i++) {
        int x1 = iround(m_line_segments[i].m_p1[0] + 0.5 * w);
        int y1 = iround(m_line_segments[i].m_p1[1] + 0.5 * h);
        int x2 = iround(m_line_segments[i].m_p2[0] + 0.5 * w);
        int y2 = iround(m_line_segments[i].m_p2[1] + 0.5 * h);

        int r = rand() % 256;
        int g = rand() % 256;
        int b = rand() % 256;

        img_draw_line(segment_img, x1, y1, x2, y2, r, g, b);
    }

    return segment_img;
}

/* Query whether a given line segment is supported by the image */
bool ImageData::LineSegmentSupportedByImage(LineSegment2D &line) 
{
    int w = GetWidth();
    int h = GetHeight();

    /* Check that the line is completely inside the image */
    if (line.m_p1[0] + 0.5 * w < 5.0 || line.m_p1[0] + 0.5 * w > w - 6 ||
	line.m_p1[1] + 0.5 * h < 5.0 || line.m_p1[1] + 0.5 * h > h - 6)
	return false;

    if (line.m_p2[0] + 0.5 * w < 5.0 || line.m_p2[0] + 0.5 * w > w - 6 ||
	line.m_p2[1] + 0.5 * h < 5.0 || line.m_p2[1] + 0.5 * h > h - 6)
	return false;

#define NUM_SEGMENT_SAMPLES 12
    int num_supporting = 0;

    /* Find the homogeneous line */
    double l[3];
    line.Homogeneous(l);
    
    double mag = sqrt(l[0] * l[0] + l[1] * l[1]);
    matrix_scale(3, 1, l, 1.0 / mag, l);

    for (int i = 0; i < NUM_SEGMENT_SAMPLES; i++) {
	/* Compute the sample */
	double t = (double) i / (double) (NUM_SEGMENT_SAMPLES - 1);
	double sample[2];
	line.Sample(t, sample);

	/* Check the angle of the local gradient with the line segment */
	double grad[2];
	Gradient(sample[0], sample[1], grad);

	double mag = sqrt(grad[0] * grad[0] + grad[1] * grad[1]);
#define SUPPORT_MAGNITUDE_THRESHOLD 20.0
	if (mag > SUPPORT_MAGNITUDE_THRESHOLD) {
	    grad[0] /= mag;
	    grad[1] /= mag;

	    double dot;
	    matrix_product(1, 2, 2, 1, grad, l, &dot);

	    double angle = acos(fabs(dot));

#define SUPPORT_ANGLE_THRESHOLD 0.0872664625997165 /* 5 degrees */
	    if (angle < SUPPORT_ANGLE_THRESHOLD)
		num_supporting++;
	}
    }
    
#define SUPPORTING_PERCENT_THRESHOLD 0.5
    if (num_supporting >= SUPPORTING_PERCENT_THRESHOLD * NUM_SEGMENT_SAMPLES)
	return true;
    
    return false;
}

/* Return true if the given line intersects the image */
bool ImageData::LineIntersectsImage(double *line, 
				    double *isect1, double *isect2,
				    ImageBoundary &b1, ImageBoundary &b2)
{
    double  upper_left[3] = { -0.5 * m_width,  0.5 * m_height, 1.0 };
    double  lower_left[3] = { -0.5 * m_width, -0.5 * m_height, 1.0 };
    double upper_right[3] = {  0.5 * m_width,  0.5 * m_height, 1.0 };
    double lower_right[3] = {  0.5 * m_width, -0.5 * m_height, 1.0 };

    double north[3], east[3], south[3], west[3];
    matrix_cross(upper_right, upper_left, north);
    matrix_cross(upper_right, lower_right, east);
    matrix_cross(lower_left, lower_right, south);
    matrix_cross(lower_left, upper_left, west);

    double north_isect[3], east_isect[3], south_isect[3], west_isect[3];
    matrix_cross(north, line, north_isect);
    matrix_cross(east, line,  east_isect);
    matrix_cross(south, line, south_isect);
    matrix_cross(west, line,  west_isect);

    matrix_scale(3, 1, north_isect, 1.0 / north_isect[2], north_isect);
    matrix_scale(3, 1, east_isect, 1.0 / east_isect[2], east_isect);
    matrix_scale(3, 1, south_isect, 1.0 / south_isect[2], south_isect);
    matrix_scale(3, 1, west_isect, 1.0 / west_isect[2], west_isect);

    int num_isect = 0;
    if (north_isect[0] >= -0.5 * m_width && north_isect[0] <= 0.5 * m_width) {
	b1 = BoundaryNorth;
	memcpy(isect1, north_isect, 3 * sizeof(double));
	num_isect++;
    }

    if (south_isect[0] >= -0.5 * m_width && south_isect[0] <= 0.5 * m_width) {
	if (num_isect == 0) {
	    b1 = BoundarySouth;
	    memcpy(isect1, south_isect, 3 * sizeof(double));	    
	} else if (num_isect == 1) {
	    b2 = BoundarySouth;
	    memcpy(isect2, south_isect, 3 * sizeof(double));	    
	}
	
	num_isect++;
    }

    if (east_isect[1] >= -0.5 * m_height && east_isect[1] <= 0.5 * m_height) {
	if (num_isect == 0) {
	    b1 = BoundaryEast;
	    memcpy(isect1, east_isect, 3 * sizeof(double));	    
	} else if (num_isect == 1) {
	    b2 = BoundaryEast;
	    memcpy(isect2, east_isect, 3 * sizeof(double));	    
	} else {
	    printf("[ImageData::LineIntersectsImage] Strange number of "
		   "intersections\n");
	}

	num_isect++;
    }

    if (west_isect[1] >= -0.5 * m_height && west_isect[1] <= 0.5 * m_height) {
	if (num_isect == 0) {
	    b1 = BoundaryWest;
	    memcpy(isect1, west_isect, 3 * sizeof(double));	    
	} else if (num_isect == 1) {
	    b2 = BoundaryWest;
	    memcpy(isect2, west_isect, 3 * sizeof(double));	    
	} else {
	    printf("[ImageData::LineIntersectsImage] Strange number of "
		   "intersections\n");
	}

	num_isect++;
    }
    
    if (num_isect == 2)
	return true;
    else if (num_isect == 0)
	return false;
    else {
	printf("[ImageData::LineIntersectsImage] num_isect = %d "
	       "doesn't make sense\n", num_isect);
	return false;
    }
}

/* Sobel kernels */
static const double sobel_x[9] = { -1, 0, 1,
				   -2, 0, 2,
				   -1, 0, 1 };

static const double sobel_y[9] = { -1, -2, -1, 
				    0,  0,  0,
				    1,  2,  1 };

/* Find the gradient direction at the given position */
void ImageData::Gradient(double x, double y, double *grad) 
{
    double gx = 0.0, gy = 0.0;

    int w = GetWidth();
    int h = GetHeight();

    double ident[9];
    GetRotationFromSpherical(-0.5 * M_PI, 0.5 * M_PI, ident);

    int idx = 0;
    for (int dy = -1; dy <= 1; dy++) {
	for (int dx = -1; dx <= 1; dx++, idx++) {
	    double x_d, y_d;
	    
	    /* Distort */
	    DistortPoint(x + dx, y + dy, ident, x_d, y_d);

	    fcolor_t c = pixel_lerp(m_img, x_d + 0.5 * w, y_d + 0.5 * h);
	    double v = fcolor_intensity(c);
	    
	    gx += sobel_x[idx] * v;
	    gy += sobel_y[idx] * v;
	}
    }

    // double mag = sqrt(gx * gx + gy * gy);
    // grad[0] = gx / mag;
    // grad[1] = gy / mag;

    grad[0] = gx;
    grad[1] = gy;
}

void ImageData::ReadKeyColors()
{
    double ident[9];
    GetRotationFromSpherical(-0.5 * M_PI, 0.5 * M_PI, ident);

    bool unload = false;
    if (!m_image_loaded) {
	LoadImage();
	unload = true;
    }
    
    int w = GetWidth();
    int h = GetHeight();

    int num_keys = (int) m_keys.size();
	
    for (int i = 0; i < num_keys; i++) {
	double x = m_keys[i].m_x;
	double y = m_keys[i].m_y;

	double x_d, y_d;
	DistortPoint(x, y, ident, x_d, y_d);

	fcolor_t col;
	col = pixel_lerp(m_img, x_d + 0.5 * w, y_d + 0.5 * h);

	m_keys[i].m_r = iround(col.r);
	m_keys[i].m_g = iround(col.g);
	m_keys[i].m_b = iround(col.b);
    }
    
    if (unload)
	UnloadImage();
}

/* Read the camera */
bool ImageData::ReadCamera()
{
    char cam_buf[256];
    strcpy(cam_buf, m_name);
    cam_buf[strlen(m_name) - 3] = 'c';
    cam_buf[strlen(m_name) - 2] = 'a';
    cam_buf[strlen(m_name) - 1] = 'm';

    FILE *f = fopen(cam_buf, "r");
    
    if (f == NULL) {
	return false;
    }

    fscanf(f, "%lf %lf %lf\n", 
           &(m_camera.m_focal), m_camera.m_k+0, m_camera.m_k+1);
    fscanf(f, "%lf %lf %lf\n", 
	   m_camera.m_R + 0, m_camera.m_R + 1, m_camera.m_R + 2);
    fscanf(f, "%lf %lf %lf\n", 
	   m_camera.m_R + 3, m_camera.m_R + 4, m_camera.m_R + 5);
    fscanf(f, "%lf %lf %lf\n", 
	   m_camera.m_R + 6, m_camera.m_R + 7, m_camera.m_R + 8);
    fscanf(f, "%lf %lf %lf\n", 
	   m_camera.m_t + 0, m_camera.m_t + 1, m_camera.m_t + 2);

    fclose(f);

    m_camera.m_width = GetWidth();
    m_camera.m_height = GetHeight();

#ifndef __BUNDLER__
    for (int i = 0; i < NUM_LINK_DIRECTIONS; i++) {
	m_camera.m_links[i] = -1;
    }
#endif

    m_camera.m_adjusted = true;
    m_camera.Finalize();

    return true;
}

/* Read the tracks */
bool ImageData::ReadTracks(int img_idx, std::vector<PointData> &pt_data)
{
    char cam_buf[256];
    strcpy(cam_buf, m_name);
    cam_buf[strlen(m_name) - 3] = 't';
    cam_buf[strlen(m_name) - 2] = 'r';
    cam_buf[strlen(m_name) - 1] = 'k';

    FILE *f = fopen(cam_buf, "r");
    
    if (f == NULL) {
	printf("[ImageData::ReadTracks] Error opening file %s for reading\n",
	       cam_buf);
	return false;
    }

    m_visible_points.clear();

    int num_tracks;
    fscanf(f, "%d\n", &num_tracks);
    
    /* Read the tracks */
    for (int i = 0; i < num_tracks; i++) {
	int key, track;
	fscanf(f, "%d %d\n", &key, &track);
	
#if 0
	if (key < 0 || key >= num_keys) 
	    printf("[ImageData::ReadTracks] Key %d out of bounds\n", key);
	else
	    m_keys[key] = track;
#endif
	
	m_visible_points.push_back(track);
        m_visible_keys.push_back(key);
        pt_data[track].m_views.push_back(ImageKey(img_idx, key));
    }
    
    fclose(f);

    return true;
}

/* Write the camera */
void ImageData::WriteCamera()
{
    char cam_buf[256];
    strcpy(cam_buf, m_name);
    cam_buf[strlen(m_name) - 3] = 'c';
    cam_buf[strlen(m_name) - 2] = 'a';
    cam_buf[strlen(m_name) - 1] = 'm';

    FILE *f = fopen(cam_buf, "w");
    
    if (f == NULL) {
	printf("[ImageData::WriteCamera] Error opening file %s for writing\n",
	       cam_buf);
	return;
    }
    
    fprintf(f, "%0.8e %0.8e %0.8e\n", 
            m_camera.m_focal, m_camera.m_k[0], m_camera.m_k[1]);
    fprintf(f, "%0.8e %0.8e %0.8e\n", 
	    m_camera.m_R[0], m_camera.m_R[1], m_camera.m_R[2]);
    fprintf(f, "%0.8e %0.8e %0.8e\n", 
	    m_camera.m_R[3], m_camera.m_R[4], m_camera.m_R[5]);
    fprintf(f, "%0.8e %0.8e %0.8e\n", 
	    m_camera.m_R[6], m_camera.m_R[7], m_camera.m_R[8]);
    fprintf(f, "%0.8e %0.8e %0.8e\n", 
	    m_camera.m_t[0], m_camera.m_t[1], m_camera.m_t[2]);

    fclose(f);
}

/* Write the camera in XML format */
void ImageData::WriteCameraXML(FILE *f)
{
    static char *spacer = "  ";

    fprintf(f, "%s<camera>\n", spacer);
    fprintf(f, "%s  <w> %d </w>\n%s  <h> %d </h>\n", 
	    spacer, GetWidth(), spacer, GetHeight());
    fprintf(f, "%s  <adj> %d </adj>\n", spacer, m_camera.m_adjusted ? 1 : 0);

    char jpeg_buf[256];
    strcpy(jpeg_buf, m_name);
    jpeg_buf[strlen(m_name) - 3] = 'j';
    jpeg_buf[strlen(m_name) - 2] = 'p';
    jpeg_buf[strlen(m_name) - 1] = 'g';

    char flickr_url[1024];
    sprintf(flickr_url, "http://www.flickr.com/photos/%s/%s/",
            m_user_name, m_flickr_index);

    fprintf(f, "%s  <name> %s </name>\n", spacer, jpeg_buf);
    fprintf(f, "%s  <flickr> %s </flickr>\n", spacer, flickr_url);
    fprintf(f, "%s  <title> %s </title>\n", spacer, "");
    fprintf(f, "%s  <takenBy> %s </takenBy>\n", spacer, "");
    fprintf(f, "%s  <takenOn> %s </takenOn>\n", spacer, "");
    fprintf(f, "%s  <perms> %d </perms>\n", spacer, m_licensed ? 1 : 0);

    if (!m_camera.m_adjusted) {
	fprintf(f, "%s</camera>\n", spacer);
	return;
    }
    
    m_camera.WriteXML(f);

    /* Compute the 3D points on the projection plane */
    double eye[3];
    m_camera.GetPosition(eye);
    int w = GetWidth();
    int h = GetHeight();
    
    double ray1[3] = { -0.5 * w, -0.5 * h, -m_camera.m_focal };
    double ray2[3] = {  0.5 * w, -0.5 * h, -m_camera.m_focal };
    double ray3[3] = { -0.5 * w,  0.5 * h, -m_camera.m_focal };
    double ray4[3] = {  0.5 * w,  0.5 * h, -m_camera.m_focal };

    double ray_world[18];
    double R[9];
    m_camera.GetPose(R);
    matrix_product(3, 3, 3, 1, R, ray1, ray_world + 0);
    matrix_product(3, 3, 3, 1, R, ray2, ray_world + 3);
    matrix_product(3, 3, 3, 1, R, ray3, ray_world + 6);
    matrix_product(3, 3, 3, 1, R, ray4, ray_world + 9);
    
    double isect[18];
    double t0 = m_fit_plane.IntersectRay(eye, ray_world + 0, isect + 0);
    double t1 = m_fit_plane.IntersectRay(eye, ray_world + 3, isect + 3);
    double t2 = m_fit_plane.IntersectRay(eye, ray_world + 6, isect + 6);
    double t3 = m_fit_plane.IntersectRay(eye, ray_world + 9, isect + 9);

    if (t0 < 0.0 || t1 < 0.0 || t2 < 0.0 || t3 < 0.0) {
        fprintf(f, "%s  <p1> 0.0 0.0 0.0 </p1>\n", spacer);
        fprintf(f, "%s  <p2> 0.0 0.0 0.0 </p2>\n", spacer);
        fprintf(f, "%s  <p3> 0.0 0.0 0.0 </p3>\n", spacer);
        fprintf(f, "%s  <p4> 0.0 0.0 0.0 </p4>\n", spacer);
    } else {
        fprintf(f, "%s  <p1> %0.6e %0.6e %0.6e </p1>\n", 
                spacer, isect[0], isect[1], isect[2]);
        fprintf(f, "%s  <p2> %0.6e %0.6e %0.6e </p2>\n", 
                spacer, isect[3], isect[4], isect[5]);
        fprintf(f, "%s  <p3> %0.6e %0.6e %0.6e </p3>\n", 
                spacer, isect[6], isect[7], isect[8]);
        fprintf(f, "%s  <p4> %0.6e %0.6e %0.6e </p4>\n", 
                spacer, isect[9], isect[10], isect[11]);        
    }

    fprintf(f, "%s</camera>\n", spacer);
}




/* Write the tracks */
void ImageData::WriteTracks()
{
    char cam_buf[256];
    strcpy(cam_buf, m_name);
    cam_buf[strlen(m_name) - 3] = 't';
    cam_buf[strlen(m_name) - 2] = 'r';
    cam_buf[strlen(m_name) - 1] = 'k';

    FILE *f = fopen(cam_buf, "w");
    
    if (f == NULL) {
	printf("[ImageData::WriteTracks] Error opening file %s for writing\n",
	       cam_buf);
	return;
    }

    /* Count the number of tracks */
    int num_keys = (int) m_keys.size();
    int num_tracks = 0;
    
    for (int i = 0; i < num_keys; i++)
	if (m_keys[i].m_extra != -1)
	    num_tracks++;

    fprintf(f, "%d\n", num_tracks);
    
    /* Write the tracks */
    for (int i = 0; i < num_keys; i++)
	if (m_keys[i].m_extra != -1)
	    fprintf(f, "%d %d\n", i, m_keys[i].m_extra);

    fclose(f);
}

void ImageData::ReadMetadata()
{
    char meta_buf[256];
    strcpy(meta_buf, m_name);
    meta_buf[strlen(m_name) - 3] = 'm';
    meta_buf[strlen(m_name) - 2] = 'e';
    meta_buf[strlen(m_name) - 1] = 't';
    
    FILE *f = fopen(meta_buf, "r");
    
    if (f == NULL) return;
    
    char buf[256];
    fgets(buf, 256, f);

    if (buf[strlen(buf) - 1] == '\n')
	buf[strlen(buf) - 1] = 0;
    
    if (buf[strlen(buf) - 1] == '\r')
	buf[strlen(buf) - 1] = 0;

    m_real_name = strdup(buf);

    fgets(buf, 256, f);

    if (buf[strlen(buf) - 1] == '\n')
	buf[strlen(buf) - 1] = 0;
    
    if (buf[strlen(buf) - 1] == '\r')
	buf[strlen(buf) - 1] = 0;

    strcpy(m_user_name, buf);
    
    fscanf(f, "%d:%d:%d %d:%d:%d\n", 
	   &(m_date.m_year), &(m_date.m_month), &(m_date.m_day),
	   &(m_date.m_hour), &(m_date.m_minute), &(m_date.m_second));

    fclose(f);
}

/* Get notes associated with this image */
std::vector<ImageNote> ImageData::GetNotes()
{
    char note_buf[256];
    strcpy(note_buf, m_name);
    note_buf[strlen(m_name) - 3] = 'n';
    note_buf[strlen(m_name) - 2] = 'o';
    note_buf[strlen(m_name) - 1] = 't';
    
    FILE *f = fopen(note_buf, "r");
    
    std::vector<ImageNote> notes;
    if (f == NULL)
	return notes;

    int num_notes;
    fscanf(f, "%d\n", &num_notes);

    for (int i = 0; i < num_notes; i++) {
	ImageNote note;
	note.Read(f, (double) GetWidth(), (double) GetHeight());
	notes.push_back(note);
    }

    return notes;
}

#ifndef __BUNDLER__
#define GRID_SIZE 25
/* Compute 3D points used for estimating the amount of distortion
 * in the image */
void ImageData::ComputeDistortionPoints(const std::vector<PointData> &pt_data,
                                        std::vector<DPoint3> &pts_world,
                                        std::vector<DPoint3> &pts_plane)
{
    std::vector<int> *grid = new std::vector<int> [GRID_SIZE * GRID_SIZE];

    int num_vis_points = (int) m_visible_points.size();
    double w_inv = 1.0 / GetWidth();
    double h_inv = 1.0 / GetHeight();

    for (int i = 0; i < num_vis_points; i++) {
        int idx = m_visible_points[i];
        double proj[2];

        m_camera.Project(pt_data[idx].m_pos, proj);

        proj[0] *= w_inv + 0.5;
        proj[1] *= h_inv + 0.5;

        proj[0] = CLAMP(proj[0], 0.0, 0.9999);
        proj[1] = CLAMP(proj[1], 0.0, 0.9999);

        int g_x = iround(floor(GRID_SIZE * proj[0]));
        int g_y = iround(floor(GRID_SIZE * proj[1]));
        
        assert(g_x >= 0 && g_x < GRID_SIZE && g_y >= 0 && g_y < GRID_SIZE);

        grid[g_y * GRID_SIZE + g_x].push_back(idx);
    }

    /* Choose one point from each cell */
    double eye[3];
    m_camera.GetPosition(eye);

    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        int num_bin_points = (int) grid[i].size();
        
        if (num_bin_points == 0)
            continue;
        
        int bin_idx = rand() % num_bin_points;
        int pt_idx = grid[i][bin_idx];
        
        pts_world.push_back(DPoint3(pt_data[pt_idx].m_pos));

        /* Project the point into the image, and reproject onto the
         * plane */
        double proj[2];
        m_camera.Project(pt_data[pt_idx].m_pos, proj);
        
        double ray[3];
        m_camera.PixelToCameraRayAbsolute(proj[0], proj[1], ray);

        double isect[3];
        m_fit_plane.IntersectRay(eye, ray, isect);

        pts_plane.push_back(DPoint3(isect));
    }

    delete [] grid;

    m_dist_pts_world = pts_world;
    m_dist_pts_plane = pts_plane;
}
#undef GRID_SIZE

#if 0
void ImageData::ComputeOrthoPlane(const std::vector<PointData> &pt_data)
{
#if 0
    int num_vis_points = (int) m_visible_points.size();
    
    double mean[3] = { 0.0, 0.0, 0.0 };

    for (int i = 0; i < num_vis_points; i++) {
        int idx = m_visible_points[i];
        double *pos = (double *) pt_data[idx].m_pos;
        matrix_sum(3, 1, 3, 1, mean, pos, mean);
    }
    
    matrix_scale(3, 1, mean, 1.0 / num_vis_points, mean);
#endif

    double mean[3] = { -0.043, 0.120, -0.500 };
    double view[3];
    m_camera.GetViewDirection(view);

#if 1
    /* Set y-axis to 0.0 */
    view[1] = 0.0;
    double norm = matrix_norm(3, 1, view);
    matrix_scale(3, 1, view, 1.0 / norm, view);
#endif

    double dot;
    matrix_product(1, 3, 3, 1, mean, view, &dot);

    m_ortho_plane.m_normal[0] = view[0];
    m_ortho_plane.m_normal[1] = view[1];
    m_ortho_plane.m_normal[2] = view[2];
    m_ortho_plane.m_dist = -dot;
}
#endif
#endif /* __BUNDLER__ */

#define MAX_DIM 256
void ImageData::ComputeFeatureWeightMap(int id, 
                                        const std::vector<PointData> &pt_data)
{
    char buf[1024];
    // GetBaseName(buf);
    strcpy(buf, m_name);

    int len = strlen(buf);
    buf[len-3] = 'w';
    buf[len-2] = 't';
    buf[len-1] = '.';
    buf[len+0] = 'b';
    buf[len+1] = 'm';
    buf[len+2] = 'p';
    buf[len+3] = 0;

    if (FileExists(buf)) {
        m_feature_weight_map = img_read_bmp_file(buf);
        return;
    }

    LoadKeysWithScaleRot(false);
    LoadImage();

    double scale;
    // img_t *scaled = RescaleImage(m_img, MAX_DIM, scale);
    int dim = MAX(m_img->w, m_img->h);
    double ratio = (double) dim / (double) MAX_DIM;
    scale = 1.0 / ratio;

    // int w = scaled->w, h = scaled->h;
    int w = iround(m_img->w * scale + 1);
    int h = iround(m_img->h * scale + 1);
    double *weights = new double[w * h];

    for (int i = 0; i < w * h; i++)
        weights[i] = 0.0;

    m_feature_weight_map = img_new(w, h);

    int num_vis_pts = (int) m_visible_points.size();
    for (int i = 0; i < num_vis_pts; i++) {
        int pt_idx = m_visible_points[i];
        int key = -1;

        int num_views = (int) pt_data[pt_idx].m_views.size();
        for (int j = 0; j < num_views; j++) {
            if (pt_data[pt_idx].m_views[j].first == id) {
                key = pt_data[pt_idx].m_views[j].second;
                break;
            }
        }

        if (key == -1) {
            printf("[ImageData::ComputeFeatureWeightMap] "
                   "Couldn't find key for image %d!\n", id);
            continue;
        }

        double fscale = /*4.0*/ 3.0 * scale * m_keys_scale_rot[key].m_scale;
        // double x = scale * (m_keys_scale_rot[key].m_x + 0.5 * m_img->w);
        // double y = scale * (m_keys_scale_rot[key].m_y + 0.5 * m_img->h);

        double proj[2];
        m_camera.Project(pt_data[pt_idx].m_pos, proj);
        double x = scale * (proj[0] + 0.5 * m_img->w);
        double y = scale * (proj[1] + 0.5 * m_img->h);

        fscale = CLAMP(fscale, 1.0, 10.0);
        
        int width;
        double *kernel = compute_gaussian_filter(fscale, 3.0, &width);
        double mid = kernel[width / 2];
        matrix_scale(width, 1, kernel, 0.5 / mid, kernel);

        for (int dy = -width / 2; dy <= width / 2; dy++) {
            for (int dx = -width / 2; dx <= width / 2; dx++) {
                double wt = kernel[dy + width / 2] * kernel[dx + width / 2];
                int x_img = iround(x + dx);
                int y_img = iround(y + dy);

                if (x_img < 0 || x_img >= w || y_img < 0 || y_img >= h)
                    continue;
                
                weights[y_img * w + x_img] += wt;
            }
        }

        delete [] kernel;
    }

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int idx = y * w + x;
            double wt = CLAMP(weights[idx], 0.0, 1.0);
            int c = iround(wt * 255.0);
            img_set_pixel(m_feature_weight_map, x, y, c, c, c);
        }
    }
    
    img_write_bmp_file(m_feature_weight_map, buf);

    delete [] weights;
    // img_free(scaled);

    UnloadKeysWithScaleRot();
    UnloadImage();
}

bool ImageData::LoadFeatureWeightMap()
{
    char buf[1024];
    // GetBaseName(buf);
    strcpy(buf, m_name);

    int len = strlen(buf);
    buf[len-3] = 'w';
    buf[len-2] = 't';
    buf[len-1] = '.';
    buf[len+0] = 'b';
    buf[len+1] = 'm';
    buf[len+2] = 'p';
    buf[len+3] = 0;

    if (FileExists(buf)) {
        m_feature_weight_map = img_read_bmp_file(buf);
        return true;
    }

    return false;
}

void ImageData::SaveFeatureWeightMap()
{
    char buf[1024];
    // GetBaseName(buf);
    strcpy(buf, m_name);

    int len = strlen(buf);
    buf[len-3] = 'w';
    buf[len-2] = 't';
    buf[len-1] = '.';
    buf[len+0] = 'b';
    buf[len+1] = 'm';
    buf[len+2] = 'p';
    buf[len+3] = 0;

    img_write_bmp_file(m_feature_weight_map, buf);
}

void ImageData::ComputeFeatureWeightMap(int id, 
                                        const std::vector<PointData> &pt_data,
                                        const std::vector<int> &used_points)
{
    if (LoadFeatureWeightMap())
        return;

    LoadKeysWithScaleRot(false);
    LoadImage();

    double scale;
    // img_t *scaled = RescaleImage(m_img, MAX_DIM, scale);
    int dim = MAX(m_img->w, m_img->h);
    double ratio = (double) dim / (double) MAX_DIM;
    scale = 1.0 / ratio;

    // int w = scaled->w, h = scaled->h;
    int w = iround(m_img->w * scale + 1);
    int h = iround(m_img->h * scale + 1);
    double *weights = new double[w * h];

    for (int i = 0; i < w * h; i++)
        weights[i] = 0.0;

    m_feature_weight_map = img_new(w, h);

    std::vector<int> visible_points = 
        GetVectorIntersection(used_points, m_visible_points);
    
    int num_vis_pts = (int) visible_points.size();

    for (int i = 0; i < num_vis_pts; i++) {
        int pt_idx = visible_points[i];
        int key = -1;

        int num_views = (int) pt_data[pt_idx].m_views.size();
        for (int j = 0; j < num_views; j++) {
            if (pt_data[pt_idx].m_views[j].first == id) {
                key = pt_data[pt_idx].m_views[j].second;
                break;
            }
        }

        if (key == -1) {
            printf("[ImageData::ComputeFeatureWeightMap] "
                   "Couldn't find key for image %d!\n", id);
            continue;
        }

        double fscale = /*4.0*/ 3.0 * scale * m_keys_scale_rot[key].m_scale;
        double x = scale * (m_keys_scale_rot[key].m_x + 0.5 * m_img->w);
        double y = scale * (m_keys_scale_rot[key].m_y + 0.5 * m_img->h);

        fscale = CLAMP(fscale, 1.0, 10.0);
        
        int width;
        double *kernel = compute_gaussian_filter(fscale, 3.0, &width);
        double mid = kernel[width / 2];
        matrix_scale(width, 1, kernel, 0.5 / mid, kernel);

        for (int dy = -width / 2; dy <= width / 2; dy++) {
            for (int dx = -width / 2; dx <= width / 2; dx++) {
                double wt = kernel[dy + width / 2] * kernel[dx + width / 2];
                int x_img = iround(x + dx);
                int y_img = iround(y + dy);

                if (x_img < 0 || x_img >= w || y_img < 0 || y_img >= h)
                    continue;
                
                weights[y_img * w + x_img] += wt;
            }
        }

        delete [] kernel;
    }

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int idx = y * w + x;
            double wt = CLAMP(weights[idx], 0.0, 1.0);
            int c = iround(wt * 255.0);
            img_set_pixel(m_feature_weight_map, x, y, c, c, c);
        }
    }
    
    SaveFeatureWeightMap();

    delete [] weights;
    // img_free(scaled);

    UnloadKeysWithScaleRot();
    UnloadImage();    
}


#ifndef __BUNDLER__
void ImageData::ComputeTPSBasis(const std::vector<PointData> &pt_data)
{
    printf("[ImageData::ComputeTPSBasis] Computing basis for image %s\n",
           m_name);

    clock_t start = clock();

    int num_vis_points = (int) m_visible_points.size();

    std::vector<DPoint> pts;
    pts.resize(num_vis_points);

    for (int i = 0; i < num_vis_points; i++) {
        /* Project points into both views */
        int pt_idx = m_visible_points[i];
        double p[2];

        m_camera.Project(pt_data[pt_idx].m_pos, p);

        pts[i] = DPoint(p[0], p[1]);
    }

    m_tps_basis.m_basis_points = 
        GetThinPlateSplineBasis(pts, 
                                m_tps_basis.m_indices,
                                0.005 * MAX(GetWidth(), GetHeight()),
                                &(m_tps_basis.m_LU), &(m_tps_basis.m_ipiv));

    clock_t end = clock();
    
    printf("[ImageData::ComputeTPSBasis] Computed basis in %0.3fs\n",
           (double) (end - start) / CLOCKS_PER_SEC);
}

void ImageData::FreeTPSBasis()
{
    m_tps_basis.Clear();
}
#endif /* __BUNDLER__ */

void ImageData::ComputeHistogram(int nbins, double *r_hist, 
                                 double *g_hist, double *b_hist)
{
    LoadThumb256();
    
    if (!LoadFeatureWeightMap()) {
        printf("[ImageData::ComputeHistogram] Error: weight map for image %s "
               "not found!\n", m_name);

        return;
    }

    int w = m_thumb256->w;
    int h = m_thumb256->h;

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            color_t c = img_get_pixel(m_thumb256, x, y);
            color_t cwt = img_get_pixel(m_feature_weight_map, x, y);
            
            double wt = cwt.r / 255.0;

            double r_norm = c.r / 256.0 * nbins;
            double g_norm = c.g / 256.0 * nbins;
            double b_norm = c.b / 256.0 * nbins;

            int r_f = iround(floor(r_norm));
            int g_f = iround(floor(g_norm));
            int b_f = iround(floor(b_norm));

            int r_c = r_f + 1;
            int g_c = g_f + 1;
            int b_c = b_f + 1;

            assert(r_f <= nbins - 1);
            assert(g_f <= nbins - 1);
            assert(b_f <= nbins - 1);

            double t = 0.0;

            t = r_norm - r_f;
            if (r_f == nbins - 1) {
                r_hist[r_f] += wt;
            } else {
                r_hist[r_f] += wt * (1.0 - t);
                r_hist[r_c] += wt * t;
            }

            t = g_norm - g_f;
            if (g_f == nbins - 1) {
                g_hist[g_f] += wt;
            } else {
                g_hist[g_f] += wt * (1.0 - t);
                g_hist[g_c] += wt * t;
            }

            t = b_norm - b_f;
            if (b_f == nbins - 1) {
                b_hist[b_f] += wt;
            } else {
                b_hist[b_f] += wt * (1.0 - t);
                b_hist[b_c] += wt * t;
            }
        }
    }

    /* Normalize bins */
    double r_sum = 0.0, g_sum = 0.0, b_sum = 0.0;
    for (int i = 0; i < nbins; i++) {
        r_sum += r_hist[i];
        g_sum += g_hist[i];
        b_sum += b_hist[i];        
    }

    for (int i = 0; i < nbins; i++) {
        r_hist[i] /= r_sum;
        g_hist[i] /= g_sum;
        b_hist[i] /= b_sum;
    }    
}


void ImageData::ComputeHistogramLUV(int nbins, double *L_hist, 
                                    double *U_hist, double *V_hist)
{
    LoadThumb256();
    
    if (!LoadFeatureWeightMap()) {
        printf("[ImageData::ComputeHistogram] Error: weight map for image %s "
               "not found!\n", m_name);

        return;
    }

    int w = m_thumb256->w;
    int h = m_thumb256->h;

    const double Lmin = 0.0;
    const double Lmax = 100.01;
    
    const double Umin = -84.0;
    const double Umax = 180.0;
    
    const double Vmin = -105.0;
    const double Vmax = 110.0;

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            color_t c = img_get_pixel(m_thumb256, x, y);
            color_t cwt = img_get_pixel(m_feature_weight_map, x, y);
 
            double L, U, V;
            color_RGBtoLUV(c.r, c.g, c.b, &L, &U, &V);
           
            if (Lmax < L || L < Lmin)
                printf("L: %0.3f\n", L);

            if (Umax < U || U < Umin)
                printf("U: %0.3f\n", U);

            if (Vmax < V || V < Vmin)
                printf("V: %0.3f\n", V);

            double wt = cwt.r / 255.0;

            double L_norm = (L - Lmin) / (Lmax - Lmin);
            double U_norm = (U - Umin) / (Umax - Umin);
            double V_norm = (V - Vmin) / (Vmax - Vmin);

            L_norm = CLAMP(L_norm, 0.0, 1.0) * nbins;
            U_norm = CLAMP(U_norm, 0.0, 1.0) * nbins;
            V_norm = CLAMP(V_norm, 0.0, 1.0) * nbins;

            int L_f = iround(floor(L_norm));
            int U_f = iround(floor(U_norm));
            int V_f = iround(floor(V_norm));

            int L_c = L_f + 1;
            int U_c = U_f + 1;
            int V_c = V_f + 1;

            assert(L_f <= nbins - 1);
            assert(U_f <= nbins - 1);
            assert(V_f <= nbins - 1);

            double t = 0.0;

            t = L_norm - L_f;
            if (L_f == nbins - 1) {
                L_hist[L_f] += wt;
            } else {
                L_hist[L_f] += wt * (1.0 - t);
                L_hist[L_c] += wt * t;
            }

            t = U_norm - U_f;
            if (U_f == nbins - 1) {
                U_hist[U_f] += wt;
            } else {
                U_hist[U_f] += wt * (1.0 - t);
                U_hist[U_c] += wt * t;
            }

            t = V_norm - V_f;
            if (V_f == nbins - 1) {
                V_hist[V_f] += wt;
            } else {
                V_hist[V_f] += wt * (1.0 - t);
                V_hist[V_c] += wt * t;
            }
        }
    }

    /* Normalize bins */
    double L_sum = 0.0, U_sum = 0.0, V_sum = 0.0;
    for (int i = 0; i < nbins; i++) {
        L_sum += L_hist[i];
        U_sum += U_hist[i];
        V_sum += V_hist[i];        
    }

    for (int i = 0; i < nbins; i++) {
        L_hist[i] /= L_sum;
        U_hist[i] /= U_sum;
        V_hist[i] /= V_sum;
    }    
}


void ImageData::ComputeHistogramLUVNorm(int nbins, double *L_hist, 
                                        double *U_hist, double *V_hist)
{
    LoadThumb256();
    
    if (!LoadFeatureWeightMap()) {
        printf("[ImageData::ComputeHistogram] Error: weight map for image %s "
               "not found!\n", m_name);

        return;
    }

    int w = m_thumb256->w;
    int h = m_thumb256->h;

    double Lsum = 0.0;
    double Usum = 0.0;
    double Vsum = 0.0;

    double wt_sum = 0.0;

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            color_t c = img_get_pixel(m_thumb256, x, y);
            color_t cwt = img_get_pixel(m_feature_weight_map, x, y);
 
            double L, U, V;
            color_RGBtoLUV(c.r, c.g, c.b, &L, &U, &V);
           
            double wt = cwt.r / 255.0;

            Lsum += wt * L;
            Usum += wt * U;
            Vsum += wt * V;

            wt_sum += wt;
        }
    }
    
    double Lmean = Lsum / wt_sum;
    double Umean = Usum / wt_sum;
    double Vmean = Vsum / wt_sum;

    Lsum = Usum = Vsum = 0.0;

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            color_t c = img_get_pixel(m_thumb256, x, y);
            color_t cwt = img_get_pixel(m_feature_weight_map, x, y);
 
            double L, U, V;
            color_RGBtoLUV(c.r, c.g, c.b, &L, &U, &V);
           
            double wt = cwt.r / 255.0;

            Lsum += wt * (L - Lmean) * (L - Lmean);
            Usum += wt * (U - Umean) * (U - Umean);
            Vsum += wt * (V - Vmean) * (V - Vmean);
        }
    }

    double Lstd = 2.0 * sqrt(Lsum / wt_sum);
    double Ustd = 2.0 * sqrt(Usum / wt_sum);
    double Vstd = 2.0 * sqrt(Vsum / wt_sum);
    
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            color_t c = img_get_pixel(m_thumb256, x, y);
            color_t cwt = img_get_pixel(m_feature_weight_map, x, y);

            double L, U, V;
            color_RGBtoLUV(c.r, c.g, c.b, &L, &U, &V);

            double L_norm = (L - Lmean) / Lstd + 0.5;
            double U_norm = (U - Umean) / Ustd + 0.5;
            double V_norm = (V - Vmean) / Vstd + 0.5;

            double wt = cwt.r / 255.0;

            L_norm = CLAMP(L_norm, 0.0, 1.0 - 1.0e-3) * nbins;
            U_norm = CLAMP(U_norm, 0.0, 1.0 - 1.0e-3) * nbins;
            V_norm = CLAMP(V_norm, 0.0, 1.0 - 1.0e-3) * nbins;

            int L_f = iround(floor(L_norm));
            int U_f = iround(floor(U_norm));
            int V_f = iround(floor(V_norm));

            int L_c = L_f + 1;
            int U_c = U_f + 1;
            int V_c = V_f + 1;

            assert(L_f <= nbins - 1);
            assert(U_f <= nbins - 1);
            assert(V_f <= nbins - 1);

            double t = 0.0;

            t = L_norm - L_f;
            if (L_f == nbins - 1) {
                L_hist[L_f] += wt;
            } else {
                L_hist[L_f] += wt * (1.0 - t);
                L_hist[L_c] += wt * t;
            }

            t = U_norm - U_f;
            if (U_f == nbins - 1) {
                U_hist[U_f] += wt;
            } else {
                U_hist[U_f] += wt * (1.0 - t);
                U_hist[U_c] += wt * t;
            }

            t = V_norm - V_f;
            if (V_f == nbins - 1) {
                V_hist[V_f] += wt;
            } else {
                V_hist[V_f] += wt * (1.0 - t);
                V_hist[V_c] += wt * t;
            }
        }
    }

    /* Normalize bins */
    double L_sum = 0.0, U_sum = 0.0, V_sum = 0.0;
    for (int i = 0; i < nbins; i++) {
        L_sum += L_hist[i];
        U_sum += U_hist[i];
        V_sum += V_hist[i];        
    }

    for (int i = 0; i < nbins; i++) {
        L_hist[i] /= L_sum;
        U_hist[i] /= U_sum;
        V_hist[i] /= V_sum;
    }    
}
