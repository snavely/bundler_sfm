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

/* BundleUtil.cpp */
/* Various utility routines */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#ifndef WIN32
#include <ext/hash_set>
#include <ext/hash_map>
#else
#include <hash_set>
#include <hash_map>
#endif

#include <vector.h>

#include "BundleUtil.h"

#include "defines.h"
#include "filter.h"
#include "matrix.h"
#include "resample.h"

img_t *RescaleImage(img_t *img, double scale) 
{
    img_t *blur = img_smooth(img, 0.35 / scale, 0);
    trans2D_t *T = new_scaling_transform(scale, scale);
    img_t *scaled = img_resample_bbox(blur, T);

    img_free(blur);
    transform_free(T);

    return scaled;
}

img_t *RescaleImage(img_t *img, int max_dim, double &scale) 
{
    int dim = MAX(img->w, img->h);
    double ratio = (double) dim / (double) max_dim;

    if (dim <= max_dim) {
        return img_copy(img);
    }

    scale = 1.0 / ratio;

    return RescaleImage(img, scale);
}

void GetRotationFromSpherical(double theta, double phi, double *R) 
{
    /* Compute the direction we're looking */
    double v[3];
    
    v[0] = cos(theta) * sin(phi);
    v[1] = cos(phi);
    v[2] = sin(theta) * sin(phi);

    /* Compute the new up vector */
    double phi_up = phi - 0.5 * M_PI;
    double theta_up = theta;
    double up[3];
    
    up[0] = cos(theta_up) * sin(phi_up);
    up[1] = cos(phi_up);
    up[2] = sin(theta_up) * sin(phi_up);

    double x_axis[3];
    matrix_cross(v, up, x_axis);
    
    memcpy(R + 0, x_axis, 3 * sizeof(double));
    memcpy(R + 3, up, 3 * sizeof(double));
    memcpy(R + 6, v, 3 * sizeof(double));
}

/* Normalize a patch for mean and variance (in place) */
void NormalizePatchMeanVariance(int w, int h, double *patch)
{
    int nelems = w * h;

    double mean = 0.0;    
    int num_used = 0;

    for (int i = 0; i < nelems; i++) {
	if (patch[i] == DBL_MAX) continue;
	mean += patch[i];
	num_used++;
    }

    if (num_used == 0.0)
	return;

    mean /= num_used;

    double variance = 0.0;
    
    for (int i = 0; i < nelems; i++) {
	if (patch[i] == DBL_MAX) continue;

	double d = patch[i] - mean;
	variance += d * d;
    }

    variance /= num_used;
    variance = sqrt(variance);

    for (int i = 0; i < nelems; i++) {
	if (patch[i] == DBL_MAX) continue;
	patch[i] = 32.0 * (patch[i] - mean) / variance;
    }
}

/* Return the intersection of two int vectors */
std::vector<int> GetVectorIntersection(const std::vector<int> &v1,
				       const std::vector<int> &v2)
{
#ifndef WIN32
	__gnu_cxx::hash_set<int> seen;
#else
	stdext::hash_set<int> seen;
#endif

    int v1_size = (int) v1.size();
    int v2_size = (int) v2.size();

    std::vector<int> intersection;

    for (int i = 0; i < v1_size; i++)
	seen.insert(v1[i]);

    for (int i = 0; i < v2_size; i++) {
	if (seen.find(v2[i]) != seen.end())
	    intersection.push_back(v2[i]);
    }

    seen.clear();

    return intersection;
}

/* Return the intersection of two int vectors */
bool VectorIntersectionNonEmpty(const std::vector<int> &v1,
									   const std::vector<int> &v2)
{
#ifndef WIN32
	__gnu_cxx::hash_set<int> seen;
#else
	stdext::hash_set<int> seen;
#endif

	int v1_size = (int) v1.size();
	int v2_size = (int) v2.size();


	for (int i = 0; i < v1_size; i++)
		seen.insert(v1[i]);

	for (int i = 0; i < v2_size; i++) {
		if (seen.find(v2[i]) != seen.end())
			return true;
	}

	seen.clear();

	return false;
}

std::vector<std::pair<int, int> > 
  GetArrayIntersection(int m, int n, const int *a1, const int *a2) 
{
    int c1 = 0, c2 = 0;
    std::vector<std::pair<int,int> > intersection;

    intersection.reserve(MIN(m, n));

    int count = 0;
    while (c1 < m && c2 < n) {
        while (a1[c1] < a2[c2])
            c1++;

        if (a1[c1] == a2[c2]) {
	    intersection.push_back(std::pair<int,int> (c1, c2));
            count++;
        }

        c2++;
    }

    return intersection;
}

std::vector<std::pair<int, int> > 
GetArrayIntersectionSorted(int m, int n, const int *a1, const int *a2) 
{
    int c1 = 0, c2 = 0;
    std::vector<std::pair<int,int> > intersection;

    intersection.reserve(MIN(m, n));

    int count = 0;
    while (c1 < m && c2 < n) {
        while (a1[c1] < a2[c2] && c1 < m)
            c1++;

        if (c1 >= m)
            break;

        if (a1[c1] == a2[c2]) {
            intersection.push_back(std::pair<int,int> (c1, c2));
            count++;
        }

        c2++;
    }

    return intersection;
}

std::vector<std::pair<int, int> > 
   GetArrayIntersectionUnsorted(int m, int n, const int *a1, const int *a2)
{
#ifndef WIN32
    __gnu_cxx::hash_map<int,int> seen;
#else
    stdext::hash_map<int,int> seen;
#endif

    std::vector<std::pair<int,int> > isect;
    isect.reserve(MIN(m, n));

    for (int i = 0; i < n; i++)
        seen[a2[i]] = i;

    for (int i = 0; i < m; i++) {
        if (seen.find(a1[i]) != seen.end())
            isect.push_back(std::pair<int,int> (i, seen[a1[i]]));
    }

    return isect;    
}

void Tokenize(const std::string& str,
              std::vector<std::string>& tokens,
              const std::string& delimiters)
{
    /* Skip delimiters at beginning */
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    /* Find first "non-delimiter" */
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        /* Found a token, add it to the vector */
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        /* Skip delimiters.  Note the "not_of" */
        lastPos = str.find_first_not_of(delimiters, pos);
        /* Find next "non-delimiter" */
        pos = str.find_first_of(delimiters, lastPos);
    }
}

/* Read a list file */
bool ReadListFile(const char *list_file, std::vector<image_t> &images, 
                  const std::string &prefix)
{
    images.clear();

    FILE *f = fopen(list_file, "r");
    
    if (f == NULL) {
        printf("[ReadListFile] ERROR opening file %s for reading!\n",
               list_file);

        return false;
    }

    char buf[1024];
    while (fgets(buf, 1024, f) != NULL) {
        std::vector<std::string> tokens;
        Tokenize(std::string(buf), tokens, std::string(" \n"));

        image_t image;
        image.name = prefix + "/" + tokens[0];
        if ((int) tokens.size() == 1) {
            image.is_fisheye = false;
            image.focal = 0.0;
        } else {
	    image.is_fisheye = (tokens[1] == "1");
            image.focal = atoi(tokens[2].c_str());
        }

        images.push_back(image);
    }

    fclose(f);

    return true;
}

bool FileExists(const char *filename)
{
    FILE *f = fopen(filename, "r");

    if (f == NULL)
        return false;
    
    fclose(f);

    return true;
}
     
void choose(int n, int k, int *arr)
{
    int i;
    
    if (k > n) {
        printf("[choose] Error: k > n\n");
        return;
    }

    for (i = 0; i < k; i++) {
        while (1) {
            int idx = rand() % n;
            int j, redo = 0;

            for (j = 0; j < i; j++) {
                if (idx == arr[j]) {
                    redo = 1;
                    break;
                }
            }

            if (!redo) {
                arr[i] = idx;
                break;
            }
        }
    }
}

static double RGBtoLMS[9] = { 0.3811, 0.5783, 0.0402, 
                              0.1967, 0.7244, 0.0782, 
                              0.0241, 0.1288, 0.8444 };

static double LMStoLAB[9] = { 0.5774, 0.5774, 0.5774, 
                              0.4082, 0.4082, -0.8165, 
                              0.7071, -0.7071, 0.000 };

static double LABtoLMS[9] = { 0.5774, 0.4082, 0.7071, 
                              0.5774, 0.4082, -0.7071, 
                              0.5774, -0.8165, 0.0000 };    

static double LMStoRGB[9] = { 4.4687, -3.5887, 0.1196, 
                              -1.2197, 2.3831, -0.1626, 
                              0.0585, -0.2611, 1.2057 };

double RGBtoLAB[9] = { 0.4845,    0.4277,   -0.2522,
                       0.4646,    0.3540,   -0.4779,
                       0.6636,   -0.5306,   -0.0913 };

double RGBtoLAB_4x4[16] = { 0.4845,    0.4277,   -0.2522, 0.0,
                            0.4646,    0.3540,   -0.4779, 0.0,
                            0.6636,   -0.5306,   -0.0913, 0.0,
                            0.0, 0.0, 0.0, 1.0 };

double LABtoRGB[9] = { 2.1237,   -1.2840,    0.8552,
                       2.0410,   -0.9147,   -0.8499,
                       3.5761,   -4.0179,    0.2018 };
    



void ConvertRGBtoLAB(double r, double g, double b, 
                     double &L, double &A, double &B)
{
    double rgb[3] = { r, g, b };
    double lms[3], lab[3];

    matrix_product(3, 3, 3, 1, RGBtoLMS, rgb, lms);
    matrix_product(3, 3, 3, 1, LMStoLAB, lms, lab);
    
    L = lab[0];
    A = lab[1];
    B = lab[2];
}

void ConvertLABtoRGB(double L, double A, double B,
                     double &r, double &g, double &b)
{
    double lab[3] = { L, A, B };
    double lms[3], rgb[3];
    
    matrix_product(3, 3, 3, 1, LABtoLMS, lab, lms);
    matrix_product(3, 3, 3, 1, LMStoRGB, lms, rgb);
    
    r = rgb[0];
    g = rgb[1];
    b = rgb[2];
}

void generate_permutation(int n, int *arr)
{
    /* Fill list from 0 to n-1 */
    std::vector<int> a;

    for (int i = 0; i < n; i++) {
        a.push_back(i);
    }

    std::vector<int> b;

    /* Shuffle array */
    for (int i = 0; i < n; i++) {
        int idx = rand() % (n - i);
        b.push_back(a[idx]);
        a.erase(a.begin() + idx);
    }

    for (int i = 0; i < n; i++) {
        arr[i] = b[i];
    }

    b.clear();
}
