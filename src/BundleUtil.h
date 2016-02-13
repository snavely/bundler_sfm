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

/* BundleUtil.h */
/* Various utility routines */

#ifndef ___bundle_util_h___
#define ___bundle_util_h___

#include "image.h"

#ifndef WIN32
#include <ext/hash_map>
#include <ext/hash_set>
#else
#include <hash_map>
#include <hash_set>
#endif

#include <string>

typedef struct  {
    std::string name;
    bool is_fisheye;
    float focal;
} image_t;

img_t *RescaleImage(img_t *img, double scale);
img_t *RescaleImage(img_t *img, int max_dim, double &scale);

void GetRotationFromSpherical(double theta, double phi, double *R);

/* Normalize a patch for mean and variance (in place) */
void NormalizePatchMeanVariance(int w, int h, double *patch);

/* Return the intersection of two int vectors */
std::vector<int> GetVectorIntersection(const std::vector<int> &v1,
				       const std::vector<int> &v2);

/* Returns true if the intersection of two int vectors is non empty*/
bool VectorIntersectionNonEmpty(const std::vector<int> &v1,
									   const std::vector<int> &v2);

std::vector<std::pair<int, int> >
    GetArrayIntersectionSorted(int m, int n, const int *a1, const int *a2);

std::vector<std::pair<int, int> >
    GetArrayIntersectionUnsorted(int m, int n, const int *a1, const int *a2);

void choose(int n, int k, int *arr);

void generate_permutation(int n, int *arr);

void Tokenize(const std::string &str,
              std::vector<std::string> &tokens,
              const std::string &delimiters);

bool ReadListFile(const char *list_file, std::vector<image_t> &images, 
                  const std::string &prefix = ".");

bool FileExists(const char *filename);

void ConvertRGBtoLAB(double r, double g, double b, 
                     double &L, double &A, double &B);

void ConvertLABtoRGB(double L, double A, double B,
                     double &r, double &g, double &b);

extern double RGBtoLAB[9];
extern double RGBtoLAB_4x4[16];
extern double LABtoRGB[9];

#ifndef WIN32
typedef __gnu_cxx::hash_set<int> HashSetInt;
#else
typedef stdext::hash_set<int> HashSetInt;
#endif

#endif /* __bundle_util_h__ */
