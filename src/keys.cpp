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

/* keys.cpp */
/* Class for SIFT keypoints */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef WIN32
#include <ext/hash_map>
#else
#include <hash_map>
#endif

#include "keys.h"

#include "defines.h"

#ifdef __BUNDLER_DISTR__
#include "ANN/ANN.h"
#else
#include "ann_1.1_char/include/ANN/ANN.h"
#endif

int GetNumberOfKeysNormal(FILE *fp)
{
    int num, len;

    if (fscanf(fp, "%d %d", &num, &len) != 2) {
	printf("Invalid keypoint file.\n");
	return 0;
    }

#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    return num;
}

int GetNumberOfKeysGzip(gzFile fp)
{
    int num, len;

    char header[256];
    gzgets(fp, header, 256);

    if (sscanf(header, "%d %d", &num, &len) != 2) {
	printf("Invalid keypoint file.\n");
	return 0;
    }

#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    return num;
}

int GetNumberOfKeysBin(FILE *f)
{
    int num;
    fread(&num, sizeof(int), 1, f);
    
#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    return num;
}

int GetNumberOfKeysBinGzip(gzFile gzf)
{
    int num;
    gzread(gzf, &num, sizeof(int));
    
#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    return num;
}

/* Returns the number of keys in a file */
int GetNumberOfKeys(const char *filename)
{
    FILE *file;

    file = fopen (filename, "r");
    if (! file) {
        /* Try to open a gzipped keyfile */
        char buf[1024];
        sprintf(buf, "%s.gz", filename);
        gzFile gzf = gzopen(buf, "rb");

        if (gzf == NULL) {
            /* Try to open a .bin file */
            sprintf(buf, "%s.bin", filename);
            file = fopen(buf, "rb");
            
            if (file == NULL) {
                /* Try to open a gzipped .bin file */
                sprintf(buf, "%s.bin.gz", filename);
                gzf = gzopen(buf, "rb");
                
                if (gzf == NULL) {
                    printf("Could not open file: %s\n", filename);
                    return 0;
                } else {
                    int n = GetNumberOfKeysBinGzip(gzf);
                    gzclose(gzf);
                    return n;
                }
            } else {
                int n = GetNumberOfKeysBin(file);
                fclose(file);
                return n;
            }
        } else {
            int n = GetNumberOfKeysGzip(gzf);
            gzclose(gzf);
            return n;
        }
    } else {
        int n = GetNumberOfKeysNormal(file);
        fclose(file);
        return n;
    }
}

/* This reads a keypoint file from a given filename and returns the list
 * of keypoints. */
std::vector<KeypointWithDesc> ReadKeyFileWithDesc(const char *filename,
                                                  bool descriptor)
{
    FILE *file;

    file = fopen (filename, "r");
    if (! file) {
        /* Try to file a gzipped keyfile */
        char buf[1024];
        sprintf(buf, "%s.gz", filename);
        gzFile gzf = gzopen(buf, "rb");

        if (gzf == NULL) {
            /* Try to open a .bin file */
            sprintf(buf, "%s.bin", filename);
            file = fopen(buf, "rb");
            
            if (file == NULL) {
                /* Try to open a gzipped .bin file */
                sprintf(buf, "%s.bin.gz", filename);
                gzf = gzopen(buf, "rb");
                
                if (gzf == NULL) {
                    std::vector<KeypointWithDesc> empty;
                    printf("Could not open file: %s\n", filename);
                    return empty;
                } else {
                    std::vector<KeypointWithDesc> kps_desc = 
                        ReadKeysFastBinGzip(gzf, descriptor);
                    gzclose(gzf);
                    return kps_desc;
                }
            } else {
                std::vector<KeypointWithDesc> kps_desc = 
                    ReadKeysFastBin(file, descriptor);
                fclose(file);
                return kps_desc;
            }
        } else {
            std::vector<KeypointWithDesc> kps_desc = 
                ReadKeysFastGzip(gzf, descriptor);
            gzclose(gzf);
            return kps_desc;
        }
    } else {
        std::vector<KeypointWithDesc> kps_desc = ReadKeysFast(file, descriptor);
        fclose(file);
        return kps_desc;
    }
}

std::vector<Keypoint> ReadKeyFile(const char *filename)
{
    std::vector<KeypointWithDesc> kps_d = ReadKeyFileWithDesc(filename, false);

    std::vector<Keypoint> kps;
    int num_keys = (int) kps_d.size();
    kps.resize(num_keys);
    for (int i = 0; i < num_keys; i++) {
        kps[i].m_x = kps_d[i].m_x;
        kps[i].m_y = kps_d[i].m_y;
    }

    kps_d.clear();

    return kps;
}

/* This reads a keypoint file from a given filename and returns the list
 * of keypoints. */
std::vector<KeypointWithScaleRot> 
    ReadKeyFileWithScaleRot(const char *filename, bool descriptor)
{
    FILE *file;
    std::vector<KeypointWithDesc> kps;
    float *scale = NULL, *orient = NULL;

    file = fopen (filename, "r");
    if (! file) {
        /* Try to file a gzipped keyfile */
        char buf[1024];
        sprintf(buf, "%s.gz", filename);
        gzFile gzf = gzopen(buf, "rb");

        if (gzf == NULL) {
            /* Try to open a .bin file */
            sprintf(buf, "%s.bin", filename);
            file = fopen(buf, "rb");
            
            if (file == NULL) {
                /* Try to open a gzipped .bin file */
                sprintf(buf, "%s.bin.gz", filename);
                gzf = gzopen(buf, "rb");
                
                if (gzf == NULL) {
                    std::vector<KeypointWithScaleRot> empty;
                    printf("Could not open file: %s\n", filename);
                    return empty;
                } else {
                    kps = ReadKeysFastBinGzip(gzf, descriptor, &scale, &orient);
                    gzclose(gzf);
                }
            } else {
                kps = ReadKeysFastBin(file, descriptor, &scale, &orient);
                fclose(file);
            }
        } else {
            kps = ReadKeysFastGzip(gzf, descriptor, &scale, &orient);
            gzclose(gzf);
        }
    } else {
        kps = ReadKeysFast(file, descriptor, &scale, &orient);
        fclose(file);
    }

    std::vector<KeypointWithScaleRot> kps_w;
    int num_keys = (int) kps.size();
    kps_w.resize(num_keys);
    for (int i = 0; i < num_keys; i++) {
        kps_w[i].m_x = kps[i].m_x;
        kps_w[i].m_y = kps[i].m_y;
        kps_w[i].m_d = kps[i].m_d;
        kps_w[i].m_scale = scale[i];
        kps_w[i].m_orient = orient[i];
    }

    kps.clear();

    if (scale != NULL)
        delete [] scale;

    if (scale != NULL)
        delete [] orient;

    return kps_w;
}

static char *strchrn(char *str, int c, int n) {
    for (int i = 0; i < n; i++) {
	str = strchr(str, c) + 1;
	if (str == NULL) return NULL;
    }

    return str - 1;
}

/* Read keypoints from the given file pointer and return the list of
 * keypoints.  The file format starts with 2 integers giving the total
 * number of keypoints and the size of descriptor vector for each
 * keypoint (currently assumed to be 128). Then each keypoint is
 * specified by 4 floating point numbers giving subpixel row and
 * column location, scale, and orientation (in radians from -PI to
 * PI).  Then the descriptor vector for each keypoint is given as a
 * list of integers in range [0,255]. */
std::vector<Keypoint> ReadKeys(FILE *fp, bool descriptor)
{
    int i, j, num, len, val;

    std::vector<Keypoint> kps;

    if (fscanf(fp, "%d %d", &num, &len) != 2) {
	printf("Invalid keypoint file beginning.");
	return kps;
    }

#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    if (len != 128) {
	printf("Keypoint descriptor length invalid (should be 128).");
	return kps;
    }

    for (i = 0; i < num; i++) {
	/* Allocate memory for the keypoint. */
	unsigned char *d = new unsigned char[len];
	float x, y, scale, ori;

	if (fscanf(fp, "%f %f %f %f", &y, &x, &scale, &ori) != 4) {
	    printf("Invalid keypoint file format.");
	    return kps;
	}

	for (j = 0; j < len; j++) {
	    if (fscanf(fp, "%d", &val) != 1 || val < 0 || val > 255) {
		printf("Invalid keypoint file value.");
		return kps;
	    }
	    d[j] = (unsigned char) val;
	}

	if (descriptor) {
	    kps.push_back(KeypointWithDesc(x, y, d));
	} else {
	    delete [] d;
	    kps.push_back(Keypoint(x, y));
	}
    }

    return kps;
}

/* Read keys more quickly */
std::vector<KeypointWithDesc> ReadKeysFast(FILE *fp, bool descriptor,
                                           float **scales, float **orients)
{
    int i, j, num, len;

    std::vector<KeypointWithDesc> kps;

    if (fscanf(fp, "%d %d", &num, &len) != 2) {
	printf("Invalid keypoint file beginning.");
	return kps;
    }

#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    if (len != 128) {
	printf("Keypoint descriptor length invalid (should be 128).");
	return kps;
    }

    kps.resize(num);

    if (num > 0 && scales != NULL) {
        *scales = new float[num];
    }

    if (num > 0 && orients != NULL) {
        *orients = new float[num];
    }

    for (i = 0; i < num; i++) {
	/* Allocate memory for the keypoint. */
	float x, y, scale, ori;

	if (fscanf(fp, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
	    printf("Invalid keypoint file format.");
	    return kps;
	}

        if (scales != NULL) {
            (*scales)[i] = scale;
        }
        
        if (orients != NULL) {
            (*orients)[i] = ori;
        }

	char buf[1024];

	/* Allocate memory for the keypoint. */
	unsigned char *d = NULL;
	
	if (descriptor)
	    d = new unsigned char[len];

	int start = 0;
	for (int line = 0; line < 7; line++) {
	    fgets(buf, 1024, fp);

	    if (!descriptor) continue;

	    short int p[20];

	    if (line < 6) {
		sscanf(buf, 
		       "%hu %hu %hu %hu %hu %hu %hu %hu %hu %hu "
		       "%hu %hu %hu %hu %hu %hu %hu %hu %hu %hu", 
		       p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, p+8, p+9, 
		       p+10, p+11, p+12, p+13, p+14, 
		       p+15, p+16, p+17, p+18, p+19);

		for (j = 0; j < 20; j++)
		    d[start + j] = p[j];

		start += 20;
	    } else {
		sscanf(buf, 
		       "%hu %hu %hu %hu %hu %hu %hu %hu",
		       p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7);

		for (j = 0; j < 8; j++)
		    d[start + j] = p[j];
	    }
	}

        // kps.push_back(KeypointWithDesc(x, y, d));
        kps[i] = KeypointWithDesc(x, y, d);
    }

    return kps;
}

std::vector<KeypointWithDesc> ReadKeysFastGzip(gzFile fp, bool descriptor,
                                               float **scales, float **orients)
{
    int i, j, num, len;

    std::vector<KeypointWithDesc> kps;
    char header[256];
    gzgets(fp, header, 256);

    if (sscanf(header, "%d %d", &num, &len) != 2) {
	printf("Invalid keypoint file.\n");
	return kps;
    }

#ifdef KEY_LIMIT
    num = MIN(num, 65536); // we'll store at most 65536 features per
                           // image
#endif /* KEY_LIMIT */

    if (len != 128) {
	printf("Keypoint descriptor length invalid (should be 128).");
	return kps;
    }

    kps.resize(num);

    if (num > 0 && scales != NULL) {
        *scales = new float[num];
    }

    if (num > 0 && orients != NULL) {
        *orients = new float[num];
    }

    for (i = 0; i < num; i++) {
	/* Allocate memory for the keypoint. */
	float x, y, scale, ori;
        char buf[1024];
        gzgets(fp, buf, 1024);

	if (sscanf(buf, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
	    printf("Invalid keypoint file format.");
	    return kps;
	}

        if (scales != NULL) {
            (*scales)[i] = scale;
        }
        
        if (orients != NULL) {
            (*orients)[i] = ori;
        }

	/* Allocate memory for the keypoint. */
	unsigned char *d = NULL;
	
	if (descriptor)
	    d = new unsigned char[len];

	int start = 0;
	for (int line = 0; line < 7; line++) {
	    gzgets(fp, buf, 1024);

	    if (!descriptor) continue;

	    short int p[20];

	    if (line < 6) {
		sscanf(buf, 
		       "%hu %hu %hu %hu %hu %hu %hu %hu %hu %hu "
		       "%hu %hu %hu %hu %hu %hu %hu %hu %hu %hu", 
		       p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, p+8, p+9, 
		       p+10, p+11, p+12, p+13, p+14, 
		       p+15, p+16, p+17, p+18, p+19);

		for (j = 0; j < 20; j++)
		    d[start + j] = p[j];

		start += 20;
	    } else {
		sscanf(buf, 
		       "%hu %hu %hu %hu %hu %hu %hu %hu",
		       p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7);

		for (j = 0; j < 8; j++)
		    d[start + j] = p[j];
	    }
	}

        // kps.push_back(KeypointWithDesc(x, y, d));
        kps[i] = KeypointWithDesc(x, y, d);
    }

    return kps;
}

/* Read keys from binary file */
std::vector<KeypointWithDesc> ReadKeysFastBin(FILE *fp, bool descriptor,
                                              float **scales, 
                                              float **orients)
{
    int num_keys;
    fread(&num_keys, sizeof(int), 1, fp);

    std::vector<KeypointWithDesc> keys;
    keys.resize(num_keys);

    keypt_t *info;
    unsigned char *d;

    info = new keypt_t[num_keys];
    
    fread(info, sizeof(keypt_t), num_keys, fp);

    if (scales != NULL)
        *scales = new float[num_keys];
    
    if (orients != NULL)
        *orients = new float[num_keys];

    for (int i = 0; i < num_keys; i++) {
        keys[i].m_x = info[i].x;
        keys[i].m_y = info[i].y;
        
        if (scales != NULL)
            (*scales)[i] = info[i].scale;
        if (orients != NULL)
            (*orients)[i] = info[i].orient;
    }

    delete [] info;

    if (!descriptor)
        return keys;
    
    d = new unsigned char [128 * num_keys];

    fread(d, sizeof(unsigned char), 128 * num_keys, fp);
    
    for (int i = 0; i < num_keys; i++) {
        keys[i].m_d = d + 128 * i;
    }

    return keys;
}

/* Read keys from gzipped binary file */
std::vector<KeypointWithDesc> ReadKeysFastBinGzip(gzFile fp, bool descriptor,
                                                  float **scales, 
                                                  float **orients)
{
    int num_keys;
    gzread(fp, &num_keys, sizeof(int));

    std::vector<KeypointWithDesc> keys;
    keys.resize(num_keys);

    keypt_t *info;
    unsigned char *d;

    info = new keypt_t[num_keys];
    
    gzread(fp, info, sizeof(keypt_t) * num_keys);

    if (scales != NULL)
        *scales = new float[num_keys];
    
    if (orients != NULL)
        *orients = new float[num_keys];

    for (int i = 0; i < num_keys; i++) {
        keys[i].m_x = info[i].x;
        keys[i].m_y = info[i].y;
        
        if (scales != NULL)
            (*scales)[i] = info[i].scale;
        if (orients != NULL)
            (*orients)[i] = info[i].orient;
    }

    delete [] info;

    if (!descriptor)
        return keys;
    
    d = new unsigned char [128 * num_keys];

    gzread(fp, d, sizeof(unsigned char) * 128 * num_keys);
    
    for (int i = 0; i < num_keys; i++) {
        keys[i].m_d = d + 128 * i;
    }

    return keys;
}

#if 0
ANNkd_tree *CreateSearchTree(const std::vector<KeypointWithDesc> &k, 
                             bool spatial, double alpha) 
{
    /* Create a new array of points */
    int num_pts = (int) k.size();

    int dim = 128;
    if (spatial) dim = 130;

    ANNpointArray pts = annAllocPts(num_pts, dim);

    int offset = 0;
    if (spatial) offset = 2;
     
    for (int i = 0; i < num_pts; i++) {
        int j;

        assert(k[i].m_d != NULL);
        
        if (spatial) {
            pts[i][0] = alpha * k[i].m_x;
            pts[i][1] = alpha * k[i].m_y;
        }
        
        for (j = 0; j < 128; j++)
            pts[i][j+offset] = k[i].m_d[j];
    }
    
    /* Create a search tree for k2 */
    ANNkd_tree *tree = new ANNkd_tree(pts, num_pts, dim, 4);

    // annDeallocPts(pts);

    return tree;
}
#endif

ann_1_1_char::ANNkd_tree
    *CreateSearchTreeChar(const std::vector<KeypointWithDesc> &k) 
{
    /* Create a new array of points */
    int num_pts = (int) k.size();

    int dim = 128;

    ann_1_1_char::ANNpointArray pts = ann_1_1_char::annAllocPts(num_pts, dim);

    int offset = 0;
     
    for (int i = 0; i < num_pts; i++) {
        int j;
        
        for (j = 0; j < 128; j++)
            pts[i][j+offset] = k[i].m_d[j];
    }
    
    /* Create a search tree for k2 */
    ann_1_1_char::ANNkd_tree *tree = 
        new ann_1_1_char::ANNkd_tree(pts, num_pts, dim, 4);

    // ann_1_1_char::annDeallocPts(pts);

    return tree;
}

/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatch> MatchKeys(const std::vector<KeypointWithDesc> &k1, 
				     const std::vector<KeypointWithDesc> &k2, 
				     bool registered, double ratio) 
{
    ann_1_1_char::annMaxPtsVisit(200);

    int num_pts = 0;
    std::vector<KeypointMatch> matches;

    int *registered_idxs = NULL;

    if (!registered) {
	num_pts = (int) k2.size();
    } else {
	registered_idxs = new int[(int) k2.size()];
	for (int i = 0; i < (int) k2.size(); i++) {
	    if (k2[i].m_extra >= 0) {
		registered_idxs[num_pts] = i;
		num_pts++;
	    }
	}
    }

    /* Create a new array of points */
    ann_1_1_char::ANNpointArray pts = ann_1_1_char::annAllocPts(num_pts, 128);

    if (!registered) {
	for (int i = 0; i < num_pts; i++) {
	    int j;

	    for (j = 0; j < 128; j++) {
		pts[i][j] = k2[i].m_d[j];
	    }
	}
    } else {
	for (int i = 0; i < num_pts; i++) {
	    int j;
	    int idx = registered_idxs[i];

	    for (j = 0; j < 128; j++) {
		pts[i][j] = k2[idx].m_d[j];
	    }
	}	
    }
    
    // clock_t start = clock();
    /* Create a search tree for k2 */
    ann_1_1_char::ANNkd_tree *tree = new ann_1_1_char::ANNkd_tree(pts, num_pts, 128, 4);
    // clock_t end = clock();
    
    // printf("Building tree took %0.3fs\n", 
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    /* Now do the search */
    ann_1_1_char::ANNpoint query = ann_1_1_char::annAllocPt(128);
    // start = clock();
    for (int i = 0; i < (int) k1.size(); i++) {
	int j;

	for (j = 0; j < 128; j++) {
	    query[j] = k1[i].m_d[j];
	}

	ann_1_1_char::ANNidx nn_idx[2];
	ann_1_1_char::ANNdist dist[2];

	tree->annkPriSearch(query, 2, nn_idx, dist, 0.0);

	if (sqrt(((double) dist[0]) / ((double) dist[1])) <= ratio) {
	    if (!registered) {
		matches.push_back(KeypointMatch(i, nn_idx[0]));
	    } else {
		KeypointMatch match = 
		    KeypointMatch(i, registered_idxs[nn_idx[0]]);
		matches.push_back(match);
	    }
	}
    }
    // end = clock();
    // printf("Searching tree took %0.3fs\n",
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    int num_matches = (int) matches.size();

    printf("[MatchKeys] Found %d matches\n", num_matches);

    /* Cleanup */
    ann_1_1_char::annDeallocPts(pts);
    ann_1_1_char::annDeallocPt(query);

    delete tree;

    return matches;
}

/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatchWithScore> 
    MatchKeysWithScore(const std::vector<KeypointWithDesc> &k1, 
                       const std::vector<KeypointWithDesc> &k2,
                       bool registered, 
                       double ratio)
{
    ann_1_1_char::annMaxPtsVisit(200);

    int num_pts = 0;
    std::vector<KeypointMatchWithScore> matches;

    int *registered_idxs = NULL;

    if (!registered) {
	num_pts = (int) k2.size();
    } else {
	registered_idxs = new int[(int) k2.size()];
	for (int i = 0; i < (int) k2.size(); i++) {
	    if (k2[i].m_extra >= 0) {
		registered_idxs[num_pts] = i;
		num_pts++;
	    }
	}
    }

    /* Create a new array of points */
    ann_1_1_char::ANNpointArray pts = ann_1_1_char::annAllocPts(num_pts, 128);

    if (!registered) {
	for (int i = 0; i < num_pts; i++) {
	    int j;

	    for (j = 0; j < 128; j++) {
		pts[i][j] = k2[i].m_d[j];
	    }
	}
    } else {
	for (int i = 0; i < num_pts; i++) {
	    int j;
	    int idx = registered_idxs[i];

	    for (j = 0; j < 128; j++) {
		pts[i][j] = k2[idx].m_d[j];
	    }
	}	
    }
    
    // clock_t start = clock();
    /* Create a search tree for k2 */
    ann_1_1_char::ANNkd_tree *tree = new ann_1_1_char::ANNkd_tree(pts, num_pts, 128, 4);
    // clock_t end = clock();
    
    // printf("Building tree took %0.3fs\n", 
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    /* Now do the search */
    ann_1_1_char::ANNpoint query = ann_1_1_char::annAllocPt(128);
    // start = clock();
    for (int i = 0; i < (int) k1.size(); i++) {
	int j;

	for (j = 0; j < 128; j++) {
	    query[j] = k1[i].m_d[j];
	}

	ann_1_1_char::ANNidx nn_idx[2];
	ann_1_1_char::ANNdist dist[2];

	tree->annkPriSearch(query, 2, nn_idx, dist, 0.0);

	if (sqrt(((double) dist[0]) / ((double) dist[1])) <= ratio) {
	    if (!registered) {
                KeypointMatchWithScore match = 
                    KeypointMatchWithScore(i, nn_idx[0], (float) dist[0]);
                matches.push_back(match);
	    } else {
		KeypointMatchWithScore match = 
		    KeypointMatchWithScore(i, registered_idxs[nn_idx[0]], 
                                           (float) dist[0]);
		matches.push_back(match);
	    }
	}
    }
    // end = clock();
    // printf("Searching tree took %0.3fs\n",
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    int num_matches = (int) matches.size();

    printf("[MatchKeysWithScore] Found %d matches\n", num_matches);

    /* Cleanup */
    ann_1_1_char::annDeallocPts(pts);
    ann_1_1_char::annDeallocPt(query);

    delete tree;

    return matches;    
}

/* Prune matches so that they are 1:1 */
std::vector<KeypointMatchWithScore> 
    PruneMatchesWithScore(const std::vector<KeypointMatchWithScore> &matches)
{
#ifndef WIN32
    __gnu_cxx::hash_map<int, float> key_hash;
    __gnu_cxx::hash_map<int, int> map;
#else
    stdext::hash_map<int, float> key_hash;
    stdext::hash_map<int, int> map;
#endif

    int num_matches = (int) matches.size();
    
    for (int i = 0; i < num_matches; i++) {
        int idx1 = matches[i].m_idx1;
        int idx2 = matches[i].m_idx2;

        if (key_hash.find(idx2) == key_hash.end()) {
            /* Insert the new element */
            key_hash[idx2] = matches[i].m_score;
            map[idx2] = idx1;
        } else {
            float old = key_hash[idx2];
            if (old > matches[i].m_score) {
                /* Replace the old entry */
                key_hash[idx2] = matches[i].m_score;
                map[idx2] = idx1;
            }
        }
    }

    std::vector<KeypointMatchWithScore> matches_new;
    /* Now go through the list again, building a new list */
    for (int i = 0; i < num_matches; i++) {
        int idx1 = matches[i].m_idx1;
        int idx2 = matches[i].m_idx2;

        if (map[idx2] == idx1) {
            matches_new.push_back(KeypointMatchWithScore(idx1, idx2, 
                                                         key_hash[idx2]));
        }
    }

    return matches_new;
}

/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatch>
    MatchKeysExhaustive(const std::vector<KeypointWithDesc> &k1, 
                        const std::vector<KeypointWithDesc> &k2, 
                        bool registered, double ratio) 
{
    int num_pts = 0;
    std::vector<KeypointMatch> matches;

    int *registered_idxs = NULL;

    if (!registered) {
	num_pts = (int) k2.size();
    } else {
	registered_idxs = new int[(int) k2.size()];
	for (int i = 0; i < (int) k2.size(); i++) {
	    if (k2[i].m_extra >= 0) {
		registered_idxs[num_pts] = i;
		num_pts++;
	    }
	}
    }

    /* Create a new array of points */
    ann_1_1_char::ANNpointArray pts = ann_1_1_char::annAllocPts(num_pts, 128);

    if (!registered) {
	for (int i = 0; i < num_pts; i++) {
	    int j;

	    for (j = 0; j < 128; j++) {
		pts[i][j] = k2[i].m_d[j];
	    }
	}
    } else {
	for (int i = 0; i < num_pts; i++) {
	    int j;
	    int idx = registered_idxs[i];

	    for (j = 0; j < 128; j++) {
		pts[i][j] = k2[idx].m_d[j];
	    }
	}	
    }
    
    // clock_t start = clock();
    /* Create a search tree for k2 */
    ann_1_1_char::ANNkd_tree *tree = 
        new ann_1_1_char::ANNkd_tree(pts, num_pts, 128, 4);
    // clock_t end = clock();
    
    // printf("Building tree took %0.3fs\n", 
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    /* Now do the search */
    ann_1_1_char::ANNpoint query = ann_1_1_char::annAllocPt(128);
    // start = clock();
    for (int i = 0; i < (int) k1.size(); i++) {
	int j;

	for (j = 0; j < 128; j++) {
	    query[j] = k1[i].m_d[j];
	}

	ann_1_1_char::ANNidx nn_idx[2];
	ann_1_1_char::ANNdist dist[2];

	tree->annkSearch(query, 2, nn_idx, dist, 0.0);

	if (sqrt(((double) dist[0]) / ((double) dist[1])) <= ratio) {
	    if (!registered) {
		matches.push_back(KeypointMatch(i, nn_idx[0]));
	    } else {
		KeypointMatch match = 
		    KeypointMatch(i, registered_idxs[nn_idx[0]]);
		matches.push_back(match);
	    }
	}
    }
    // end = clock();
    // printf("Searching tree took %0.3fs\n",
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    int num_matches = (int) matches.size();

    printf("[MatchKeys] Found %d matches\n", num_matches);

    /* Cleanup */
    ann_1_1_char::annDeallocPts(pts);
    ann_1_1_char::annDeallocPt(query);

    delete tree;

    return matches;
}
