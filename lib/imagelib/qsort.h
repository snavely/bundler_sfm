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

#ifndef __qsort_h__
#define __qsort_h__

#ifdef __cplusplus
extern "C" {
#endif

/* Set whether we should sort in ascending or descending order */
void qsort_ascending();
void qsort_descending();

/* Sorts the array of doubles `arr' (of length n) and puts the
 * corresponding permutation in `perm' */
void qsort_perm(int n, double *arr, int *perm);

/* Permute the array `arr' given permutation `perm' */
void permute_dbl(int n, double *arr, int *perm);
void permute(int n, int size, void *arr, int *perm);

/* Find the median in a set of doubles */
double median(int n, double *arr);
double median_copy(int n, double *arr);

/* Find the kth element in an unordered list of doubles (changes the
 * array) */
double kth_element(int n, int k, double *arr);
/* Same as above, doesn't change the array */
double kth_element_copy(int n, int k, double *arr);

#ifdef __cplusplus
}
#endif

#endif /* __qsort_h__ */
