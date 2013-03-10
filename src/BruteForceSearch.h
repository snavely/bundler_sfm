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

/* BruteForceSearch.h */
/* Brute force nearest-neighbor search data structure */

#ifndef __brute_force_search_h__
#define __brute_force_search_h__

#include <stdlib.h>
#include "vector.h"

class BruteForceSearch 
{
public:
    BruteForceSearch() : m_num_points(0), m_points(NULL) { }
    BruteForceSearch(int n, v3_t *pts);
    ~BruteForceSearch() { if (m_points != NULL) delete [] m_points; }
    
    void GetClosestPoints(v3_t query, int num_points, 
			  int *idxs, double *dists);

    int m_num_points;
    v3_t *m_points;
};

#endif /* __brute_force_search_h__ */
