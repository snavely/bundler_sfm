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

/* KeyMatchFull.cpp */
/* Read in keys, match, write results to a file */

#include <assert.h>
#include <time.h>
#include <string.h>

#include "keys2a.h"

int main(int argc, char **argv) {
    char *list_in;
    char *file_out;
    double ratio;
    
    if (argc != 3) {
	printf("Usage: %s <list.txt> <outfile>\n", argv[0]);
	return -1;
    }
    
    list_in = argv[1];
    ratio = 0.6;
    file_out = argv[2];

    clock_t start = clock();

    unsigned char **keys;
    int *num_keys;

    /* Read the list of files */
    std::vector<std::string> key_files;
    
    FILE *f = fopen(list_in, "r");
    if (f == NULL) {
        printf("Error opening file %s for reading\n", list_in);
        return 1;
    }

    char buf[512];
    while (fgets(buf, 512, f)) {
        /* Remove trailing newline */
        if (buf[strlen(buf) - 1] == '\n')
            buf[strlen(buf) - 1] = 0;
        
        key_files.push_back(std::string(buf));
    }

    fclose(f);

    f = fopen(file_out, "w");
    assert(f != NULL);

    int num_images = (int) key_files.size();

    keys = new unsigned char *[num_images];
    num_keys = new int[num_images];

    /* Read all keys */
    for (int i = 0; i < num_images; i++) {
        keys[i] = NULL;
        num_keys[i] = ReadKeyFile(key_files[i].c_str(), keys+i);
    }

    clock_t end = clock();    
    printf("[KeyMatchFull] Reading keys took %0.3fs\n", 
           (end - start) / ((double) CLOCKS_PER_SEC));
    
    for (int i = 0; i < num_images; i++) {
        if (num_keys[i] == 0)
            continue;

        printf("[KeyMatchFull] Matching to image %d\n", i);

        start = clock();

        /* Create a tree from the keys */
        ANNkd_tree *tree = CreateSearchTree(num_keys[i], keys[i]);

        for (int j = 0; j < i; j++) {
            if (num_keys[j] == 0)
                continue;

            /* Compute likely matches between two sets of keypoints */
            std::vector<KeypointMatch> matches = 
                MatchKeys(num_keys[j], keys[j], tree, ratio);
            
            int num_matches = (int) matches.size();

            if (num_matches >= 16) {
                /* Write the pair */
                fprintf(f, "%d %d\n", j, i);

                /* Write the number of matches */
                fprintf(f, "%d\n", (int) matches.size());

                for (int i = 0; i < num_matches; i++) {
                    fprintf(f, "%d %d\n", 
                            matches[i].m_idx1, matches[i].m_idx2);
                }
            }
        }

        end = clock();    
        printf("[KeyMatchFull] Matching took %0.3fs\n", 
               (end - start) / ((double) CLOCKS_PER_SEC));
        fflush(stdout);

        // annDeallocPts(tree->pts);
        delete tree;
    }
    
    /* Free keypoints */
    for (int i = 0; i < num_images; i++) {
        if (keys[i] != NULL)
            delete [] keys[i];
    }
    delete [] keys;
    delete [] num_keys;
    
    fclose(f);
    return 0;
}
