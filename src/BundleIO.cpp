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

/* BundleIO.cpp */
/* Routines for reading / writing data structures */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "BaseApp.h"
#include "BundleUtil.h"
#include "LoadJPEG.h"

#include "defines.h"
#include "horn.h"
#include "matrix.h"
#include "qsort.h"
#include "util.h"



#define MIN_MATCHES 10

/* Load keys from a file */
void BaseApp::LoadKeys(bool descriptor)
{
    printf("[LoadKeys] Loading keys...\n");

    clock_t start = clock();

    int num_images = GetNumImages();

    for (int i = 0; i < num_images; i++) {
        printf("[LoadKeys] Loading keys from image %d...\n", i);
        fflush(stdout);
	m_image_data[i].LoadKeys(descriptor);
    }

    clock_t end = clock();
    
    printf("[LoadKeys] Loaded keys in %0.3fs\n", 
	   (end - start) / (double) CLOCKS_PER_SEC);
}

void BaseApp::ReadMatchFile(int i, int j)
{
    /* Read matches from a file */
    char buf[256];
    sprintf(buf, "%s/match-%03d-%03d.txt", m_match_directory, i, j);
    
    if (FileExists(buf)) {
        FILE *f = fopen(buf, "r");
	
        /* Read the number of matches */
        int num_matches;
        fscanf(f, "%d", &num_matches);
		
        if (num_matches < MIN_MATCHES) {
            fclose(f);
            return;
        }

        SetMatch(i, j);

        std::vector<KeypointMatch> matches;
		
        for (int k = 0; k < num_matches; k++) {
            int idx1, idx2;
            fscanf(f, "%d %d", &idx1, &idx2);
		   
#ifdef KEY_LIMIT
            if (idx1 > KEY_LIMIT || idx2 > KEY_LIMIT)
                continue;
#endif /* KEY_LIMIT */
 
            KeypointMatch m;

            m.m_idx1 = idx1;
            m.m_idx2 = idx2;
		    
            matches.push_back(m);
        }
        
        MatchIndex idx = GetMatchIndex(i, j);
        // m_match_lists[idx] = matches;
        m_matches.GetMatchList(idx) = matches;
    
        fclose(f);
    } else {
        // RemoveMatch(i, j);
        // RemoveMatch(j, i);
    }
}

void BaseApp::LoadMatchTable(const char *filename) {
    if (m_matches_loaded)
        return;  /* we already loaded the matches */

    printf("[LoadMatchTable] Loading matches\n");

    /* Set all matches to false */
    // ClearMatches();
    RemoveAllMatches();

    FILE *f = fopen(filename, "r");

    if (f == NULL) {
        printf("[LoadMatchTable] Error opening file %s for reading\n",
               filename);
        exit(1);
    }

    char buf[256];
    while (fgets(buf, 256, f)) {
        /* Read the images */
        int i1, i2;
        sscanf(buf, "%d %d\n", &i1, &i2);

        SetMatch(i1, i2);

        /* Read the number of matches */
        int nMatches;
        fscanf(f, "%d\n", &nMatches);

        /* Read the matches */
        std::vector<KeypointMatch> matches;
        for (int i = 0; i < nMatches; i++) {
            int k1, k2;
            fscanf(f, "%d %d\n", &k1, &k2);

#ifdef KEY_LIMIT
            if (k1 > KEY_LIMIT || k2 > KEY_LIMIT)
                continue;
#endif /* KEY_LIMIT */

            KeypointMatch m;

            m.m_idx1 = k1;
            m.m_idx2 = k2;
		    
            matches.push_back(m);
        }

        MatchIndex idx = GetMatchIndex(i1, i2);
        m_matches.GetMatchList(idx) = matches;
    }
    
    fclose(f);
}    

void BaseApp::LoadMatchIndexes(const char *index_dir)
{
    int num_images = GetNumImages();
    
    for (int i = 0; i < num_images; i++) {
        char buf[256];
        sprintf(buf, "%s/match-%03d.txt", index_dir, i);

        FILE *f = fopen(buf, "r");
        
        if (f == NULL)
            continue;
        
        printf("[LoadMatchIndexes] Loading matches for image %d... ", i);

        int num_matched_images = 0;
        int index;
        while(fgets(buf, 256, f) != NULL) {
            sscanf(buf, "%d\n", &index);

            int num_matches;
            fscanf(f, "%d\n", &num_matches);
            
            std::vector<KeypointMatch> matches;
		
            for (int k = 0; k < num_matches-1; k++) {//-1 is a temp hack
                int idx1, idx2;
                fscanf(f, "%d %d\n", &idx1, &idx2);

#ifdef KEY_LIMIT
                if (idx1 > KEY_LIMIT || idx2 > KEY_LIMIT)
                    continue;
#endif /* KEY_LIMIT */

                KeypointMatch m;

                m.m_idx1 = idx1;
                m.m_idx2 = idx2;
                
                matches.push_back(m);
            }

            if (num_matches < MIN_MATCHES || index >= num_images) {
                if (index >= num_images) 
                    printf("[LoadMatchIndexes] image index %d > num_images\n",
                           index);

                // RemoveMatch(i, index);
                matches.clear();
            } else { 
                SetMatch(i, index);
                MatchIndex idx = GetMatchIndex(i, index);
                // m_match_lists[idx] = matches;
                m_matches.GetMatchList(idx) = matches;

                num_matched_images++;
            }
        }

        printf("%d match files loaded.\n", num_matched_images);
        fflush(stdout);

        fclose(f);
    }
}

/* Load matches from files */
void BaseApp::LoadMatches() {
    if (m_matches_loaded)
	return;  /* we already loaded the matches */

    if (m_match_table != NULL) {
        LoadMatchTable(m_match_table);
    } else if (m_match_index_dir != NULL) {
        LoadMatchIndexes(m_match_index_dir);
    } else {
        printf("[LoadMatches] Loading matches\n");

        int num_images = GetNumImages();

        FILE *f = fopen("match-index.txt", "r");
        if (f == NULL) {
            /* Try all pairs */
            for (int i = 0; i < num_images; i++) {
                for (int j = i+1; j < num_images; j++) {
                    ReadMatchFile(i, j);
                }
            }
        } else {
            printf("[LoadMatches] Reading matches from 'match-index.txt'\n");
            fflush(stdout);

            /* Set all matches to false */
            RemoveAllMatches();

            char buf[256];
            unsigned int count = 0;
            while (fgets(buf, 256, f)) {
                int i1, i2;
                sscanf(buf, "%d %d\n", &i1, &i2);
                ReadMatchFile(i1, i2);
                count++;
            }

            printf("[LoadMatches] Read %d match files\n", count);
            fflush(stdout);

            fclose(f);
        }
    }

    PruneDoubleMatches();

    m_matches_loaded = true;
}

void BaseApp::RemoveAllMatches() 
{
    m_matches.RemoveAll();
}

/* Load a list of image names from a file */
void BaseApp::LoadImageNamesFromFile(FILE *f)
{
    m_image_data.clear();

    char buf[256];
    int idx = 0;

    while (fgets(buf, 256, f)) {
	ImageData data;
	data.InitFromString(buf, m_image_directory, m_fisheye);

        /* Try to find a keypoint file */

        if (strcmp(m_key_directory, ".") != 0) {
            char key_buf[256];
            data.GetBaseName(key_buf);

            char key_path[512];
            sprintf(key_path, "%s/%s.key", m_key_directory, key_buf);
            data.m_key_name = strdup(key_path);
        } else {
            /* FIXME: I think this causes a memory leak */
            char key_buf[256];
            strcpy(key_buf, data.m_name);
            int len = strlen(key_buf);
            key_buf[len - 3] = 'k';
            key_buf[len - 2] = 'e';
            key_buf[len - 1] = 'y';

            data.m_key_name = strdup(key_buf);
        }

#if 0
	if (log != NULL) {
	    log->AppendText("  ");
	    log->AppendText(buf);
	    log->AppendText("\n");
	}
#endif

#if 0
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
	wxStringTokenizer t(str, wxT(" "));

	while (t.HasMoreTokens()) {
	    wxString tok = t.GetNextToken();
	    toks.push_back(tok);
	}
#endif

#if 0
	if (log != NULL) {
	    log->AppendText("  ");
	    log->AppendText(buf);
	    log->AppendText("\n");
	}
#endif


	int num_toks = (int) toks.size();

	bool fisheye = m_fisheye; // false;
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

	// m_imgs.push_back(img);
	// printf("Adding image %s\n", toks[0].c_str());
        char buf[512];
        sprintf(buf, "%s/%s", m_image_directory, toks[0].c_str());

	data.m_name = strdup(buf);
	data.m_img = NULL;
	data.m_thumb = NULL;
	data.m_thumb8 = NULL;
	data.m_wximage = NULL;
	data.m_image_loaded = false;
	data.m_keys_loaded = false;

	data.m_fisheye = fisheye;
	data.m_has_init_focal = has_init_focal;
	data.m_init_focal = init_focal;
	data.m_camera.m_adjusted = false;
	data.m_texture_index = -1;
#endif

	m_image_data.push_back(data);

	idx++;

	// wxSafeYield();
    }

    /* Create the match table */
    m_matches = MatchTable(GetNumImages());

    RemoveAllMatches();

    m_matches_computed = true;
    m_num_original_images = GetNumImages();

    if (m_use_intrinsics)
        ReadIntrinsicsFile();
}

/* Read in information about the world */
void BaseApp::ReadBundleFile(const char *filename)
{
    printf("[ReadBundleFile] Reading file...\n");

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Error opening file %s for reading\n", filename);
        return;
    }

    int num_images, num_points;

    char first_line[256];
    fgets(first_line, 256, f);
    if (first_line[0] == '#') {
        double version;
        sscanf(first_line, "# Bundle file v%lf", &version);

        m_bundle_version = version;
        printf("[ReadBundleFile] Bundle version: %0.3f\n", version);

        fscanf(f, "%d %d\n", &num_images, &num_points);
    } else if (first_line[0] == 'v') {
        double version;
        sscanf(first_line, "v%lf", &version);
        m_bundle_version = version;
        printf("[ReadBundleFile] Bundle version: %0.3f\n", version);

        fscanf(f, "%d %d\n", &num_images, &num_points);
    } else {
        m_bundle_version = 0.1;
        sscanf(first_line, "%d %d\n", &num_images, &num_points);
    }

    printf("[ReadBundleFile] Reading %d images and %d points...\n",
        num_images, num_points);

    if (num_images != GetNumImages()) {
        printf("Error: number of images doesn't match file!\n");
        return;
    }

    /* Read cameras */
    for (int i = 0; i < num_images; i++) {
        double focal_length;
        double R[9];
        double t[3];
        double k[2] = { 0.0, 0.0 };

        if (m_bundle_version >= 0.4) {
            char name[512];
            int w, h;
            fscanf(f, "%s %d %d\n", name, &w, &h);
        }
        
        /* Focal length */
        if (m_bundle_version > 0.1) {
            fscanf(f, "%lf %lf %lf\n", &focal_length, k+0, k+1);
        } else {
            fscanf(f, "%lf\n", &focal_length);
        }

        /* Rotation */
        fscanf(f, "%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n", 
            R+0, R+1, R+2, R+3, R+4, R+5, R+6, R+7, R+8);
        /* Translation */
        fscanf(f, "%lf %lf %lf\n", t+0, t+1, t+2);

#if 0
        if (m_bundle_version < 0.3) {
            R[2] = -R[2];
            R[5] = -R[5];
            R[6] = -R[6];
            R[7] = -R[7];
            t[2] = -t[2];
        }
#endif

        if (focal_length <= 100.0 || m_image_data[i].m_ignore_in_bundle) {
            /* No (or bad) information about this camera */
            m_image_data[i].m_camera.m_adjusted = false;
        } else {
            CameraInfo cd;

            cd.m_adjusted = true;
            cd.m_width = m_image_data[i].GetWidth();
            cd.m_height = m_image_data[i].GetHeight();
            cd.m_focal = focal_length;
            cd.m_k[0] = k[0];
            cd.m_k[1] = k[1];
            memcpy(cd.m_R, R, sizeof(double) * 9);
            memcpy(cd.m_t, t, sizeof(double) * 3);

            cd.Finalize();

            m_image_data[i].m_camera = cd;
        }
    }

    /* Read points */
    m_point_data.clear();
    m_point_data.resize(num_points);

    int num_min_views_points = 0;
    for (int i = 0; i < num_points; i++) {
        PointData &pt = m_point_data[i];

        /* Position */
        fscanf(f, "%lf %lf %lf\n", 
            pt.m_pos + 0, pt.m_pos + 1, pt.m_pos + 2);

        // if (m_bundle_version < 0.3)
        //     pt.m_pos[2] = -pt.m_pos[2];

        /* Color */
        fscanf(f, "%f %f %f\n", 
            pt.m_color + 0, pt.m_color + 1, pt.m_color + 2);

        int num_visible;
        fscanf(f, "%d", &num_visible);
		pt.m_num_vis=num_visible;
        if (num_visible >=3)
            num_min_views_points++;

        // pt.m_views.resize(num_visible);
        for (int j = 0; j < num_visible; j++) {
            int view, key;
            fscanf(f, "%d %d", &view, &key);

            if (!m_image_data[view].m_camera.m_adjusted) {
                // printf("[ReadBundleFile] "
                //        "Removing view %d from point %d\n", view, i);
            } else {
                /* Check cheirality */
                bool val = (m_bundle_version >= 0.3);

                double proj_test[2];
                if (m_image_data[view].m_camera.
                    Project(pt.m_pos, proj_test) == val) {
                    
                    pt.m_views.push_back(ImageKey(view, key));
                } else {
                    printf("[ReadBundleFile] "
                           "Removing view %d from point %d [cheirality]\n", 
                               view, i);
                    // pt.m_views.push_back(ImageKey(view, key));
                }   
            }
            // pt.m_views.push_back(ImageKey(view, key));
            
            if (m_bundle_version >= 0.3) {
                double x, y;
                fscanf(f, "%lf %lf", &x, &y);
            }
        }

        // #define CROP_POINT_CLOUD
#ifdef CROP_POINT_CLOUD
        const double x_min = 1.327;
        const double x_max = 3.556;
        const double y_min = -1.414;
        const double y_max = 1.074;
        const double z_min = -5.502;
        const double z_max = -3.288;
        
        if (pt.m_pos[0] < x_min || pt.m_pos[0] > x_max ||
            pt.m_pos[1] < y_min || pt.m_pos[1] > y_max ||
            pt.m_pos[2] < z_min || pt.m_pos[2] > z_max) {
         
            pt.m_views.clear();
        }
#endif /* CROP_POINT_CLOUD */
    }

#if 0    
    /* Read outliers */
    int num_outliers;
    fscanf(f, "%d", &num_outliers);

    for (int i = 0; i < num_outliers; i++) {
        ImageKey ik;
        fscanf(f, "%d %d", &(ik.first), &(ik.second));
        m_outliers.push_back(ik);
    }
#endif

    fclose(f);

    printf("[ReadBundleFile] %d / %d points visible to more than 2 cameras!\n", 
        num_min_views_points, num_points);
}

void BaseApp::ReloadBundleFile(const char *filename)
{
#ifndef __DEMO__
    /* Count the old number of cameras */
    int num_images = GetNumImages();

    int old_num_cameras = 0;
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted)
            old_num_cameras++;
    }


    /* Save the previous model */
    std::vector<PointData> old_points = m_point_data;
    std::vector<ImageData> old_images = m_image_data;

    /* Load the new model */
    ClearModel();
    ReadBundleFile(filename);

    if (m_bundle_version < 0.3)
        FixReflectionBug();

    /* Count the new number of cameras */
    int num_cameras = 0;
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted)
            num_cameras++;
    }

    int old_num_points = old_points.size();
    int new_num_points = m_point_data.size();

    std::vector<v3_t> points_old_csp, points_new_csp;
    /* Find point correspondences */
    for (int i = 0; i < old_num_points; i++) {
        for (int j = i - 5; j < i + 5; j++) {
            if (j < 0 || j >= new_num_points) continue;

            float *old_col = old_points[i].m_color;
            float *col = m_point_data[j].m_color;

            if (old_col[0] == col[0] && 
                old_col[1] == col[1] && 
                old_col[2] == col[2]) {

                    double *old_pos = old_points[i].m_pos;
                    double *pos = m_point_data[i].m_pos;

                    points_old_csp.push_back(v3_new(old_pos[0], 
                        old_pos[1], 
                        old_pos[2]));

                    points_new_csp.push_back(v3_new(pos[0], pos[1], pos[2]));

                    goto Next;
            }
        }
Next: ;
    }

    int num_csp_points = points_old_csp.size();

    int num_points = old_num_cameras + num_csp_points;
    v3_t *left_points = new v3_t[num_points];
    v3_t *right_points = new v3_t[num_points];

    int count = 0;
    for (int i = 0; i < num_images; i++) {
        if (old_images[i].m_camera.m_adjusted) {
            double left_pos[3], right_pos[3];

            m_image_data[i].m_camera.GetPosition(left_pos);
            old_images[i].m_camera.GetPosition(right_pos);

            left_points[count] = 
                v3_new(left_pos[0], left_pos[1], left_pos[2]);
            right_points[count] = 
                v3_new(right_pos[0], right_pos[1], right_pos[2]);

            count++;
        }
    }

    for (int i = 0; i < num_csp_points; i++) {
        left_points[count] = points_new_csp[i];
        right_points[count] = points_old_csp[i];
        count++;	
    }

    /* Do the registration */
    double T[16];
    align_horn_3D(num_points, right_points, left_points, 1, T);

    /* Transform the world */
    memcpy(m_xform, T, 16 * sizeof(double));
    TransformWorldReal();
    // TransformWorld();

    delete [] left_points;
    delete [] right_points;
#endif /* __DEMO__ */
}

/* Clear the current model */
void BaseApp::ClearModel()
{
    int num_images = GetNumImages();

    for (int i = 0; i < num_images; i++) {
        m_image_data[i].m_camera.m_adjusted = false;
    }

    m_point_data.clear();
}


#ifndef __DEMO__
/* Dump an output file containing information about the current
* state of the world */
void BaseApp::DumpOutputFile(const char *output_dir, const char *filename, 
                             int num_images, int num_cameras, int num_points,
                             int *added_order, 
                             camera_params_t *cameras, 
                             v3_t *points, v3_t *colors,
                             std::vector<ImageKeyVector> &pt_views
                             /*bool output_radial_distortion*/)
{
    clock_t start = clock();

    int num_visible_points = 0;
    
    for (int i = 0; i < num_points; i++) {
        if (pt_views[i].size() > 0)
            num_visible_points++;
    }

    char buf[256];
    sprintf(buf, "%s/%s", output_dir, filename);

    FILE *f = fopen(buf, "w");
    if (f == NULL) {
        printf("Error opening file %s for writing\n", buf);
        return;
    }

    /* Print version number */
    fprintf(f, "# Bundle file v0.3\n");
    /* Print number of cameras and points */
    fprintf(f, "%d %d\n", num_images, num_visible_points);

    /* Dump cameras */
    for (int i = 0; i < num_images; i++) {

#if 0
        /* Print the name of the file */
        fprintf(f, "%s %d %d\n", 
                m_image_data[i].m_name, 
                m_image_data[i].GetWidth(), m_image_data[i].GetHeight());
#endif

        int idx = -1;
        for (int j = 0; j < num_cameras; j++) {
            if (added_order[j] == i) {
                idx = j;
                break;
            }
        }

        if (idx == -1) {
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n0 0 0\n0 0 0\n0 0 0\n");
        } else {
            fprintf(f, "%0.10e %0.10e %0.10e\n", 
                    cameras[idx].f, cameras[idx].k[0], cameras[idx].k[1]);

            fprintf(f, "%0.10e %0.10e %0.10e\n", 
                cameras[idx].R[0], 
                cameras[idx].R[1], 
                cameras[idx].R[2]);
            fprintf(f, "%0.10e %0.10e %0.10e\n", 
                cameras[idx].R[3], 
                cameras[idx].R[4], 
                cameras[idx].R[5]);
            fprintf(f, "%0.10e %0.10e %0.10e\n", 
                cameras[idx].R[6], 
                cameras[idx].R[7], 
                cameras[idx].R[8]);

            double t[3];
            matrix_product(3, 3, 3, 1, cameras[idx].R, cameras[idx].t, t);
            matrix_scale(3, 1, t, -1.0, t);
            fprintf(f, "%0.10e %0.10e %0.10e\n", t[0], t[1], t[2]);
        }
    }

    /* Dump points */
    for (int i = 0; i < num_points; i++) {
        int num_visible = (int) pt_views[i].size();

        if (num_visible > 0) {

            /* Position */
            fprintf(f, "%0.10e %0.10e %0.10e\n", 
                    Vx(points[i]), Vy(points[i]), Vz(points[i]));

            /* Color */
            fprintf(f, "%d %d %d\n", 
                    iround(Vx(colors[i])), 
                    iround(Vy(colors[i])), 
                    iround(Vz(colors[i])));

            int num_visible = (int) pt_views[i].size();
            fprintf(f, "%d", num_visible);
            for (int j = 0; j < num_visible; j++) {
                int img = added_order[pt_views[i][j].first];
                int key = pt_views[i][j].second;
                
                double x = m_image_data[img].m_keys[key].m_x;
                double y = m_image_data[img].m_keys[key].m_y;
                
                fprintf(f, " %d %d %0.4f %0.4f", img, key, x, y);
            }
            
            fprintf(f, "\n");
        }
    }

#if 0
    /* Finally, dump all outliers */
    ImageKeyVector outliers;
    for (int i = 0; i < num_images; i++) {
        /* Find the index of this camera in the ordering */
        int idx = -1;
        for (int j = 0; j < num_cameras; j++) {
            if (added_order[j] == i) {
                idx = j;
                break;
            }
        }

        if (idx == -1) continue;

        int num_keys = GetNumKeys(i);
        for (int j = 0; j < num_keys; j++) {
            if (GetKey(i,j).m_extra == -2) {
                outliers.push_back(ImageKey(i,j));
            }
        }
    }

    int num_outliers = (int) outliers.size();
    fprintf(f, "%d\n", num_outliers);

    for (int i = 0; i < num_outliers; i++) {
        fprintf(f, "%d %d\n", outliers[i].first, outliers[i].second);
    }
#endif

    fclose(f);

    clock_t end = clock();

    printf("[DumpOutputFile] Wrote file in %0.3fs\n",
        (double) (end - start) / (double) CLOCKS_PER_SEC);
}
#endif


/* Write XML files */

/* Camera I/O */
void BaseApp::WriteCamerasXML(const char *filename)
{
    FILE *f = fopen(filename, "w");

    if (f == NULL) {
	printf("[WriteCamerasXML] "
	       "Error opening file %s for writing\n", filename);
	return;
    }
    
    fprintf(f, "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\n");
    const char *url_base = 
	"http://grail.cs.cornell.edu/projects/phototour/trevi/images";

    fprintf(f, "<url_base> %s </url_base>\n", url_base);
    fprintf(f, "<cameras>\n");
    
    int num_images = (int) m_image_data.size();
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted /*&& m_image_data[i].m_licensed*/) {
            m_image_data[i].WriteCameraXML(f);
        }
    }

    fprintf(f, "</cameras>\n");
    fclose(f);
}



/* Point I/O */
void BaseApp::WritePointsXML(const char *filename) 
{
    FILE *f = fopen(filename, "w");
    int min_views = 3;

    if (f == NULL) {
	printf("[WritePointsXML] "
	       "Error opening file %s for writing\n", filename);
	return;
    }
    
    fprintf(f, "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\n");
    fprintf(f, "<points>\n");
    
    int num_points = (int) m_point_data.size();
    int num_ge2 = 0;

    for (int i = 0; i < num_points; i++) {
        if (m_num_views_orig[i] >= min_views) {
            m_point_data[i].WriteXML(f);
	    num_ge2++;
        }
    }

    fprintf(f, "</points>\n");
    fclose(f);

    printf("[WritePointsXML] %d / %d points seen by >= %d views\n",
	   num_ge2, num_points, min_views);
}

/* Point I/O */
void BaseApp::WritePointsGeoXML(const char *filename) 
{
    FILE *f = fopen(filename, "w");
    int min_views = 2;

    if (f == NULL) {
	printf("[WritePointsXML] "
	       "Error opening file %s for writing\n", filename);
	return;
    }
    
    fprintf(f, "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\n");
    fprintf(f, "<points>\n");
    
    int num_points = (int) m_point_data.size();
    int num_ge2 = 0;

    for (int i = 0; i < num_points; i++) {
        if ((int) m_point_data[i].m_views.size() >= min_views) {
            m_point_data[i].WriteGeoXML(f);
            num_ge2++;
        }
    }

    fprintf(f, "</points>\n");
    fclose(f);

    printf("[WritePointsXML] %d / %d points seen by >= %d views\n",
	   num_ge2, num_points, min_views);
}

void BaseApp::ReadMatchTable(const char *append) 
{
    int num_images = GetNumImages();
    unsigned long int num_matches_total = 0;

    char buf[256];
    sprintf(buf, "nmatches%s.txt", append);
    FILE *f0 = fopen(buf, "r");

    sprintf(buf, "matches%s.txt", append);
    FILE *f1 = fopen(buf, "r");
    
    if (f0 == NULL || f1 == NULL) {
        printf("[ReadMatchTable] "
               "Error opening files for reading.\n");
        return;
    }

    int num_images_check;
    fscanf(f0, "%d\n", &num_images_check);

    assert(num_images == num_images_check);

    RemoveAllMatches();

    for (int i = 0; i < num_images; i++) {
        for (int j = 0; j < num_images; j++) {
            MatchIndex idx = GetMatchIndex(i, j);

            int num_matches;
            fscanf(f0, "%d", &num_matches);

            if (num_matches > 0) {
                // m_match_lists[idx].clear();
                SetMatch(i, j);
                std::vector<KeypointMatch> &list = m_matches.GetMatchList(idx);
                // m_matches.ClearMatch(idx);

                for (int k = 0; k < num_matches; k++) {
                    KeypointMatch m;

                    int idx1, idx2;
                    fscanf(f1, "%d %d", &(idx1), &(idx2));

#ifdef KEY_LIMIT
                    if (idx1 > KEY_LIMIT || idx2 > KEY_LIMIT)
                        continue;
#endif /* KEY_LIMIT */

                    m.m_idx1 = idx1;
                    m.m_idx2 = idx2;

                    // m_match_lists[idx].push_back(m);
                    list.push_back(m);
                }

                num_matches_total += num_matches;
            }
        }
    }

    printf("[ReadMatchTable] Read %lu matches in total\n",
           num_matches_total);

    fclose(f0);
    fclose(f1);
}

void BaseApp::WriteMatchTable(const char *append) 
{
    int num_images = GetNumImages();

    char buf[256];
    sprintf(buf, "nmatches%s.txt", append);
    FILE *f0 = fopen(buf, "w");

    sprintf(buf, "matches%s.txt", append);
    FILE *f1 = fopen(buf, "w");
    
    if (f0 == NULL || f1 == NULL) {
        printf("[WriteMatchTable] "
               "Error opening files for writing.\n");
        return;
    }

    fprintf(f0, "%d\n", num_images);

    for (int i = 0; i < num_images; i++) {
        for (int j = 0; j < num_images; j++) {
            if (i >= j) {
                fprintf(f0, "0 ");
                fprintf(f1, "\n");
            } else {
                if (ImagesMatch(i, j)) {
                    MatchIndex idx = GetMatchIndex(i, j);
                    std::vector<KeypointMatch> &list = 
                        m_matches.GetMatchList(idx);

                    unsigned int num_matches = list.size();

                    fprintf(f0, "%d ", num_matches);
                    
                    for (unsigned int k = 0; k < num_matches; k++) {
                        KeypointMatch m = list[k];
                        fprintf(f1, "%d %d ", m.m_idx1, m.m_idx2);    
                    }
                    fprintf(f1, "\n");
                } else {
                    fprintf(f0, "0 ");
                }
            }
        }

        fprintf(f0, "\n");
    }

    fclose(f0);
    fclose(f1);
}

#ifndef __DEMO__
static char ply_header[] = 
"ply\n"
"format ascii 1.0\n"
"element face 0\n"
"property list uchar int vertex_indices\n"
"element vertex %d\n"
"property float x\n"
"property float y\n"
"property float z\n"
"property uchar diffuse_red\n"
"property uchar diffuse_green\n"
"property uchar diffuse_blue\n"
"end_header\n";

/* Write point files to a ply file */
void BaseApp::DumpPointsToPly(const char *output_directory, 
                              const char *filename, 
                              int num_points, int num_cameras, 
                              v3_t *points, v3_t *colors,
                              camera_params_t *cameras) 
{
    int num_good_pts = 0;

    for (int i = 0; i < num_points; i++) {
	if (Vx(colors[i]) == 0x0 && 
	    Vy(colors[i]) == 0x0 && 
	    Vz(colors[i]) == 0xff) 
	    continue;
	num_good_pts++;
    }
    
    char ply_out[256];
    sprintf(ply_out, "%s/%s", output_directory, filename);

    FILE *f = fopen(ply_out, "w");

    if (f == NULL) {
	printf("Error opening file %s for writing\n", ply_out);
	return;
    }

    /* Print the ply header */
    fprintf(f, ply_header, num_good_pts + 2 * num_cameras);

    /* Now triangulate all the correspondences */
    for (int i = 0; i < num_points; i++) {
	if (Vx(colors[i]) == 0x0 && 
	    Vy(colors[i]) == 0x0 && 
	    Vz(colors[i]) == 0xff) 
	    continue;

	/* Output the vertex */
	fprintf(f, "%0.6e %0.6e %0.6e %d %d %d\n", 
		Vx(points[i]), Vy(points[i]), Vz(points[i]),
		iround(Vx(colors[i])), 
		iround(Vy(colors[i])), 
		iround(Vz(colors[i])));
    }

    for (int i = 0; i < num_cameras; i++) {
	double c[3];

	double Rinv[9];
	matrix_invert(3, cameras[i].R, Rinv);

        memcpy(c, cameras[i].t, 3 * sizeof(double));
	
	if ((i % 2) == 0)
	    fprintf(f, "%0.6e %0.6e %0.6e 0 255 0\n", c[0], c[1], c[2]);
	else
	    fprintf(f, "%0.6e %0.6e %0.6e 255 0 0\n", c[0], c[1], c[2]);

	double p_cam[3] = { 0.0, 0.0, -0.05 };
	double p[3];

	matrix_product(3, 3, 3, 1, Rinv, p_cam, p);

	p[0] += c[0];
	p[1] += c[1];
	p[2] += c[2];

	fprintf(f, "%0.6e %0.6e %0.6e 255 255 0\n",
		p[0], p[1], p[2]);
    }

    fclose(f);
}
#endif

/* Grab the color of each keypoint */
void BaseApp::ReadKeyColors() 
{
    int num_images = GetNumImages();

    for (int i = 0; i < num_images; i++) {
	m_image_data[i].ReadKeyColors();
    }
}

/* Read camera constraints */
void BaseApp::ReadCameraConstraints() 
{
    if (FileExists("camera-constraints.txt")) {
	printf("[ReadCameraConstraints] Reading constraints\n");

	FILE *f = fopen("camera-constraints.txt", "r");
	char buf[256];

	while (fgets(buf, 256, f) != NULL) {
	    if (isspace(buf[0]) || buf[0] == '%')
		continue;  /* comment or whitespace */

	    int cam_idx;
	    double x, y, z, xw, yw, zw;
	    sscanf(buf, "%d %lf %lf %lf %lf %lf %lf", 
		   &cam_idx, &x, &y, &z, &xw, &yw, &zw);

	    printf("  Constraints on camera %d: %0.3f, %0.3f, %0.3f\n"
		   "    (weights %0.3f, %0.3f, %0.3f)\n", 
		   cam_idx, x, y, z, xw, yw, zw);

	    if (x != -999.0) {
		m_image_data[cam_idx].m_camera.m_constrained[0] = true;
		m_image_data[cam_idx].m_camera.m_constraints[0] = x;
		m_image_data[cam_idx].m_camera.m_constraint_weights[0] = xw;
	    }

	    if (y != -999.0) {
		m_image_data[cam_idx].m_camera.m_constrained[1] = true;
		m_image_data[cam_idx].m_camera.m_constraints[1] = y;
		m_image_data[cam_idx].m_camera.m_constraint_weights[1] = yw;
	    }
	    
	    if (z != -999.0) {
		m_image_data[cam_idx].m_camera.m_constrained[2] = true;
		m_image_data[cam_idx].m_camera.m_constraints[2] = z;
		m_image_data[cam_idx].m_camera.m_constraint_weights[2] = zw;
	    }
	}

	fclose(f);
    }
}

void BaseApp::ReadPointConstraints()
{
    FILE *f = fopen(m_point_constraint_file, "r");
    
    if (f == NULL) {
	printf("[ReadPointConstraints] Error opening file %s "
	       "for reading\n", m_point_constraint_file);
	return;
    }

    int num_points = (int) m_point_data.size();
    m_point_constraints = new v3_t[num_points];

    for (int i = 0; i < num_points; i++) {
	m_point_constraints[i] = v3_new(0.0, 0.0, 0.0);
    }

    char buf[256];
    while (fgets(buf, 256, f) != NULL) {
	double x0, y0, z0;
	double x, y, z;
	sscanf(buf, "%lf %lf %lf %lf %lf %lf", &x0, &y0, &z0, &x, &y, &z);

	int pt_idx = -1;
	double min_dist = DBL_MAX;
	for (int i = 0; i < num_points; i++) {
	    double dx = m_point_data[i].m_pos[0] - x0;
	    double dy = m_point_data[i].m_pos[1] - y0;
	    double dz = m_point_data[i].m_pos[2] - z0;

	    double dsq = dx * dx + dy * dy + dz * dz;
	    
	    if (dsq < min_dist) {
		pt_idx = i;
		min_dist = dsq;
	    }
	}

	m_point_constraints[pt_idx] = v3_new(x, y, -z);

	printf("[ReadPointConstraints] Constraining %d: "
	       "%0.3f %0.3f %0.3f (%0.3f %0.3f %0.3f) => %0.3f %0.3f %0.3f\n",
	       pt_idx, 
	       m_point_data[pt_idx].m_pos[0], 
	       m_point_data[pt_idx].m_pos[1], 
	       m_point_data[pt_idx].m_pos[2], 
	       x0, y0, z0, x, y, z);
    }
}

typedef struct {
    double K[9];
    double k[5];
} intrinsics_t;

/* Read intrinsics */
void BaseApp::ReadIntrinsicsFile()
{
    printf("[ReadIntrinsicsFile] Reading intrinsics...\n");

    assert(m_intrinsics_file != NULL);
    
    FILE *f = fopen(m_intrinsics_file, "r");
    assert(f != NULL);
    
    int num_intrinsics = 0;
    fscanf(f, "%d\n", &num_intrinsics);

    std::vector<intrinsics_t> Ks;
    for (int i = 0; i < num_intrinsics; i++) {
        intrinsics_t I;

        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
               I.K + 0, I.K + 1, I.K + 2, 
               I.K + 3, I.K + 4, I.K + 5, 
               I.K + 6, I.K + 7, I.K + 8);

        fscanf(f, "%lf %lf %lf %lf %lf\n",
               I.k + 0, I.k + 1, I.k + 2, I.k + 3, I.k + 4);

        Ks.push_back(I);
    }

    int num_images = GetNumImages();
    for (int i = 0; i < num_images; i++) {
        assert(m_image_data[i].m_has_init_focal);

        double f = m_image_data[i].m_init_focal;
        
        double min_dist = DBL_MAX;
        int best_K = -1;
        for (int j = 0; j < num_intrinsics; j++) {
            double f_j = 0.5 * (Ks[j].K[0] + Ks[j].K[4]);
            double dist = fabs(f_j - f);
            if (dist < min_dist) {
                best_K = j;
                min_dist = dist;
            }
        }

        printf("  image %d has intrinsics %d\n", i, best_K);

        memcpy(m_image_data[i].m_K, Ks[best_K].K, 9 * sizeof(double));
        memcpy(m_image_data[i].m_k, Ks[best_K].k, 5 * sizeof(double));
        m_image_data[i].m_known_intrinsics = true;
    }

    fclose(f);
}

/* Read the ignore file */
void BaseApp::ReadIgnoreFile()
{
    if (m_ignore_file == NULL)
	return;

    FILE *f = fopen(m_ignore_file, "r");

    if (f == NULL) {
	printf("[ReadIgnoreFile] Error opening file %s "
	       "for reading\n", m_ignore_file);
	return;
    }
    
    char buf[256];
    int num_images = GetNumImages();
    while (fgets(buf, 255, f)) {
	int img = atoi(buf);
	
	if (img < 0 || img >= num_images) {
	    printf("[ReadIgnoreFile] "
		   "Error: image %d out of range\n", img);
	    continue;
	}
	
	printf("[ReadIgnoreFile] Ignoring image %d\n", img);
	m_image_data[img].m_ignore_in_bundle = true;
    }

    fclose(f);

    fflush(stdout);
}

/* Initialize images read from a file without performing bundle
 * adjustment */
void BaseApp::InitializeImagesFromFile(FILE *f) 
{
    char buf[256];
    int num_images = GetNumImages();

    while (fgets(buf, 256, f)) {
	ImageData data;
	data.InitFromString(buf, m_image_directory, false);
	data.m_licensed = true;

	printf("[InitializeImagesFromFile] Initializing image %s\n",
	       data.m_name);

	/* Read the extra data */
        int img_idx = (int) m_image_data.size();
	if (data.ReadCamera() && data.ReadTracks(img_idx, m_point_data)) {
	    data.ReadMetadata();
            data.m_added = true;
	    m_image_data.push_back(data);
	}
    }

    UnscaleCameras(num_images);
}

void BaseApp::ReadLines3D(const char *filename) 
{
    /* Not implemented... */
}

void BaseApp::WriteLines3D(const char *filename)
{
    /* Not implemented... */
}
