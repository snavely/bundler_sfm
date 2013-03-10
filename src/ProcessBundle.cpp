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

/* ProcessBundle.cpp */
/* Perform operations on bundle files */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "BundlerApp.h"

#include "defines.h"
#include "matrix.h"
#include "util.h"

void BundlerApp::RotateCameras(char *rotate_file)
{
    int num_images = GetNumImages();
    FILE *ff = fopen(rotate_file, "r");
    if (ff == NULL) {
	printf("[SifterApp::RotateCameras] Error opening file %s "
	       "for reading\n", rotate_file);
	return;
    }
    
    for (int i = 0; i < num_images; i++) {
        char buf[512];
        double rot;
        fscanf(ff, "%s %lf\n", buf, &rot);

        if (rot != 0.0) {
            double rad = DEG2RAD(rot);
            double R[9] = { cos(rad), -sin(rad), 0.0,
                            sin(rad),  cos(rad), 0.0,
                            0.0, 0.0, 1.0 };

            double Rtmp[9];
            double ttmp[3];
        
            matrix_product(3, 3, 3, 3, R, m_image_data[i].m_camera.m_R, Rtmp);
            matrix_product(3, 3, 3, 1, R, m_image_data[i].m_camera.m_t, ttmp);
        
            memcpy(m_image_data[i].m_camera.m_R, Rtmp, 9 * sizeof(double));
            memcpy(m_image_data[i].m_camera.m_t, ttmp, 3 * sizeof(double));
        }
    }

    /* Output the new bundle.out file */
    char buf[256];
    sprintf(buf, "bundle.rotated.out");
    FILE *f = fopen(buf, "w");
    if (f == NULL) {
	printf("[SifterApp::RotateCameras] Error opening file %s "
	       "for writing\n", buf);
	return;
    }
    
    int num_points = (int) m_point_data.size();
    fprintf(f, "%d %d\n", num_images, num_points);

    /* Dump cameras */
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted) {
            fprintf(f, "0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            
	    continue;
        }
        
	fprintf(f, "%0.9e\n", m_image_data[i].m_camera.m_focal);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[0], 
		m_image_data[i].m_camera.m_R[1], 
		m_image_data[i].m_camera.m_R[2]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[3], 
		m_image_data[i].m_camera.m_R[4], 
		m_image_data[i].m_camera.m_R[5]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[6], 
		m_image_data[i].m_camera.m_R[7], 
		m_image_data[i].m_camera.m_R[8]);

	double *t = m_image_data[i].m_camera.m_t;
	fprintf(f, "%0.9e %0.9e %0.9e\n", t[0], t[1], t[2]);
    }
    
    /* Dump points */
    for (int i = 0; i < num_points; i++) {
	PointData &p = m_point_data[i];

	/* Position */
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		p.m_pos[0], p.m_pos[1], p.m_pos[2]);

	/* Color */
	fprintf(f, "%0.5e %0.5e %0.5e\n", 
		p.m_color[0], p.m_color[1], p.m_color[2]);

	int num_visible = (int) p.m_views.size();
	fprintf(f, "%d", num_visible);
	for (int j = 0; j < num_visible; j++) {
	    int view = p.m_views[j].first;
            // int key = p.m_views[j].second;

#if 0
            if (!m_image_data[view].m_keys_loaded)
                m_image_data[view].LoadKeys(false);

            double x = m_image_data[view].m_keys[key].m_x;
            double y = m_image_data[view].m_keys[key].m_y;

	    fprintf(f, " %d %d %.10f %.10f", 
                    view_new, p.m_views[j].second, x, y);
#else
	    fprintf(f, " %d %d", view, p.m_views[j].second);
#endif
	}

	fprintf(f, "\n");
    }
    
    fclose(f);    
}


void BundlerApp::ScaleFocalLengths(char *focal_file)
{
    int num_images = GetNumImages();
    FILE *ff = fopen(focal_file, "r");
    if (ff == NULL) {
	printf("[SifterApp::ScaleFocal] Error opening file %s "
	       "for reading\n", focal_file);
	return;
    }
    
    for (int i = 0; i < num_images; i++) {
        char buf[512];
        double scale;
        fscanf(ff, "%s %lf\n", buf, &scale);
        if (m_image_data[i].m_camera.m_adjusted) {
            m_image_data[i].m_camera.m_focal *= scale;
        }
    }

    /* Output the new bundle.out file */
    char buf[256];
    sprintf(buf, "bundle.scale.out");
    FILE *f = fopen(buf, "w");
    if (f == NULL) {
	printf("[SifterApp::ScaleFocal] Error opening file %s "
	       "for writing\n", buf);
	return;
    }
    
    int num_points = (int) m_point_data.size();
    fprintf(f, "%d %d\n", num_images, num_points);

    /* Dump cameras */
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted) {
            fprintf(f, "0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            
	    continue;
        }
        
	fprintf(f, "%0.9e\n", m_image_data[i].m_camera.m_focal);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[0], 
		m_image_data[i].m_camera.m_R[1], 
		m_image_data[i].m_camera.m_R[2]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[3], 
		m_image_data[i].m_camera.m_R[4], 
		m_image_data[i].m_camera.m_R[5]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[6], 
		m_image_data[i].m_camera.m_R[7], 
		m_image_data[i].m_camera.m_R[8]);

	double *t = m_image_data[i].m_camera.m_t;
	fprintf(f, "%0.9e %0.9e %0.9e\n", t[0], t[1], t[2]);
    }
    
    /* Dump points */
    for (int i = 0; i < num_points; i++) {
	PointData &p = m_point_data[i];

	/* Position */
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		p.m_pos[0], p.m_pos[1], p.m_pos[2]);

	/* Color */
	fprintf(f, "%0.5e %0.5e %0.5e\n", 
		p.m_color[0], p.m_color[1], p.m_color[2]);

	int num_visible = (int) p.m_views.size();
	fprintf(f, "%d", num_visible);
	for (int j = 0; j < num_visible; j++) {
	    int view = p.m_views[j].first;
            // int key = p.m_views[j].second;

#if 0
            if (!m_image_data[view].m_keys_loaded)
                m_image_data[view].LoadKeys(false);

            double x = m_image_data[view].m_keys[key].m_x;
            double y = m_image_data[view].m_keys[key].m_y;

	    fprintf(f, " %d %d %.10f %.10f", 
                    view_new, p.m_views[j].second, x, y);
#else
	    fprintf(f, " %d %d", view, p.m_views[j].second);
#endif
	}

	fprintf(f, "\n");
    }
    
    fclose(f);    
}

void BundlerApp::ScaleFocalLengths(double focal) 
{
    int num_images = GetNumImages();
    
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_camera.m_adjusted) {
            m_image_data[i].m_camera.m_focal *= focal;
        }
    }

    /* Output the new bundle.out file */
    char buf[256];
    sprintf(buf, "bundle.scale.out");
    FILE *f = fopen(buf, "w");
    if (f == NULL) {
	printf("[SifterApp::ScaleFocal] Error opening file %s "
	       "for writing\n", buf);
	return;
    }
    
    int num_points = (int) m_point_data.size();
    fprintf(f, "%d %d\n", num_images, num_points);

    /* Dump cameras */
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted) {
            fprintf(f, "0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            fprintf(f, "0 0 0\n");
            
	    continue;
        }
        
	fprintf(f, "%0.9e\n", m_image_data[i].m_camera.m_focal);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[0], 
		m_image_data[i].m_camera.m_R[1], 
		m_image_data[i].m_camera.m_R[2]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[3], 
		m_image_data[i].m_camera.m_R[4], 
		m_image_data[i].m_camera.m_R[5]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[6], 
		m_image_data[i].m_camera.m_R[7], 
		m_image_data[i].m_camera.m_R[8]);

	double *t = m_image_data[i].m_camera.m_t;
	fprintf(f, "%0.9e %0.9e %0.9e\n", t[0], t[1], t[2]);
    }
    
    /* Dump points */
    for (int i = 0; i < num_points; i++) {
	PointData &p = m_point_data[i];

	/* Position */
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		p.m_pos[0], p.m_pos[1], p.m_pos[2]);

	/* Color */
	fprintf(f, "%0.5e %0.5e %0.5e\n", 
		p.m_color[0], p.m_color[1], p.m_color[2]);

	int num_visible = (int) p.m_views.size();
	fprintf(f, "%d", num_visible);
	for (int j = 0; j < num_visible; j++) {
	    int view = p.m_views[j].first;
            // int key = p.m_views[j].second;

#if 0
            if (!m_image_data[view].m_keys_loaded)
                m_image_data[view].LoadKeys(false);

            double x = m_image_data[view].m_keys[key].m_x;
            double y = m_image_data[view].m_keys[key].m_y;

	    fprintf(f, " %d %d %.10f %.10f", 
                    view_new, p.m_views[j].second, x, y);
#else
	    fprintf(f, " %d %d", view, p.m_views[j].second);
#endif
	}

	fprintf(f, "\n");
    }
    
    fclose(f);    
}

void BundlerApp::OutputCompressed(const char *ext) 
{
    // EstimatePointNormals();

    /* Output the new list of images */
    char buf[256];
    sprintf(buf, "list.%s.txt", ext);

    FILE *f = fopen(buf, "w");
    
    if (f == NULL) {
	printf("[SifterApp::OutputCompress] "
               "Error opening file %s for writing\n", buf);
	
	return;
    }

    int num_images = GetNumImages();
    int num_adj_images = 0;
    int *map = new int[num_images];

    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;

        if (m_image_data[i].m_ignore_in_bundle) {
            printf("[OutputCompressed] Ignoring image %d\n", i);
            continue;
        }

        if (m_image_data[i].m_has_init_focal)
            fprintf(f, "%s %d %0.5f\n", m_image_data[i].m_name,
                    m_image_data[i].m_fisheye, m_image_data[i].m_init_focal);
        else
            fprintf(f, "%s\n", m_image_data[i].m_name);

	map[i] = num_adj_images;
	num_adj_images++;
    }

    fclose(f);
    
    /* Output the new bundle.out file */
    sprintf(buf, "bundle.%s.out", ext);
    f = fopen(buf, "w");
    if (f == NULL) {
	printf("[SifterApp::OutputCompress] Error opening file %s "
	       "for writing\n", buf);
	return;
    }

#if 0
    if (m_estimate_distortion && m_bundle_version > 0.1) {
        fprintf(f, "v%lf\n", m_bundle_version);
    }
#endif

    fprintf(f, "# Bundle file v0.3\n");
    
    int num_points = (int) m_point_data.size();
    int num_good_points = 0;
    
    for (int i = 0; i < num_points; i++) {
        if (m_point_data[i].m_views.size() >= 2)
            num_good_points++;
    }

    fprintf(f, "%d %d\n", num_adj_images, num_good_points);

    /* Dump cameras */
    for (int i = 0; i < num_images; i++) {
	if (!m_image_data[i].m_camera.m_adjusted)
	    continue;

        // if (m_estimate_distortion && m_bundle_version > 0.1) {

        fprintf(f, "%0.9e %0.9e %0.9e\n", 
                m_image_data[i].m_camera.m_focal, 
                m_image_data[i].m_camera.m_k[0],
                m_image_data[i].m_camera.m_k[1]);
        // } else {
        //     fprintf(f, "%0.9e\n", m_image_data[i].m_camera.m_focal);
        // }

	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[0], 
		m_image_data[i].m_camera.m_R[1], 
		m_image_data[i].m_camera.m_R[2]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[3], 
		m_image_data[i].m_camera.m_R[4], 
		m_image_data[i].m_camera.m_R[5]);
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		m_image_data[i].m_camera.m_R[6], 
		m_image_data[i].m_camera.m_R[7], 
		m_image_data[i].m_camera.m_R[8]);

	double *t = m_image_data[i].m_camera.m_t;
	fprintf(f, "%0.9e %0.9e %0.9e\n", t[0], t[1], t[2]);
    }
    
    /* Dump points */
    for (int i = 0; i < num_points; i++) {
	PointData &p = m_point_data[i];

	int num_visible = (int) p.m_views.size();
        if (num_visible < 2)
            continue;

	/* Position */
	fprintf(f, "%0.9e %0.9e %0.9e\n", 
		p.m_pos[0], p.m_pos[1], p.m_pos[2]);

	/* Color */
	fprintf(f, "%d %d %d\n", 
		iround(p.m_color[0]), 
                iround(p.m_color[1]), 
                iround(p.m_color[2]));

        /* Normal */
        // fprintf(f, "%0.5e %0.5e %0.5e\n",
        //         p.m_norm[0], p.m_norm[1], p.m_norm[2]);

	fprintf(f, "%d", num_visible);
	for (int j = 0; j < num_visible; j++) {
	    int view = p.m_views[j].first;
	    int view_new = map[view];
            int key = p.m_views[j].second;

#if 1
            if (!m_image_data[view].m_keys_loaded)
                m_image_data[view].LoadKeys(false);

            double x = 0.0, y = 0.0;

            if (key < (int) m_image_data[view].m_keys.size()) {
                x = m_image_data[view].m_keys[key].m_x;
                y = m_image_data[view].m_keys[key].m_y;
            }       

	    fprintf(f, " %d %d %0.4f %0.4f", view_new, key, x, y);
#else
	    fprintf(f, " %d %d", view_new, key);
#endif
	}

	fprintf(f, "\n");
    }
    
    for (int i = 0; i < num_images; i++) {
        if (m_image_data[i].m_keys_loaded)
            m_image_data[i].UnloadKeys();
    }

    fclose(f);

    delete [] map;
}

void BundlerApp::PruneBadPoints()
{
    int num_points = (int) m_point_data.size();
    const double MIN_ANGLE_THRESHOLD = 1.5;
    int num_pruned = 0;

    for (int i = 0; i < num_points; i++) {
        double *pos = m_point_data[i].m_pos;
        int num_views = (int) m_point_data[i].m_views.size();

        double max_angle = 0.0;
        for (int j = 0; j < num_views; j++) {
            int v1 = m_point_data[i].m_views[j].first;

            double p1[3], r1[3];
            m_image_data[v1].m_camera.GetPosition(p1);
            matrix_diff(3, 1, 3, 1, pos, p1, r1);
            double norm = matrix_norm(3, 1, r1);
            matrix_scale(3, 1, r1, 1.0 / norm, r1);

            for (int k = j+1; k < num_views; k++) {
                int v2 = m_point_data[i].m_views[k].first;

                double p2[3], r2[3];
                m_image_data[v2].m_camera.GetPosition(p2);
                matrix_diff(3, 1, 3, 1, pos, p2, r2);
                double norm = matrix_norm(3, 1, r2);
                matrix_scale(3, 1, r2, 1.0 / norm, r2);
                
                double dot;
                matrix_product(1, 3, 3, 1, r1, r2, &dot);

                double angle = 
                    acos(CLAMP(dot, -1.0 + 1.0e-8, 1.0 - 1.0e-8));

                if (angle > max_angle) {
                    max_angle = angle;
                }
            }
        }

        if (num_views < 3 || RAD2DEG(max_angle) < MIN_ANGLE_THRESHOLD) {
            printf("[PruneBadPoints] Removing point %d with angle %0.3f\n",
                   i, RAD2DEG(max_angle));

            m_point_data[i].m_views.clear();
            m_point_data[i].m_color[0] = 0x0;
            m_point_data[i].m_color[1] = 0x0;
            m_point_data[i].m_color[2] = 0xff;

            num_pruned++;
        }
    }

    printf("[PruneBadPoints] Pruned %d points\n", num_pruned);
}

void BundlerApp::ZeroDistortionParams()
{
    int num_cameras = GetNumImages();
    
    for (int i = 0; i < num_cameras; i++) {
        m_image_data[i].m_camera.m_k[0] = 0.0;
        m_image_data[i].m_camera.m_k[1] = 0.0;
    }
}

void BundlerApp::OutputRelativePoses2D(const char *outfile)
{
    FILE *f = fopen(outfile, "w");

    if (f == NULL) {
        printf("[OutputRelativePose2D] Error opening file %s for writing\n",
               outfile);
        return;
    }
    
    int num_images = GetNumImages();
    SetMatchesFromPoints(16);
    int num_pairs = 0;
    
    for (int i = 0; i < num_images; i++) {
        // std::list<unsigned int> &nbrs = m_matches.GetNeighbors(i);
        // int num_nbrs = (int) nbrs.size();
        int num_nbrs = (int) m_matches.GetNumNeighbors(i);
        num_pairs += num_nbrs;
    }

    // num_pairs = num_pairs / 2;

    fprintf(f, "%d %d\n", num_images, num_pairs);

    FILE *f_gt = fopen("groundtruth.txt", "w");
    fprintf(f_gt, "%d\n", num_images);

    for (int i = 0; i < num_images; i++) {
        if (!m_image_data[i].m_camera.m_adjusted) {
            fprintf(f_gt, "0.0 0.0 0.0\n");
            continue;
        }

        /* Get the position and view direction */
        double pos[3], view[3];

        m_image_data[i].m_camera.GetPosition(pos);
        m_image_data[i].m_camera.GetViewDirection(view);

        v2_t view2D = v2_new(view[0], -view[2]);
        view2D = v2_unit(view2D);

        double view_angle = atan2(Vy(view2D), Vx(view2D));
        fprintf(f_gt, "%0.6f %0.6f %0.6f\n", pos[0], -pos[2], view_angle);

        /* Get the neighbors */
        // std::list<unsigned int> &nbrs = m_matches.GetNeighbors(i);
        MatchAdjList &nbrs = m_matches.GetNeighbors(i);
        int num_nbrs = (int) nbrs.size();

        // std::list<unsigned int>::iterator iter;
        MatchAdjList::iterator iter;
        for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
            int nbr = iter->m_index; // *iter;

            // if (nbr <= i)
            //     continue; /* only print neighbors once */

            double pos2[3], view2[3];

            m_image_data[nbr].m_camera.GetPosition(pos2);
            m_image_data[nbr].m_camera.GetViewDirection(view2);            

            v2_t view2D_2 = v2_new(view2[0], -view2[2]);
            view2D_2 = v2_unit(view2D_2);
            
            double dx = pos2[0] - pos[0];
            double dz = pos[2] - pos2[2];
            
            double angle = v2_angle(view2D, view2D_2);
            // double angle2 = v2_angle(view2D_2, view2D);

            /* Rotate dx,dz so that the first camera is pointed in the
             * positive z direction */
            
            double R[4] = { Vy(view2D), -Vx(view2D),
                            Vx(view2D), Vy(view2D) };
            double d[2] = { dx, dz };
            double Rd[2];
            matrix_product(2, 2, 2, 1, R, d, Rd);
            // double norm = matrix_norm(2, 1, Rd);
            // matrix_scale(2, 1, Rd, 1.0 / norm, Rd);

            fprintf(f, "%d %d %0.6f %0.6f %0.6f\n", 
                    i, nbr, Rd[0], Rd[1], angle);
        }
    }

    fclose(f_gt);
}

static void WriteVector(FILE *f, int n, const double *v)
{
    for (int i = 0; i < n; i++) {
        fprintf(f, "%0.16e ", v[i]);
    }
    fprintf(f, "\n");
}

static double GetTwist(double *R)
{
    double c_twist = 
        (R[0] * R[8] - R[6] * R[2]) / sqrt(1 - R[5] * R[5]);
    
    c_twist = CLAMP(c_twist, -1.0 + 1.0e-8, 1.0 - 1.0e-8);

    double angle = acos(c_twist);  

    if (R[3] < 0.0)
        return -angle;
    else 
        return angle;    
}

void BundlerApp::OutputRelativePoses3D(const char *outfile)
{
    FILE *f = fopen(outfile, "w");

    if (f == NULL) {
        printf("[OutputRelativePose3D] Error opening file %s for writing\n",
               outfile);
        return;
    }
    
    int num_images = GetNumImages();
    SetMatchesFromPoints(16);
    int num_pairs = 0;
    
    for (int i = 0; i < num_images; i++) {
        // std::list<unsigned int> &nbrs = m_matches.GetNeighbors(i);
        // int num_nbrs = (int) nbrs.size();
        int num_nbrs = (int) m_matches.GetNumNeighbors(i);
        num_pairs += num_nbrs;
    }

    // num_pairs = num_pairs / 2;

    fprintf(f, "%d %d\n", num_images, num_pairs);

    FILE *f_gt = fopen("groundtruth.txt", "w");
    fprintf(f_gt, "%d\n", num_images);

    for (int i = 0; i < num_images; i++) {
        if (!m_image_data[i].m_camera.m_adjusted) {
            fprintf(f_gt, "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n");
            fprintf(f_gt, "0.0 0.0 0.0\n");
            continue;
        }

        /* Get the position and view direction */
        double pos[3], pose[9];

        m_image_data[i].m_camera.GetPosition(pos);
        m_image_data[i].m_camera.GetPose(pose);

        WriteVector(f_gt, 9, pose);
        WriteVector(f_gt, 3, pos);

        /* Get the neighbors */
        MatchAdjList &nbrs = m_matches.GetNeighbors(i);
        int num_nbrs = (int) nbrs.size();

        MatchAdjList::iterator iter;
        for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
            int nbr = iter->m_index; // *iter;

            if (nbr <= i)
                continue; /* only print neighbors once */

            double pos2[3], pose2[9];

            m_image_data[nbr].m_camera.GetPosition(pos2);
            m_image_data[nbr].m_camera.GetPose(pose2);

            double diff[3];
            matrix_diff(3, 1, 3, 1, pos2, pos, diff);

            /* Put everything in the right coordinate system */
            double R[9], t[3];
            matrix_transpose_product(3, 3, 3, 3, pose, pose2, R);
            // matrix_transpose_product2(3, 3, 3, 3, pose2, pose, R);
            matrix_transpose_product(3, 3, 3, 1, pose, diff, t);

            double norm = matrix_norm(3, 1, t);
            matrix_scale(3, 1, t, 1.0 / norm, t);

            double viewdir[3] = { -R[2], -R[5], -R[8] };
            double twist_angle = GetTwist(R);

            fprintf(f, "%d %d\n", i, nbr);
            
            WriteVector(f, 9, R);
            // WriteVector(f, 3, viewdir);
            // fprintf(f, "%0.8f\n", twist_angle);
            WriteVector(f, 3, t);
        }
    }

    fclose(f_gt);
}
