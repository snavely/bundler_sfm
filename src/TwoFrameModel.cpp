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

/* TwoFrameModel.cpp */

#include <float.h>

#include "TwoFrameModel.h"
#include "ImageData.h"

#include "matrix.h"
#include "qsort.h"
#include "util.h"

#ifdef WIN32
#define isnan _isnan
#endif

static void ReadVector(FILE *f, int n, double *v)
{
    for (int i = 0; i < n; i++) {
        fscanf(f, "%lf ", v + i);
    }
}

static void ReadCamera(FILE *f, camera_params_t &camera)
{
    ReadVector(f, 9, camera.R);
    ReadVector(f, 3, camera.t);
    fscanf(f, "%lf\n", &(camera.f));
}

static void WriteVector(FILE *f, int n, const double *v)
{
    for (int i = 0; i < n; i++) {
        fprintf(f, "%0.16e ", v[i]);
    }
    fprintf(f, "\n");
}

static void WriteCamera(FILE *f, const camera_params_t &camera)
{
    WriteVector(f, 9, camera.R);
    WriteVector(f, 3, camera.t);
    fprintf(f, "%0.16e\n", camera.f);
}

static void WriteCameraPose(FILE *f, const camera_params_t &camera)
{
    WriteVector(f, 9, camera.R);
    WriteVector(f, 3, camera.t);
}

void TwoFrameModel::Read(FILE *f)
{
    fscanf(f, "%d\n", &m_num_points);

    fscanf(f, "%lf\n", &m_angle);
    fscanf(f, "%lf\n", &m_error);

    v3_t *points_tmp = new v3_t[m_num_points];
    double *tracks_tmp = new double[m_num_points];
    int *k1_tmp = new int[m_num_points];
    int *k2_tmp = new int[m_num_points];

    m_points = new v3_t[m_num_points];
    m_tracks = new int[m_num_points];

    for (int i = 0; i < m_num_points; i++) {
        int tr, k1, k2;
        fscanf(f, "%d %d %d %lf %lf %lf\n", &tr, &k1, &k2,
               &(Vx(points_tmp[i])), 
               &(Vy(points_tmp[i])), 
               &(Vz(points_tmp[i])));

        tracks_tmp[i] = (double) tr;
        k1_tmp[i] = k1;
        k2_tmp[i] = k2;
    }

    bool use_tracks = true;
    if (tracks_tmp[0] < 0)
        use_tracks = false;

    if (use_tracks) {
        int *perm = new int[m_num_points];
    
        qsort_ascending();
        qsort_perm(m_num_points, tracks_tmp, perm);

        for (int i = 0; i < m_num_points; i++) {
            m_tracks[i] = iround(tracks_tmp[i]);
            m_points[i] = points_tmp[perm[i]];
        }

        m_keys1 = m_keys2 = NULL;

        delete [] perm;
    } else {
        m_keys1 = new int[m_num_points];
        m_keys2 = new int[m_num_points];

        for (int i = 0; i < m_num_points; i++) {
            m_keys1[i] = k1_tmp[i];
            m_keys2[i] = k2_tmp[i];
            m_points[i] = points_tmp[i];
        }
     
        m_tracks = NULL;
    }

    delete [] points_tmp;
    delete [] tracks_tmp;
    delete [] k1_tmp;
    delete [] k2_tmp;

    ReadCamera(f, m_camera0);
    ReadCamera(f, m_camera1);
    
    ReadVector(f, 9, m_C0);
    ReadVector(f, 9, m_C1);
}

void TwoFrameModel::Write(FILE *f) const
{
    fprintf(f, "%d\n", m_num_points);

    fprintf(f, "%0.9f\n", m_angle);

    fprintf(f, "%0.9f\n", m_error);

    for (int i = 0; i < m_num_points; i++) {
        int tr = -1, k1 = -1, k2 = -1;
            
        if (m_tracks != NULL)
            tr = m_tracks[i];
        
        if (m_keys1 != NULL)
            k1 = m_keys1[i];
        
        if (m_keys2 != NULL)
            k2 = m_keys2[i];
        
        fprintf(f, "%d %d %d %0.16e %0.16e %0.16e\n", tr, k1, k2,
                Vx(m_points[i]), Vy(m_points[i]), Vz(m_points[i]));
    }
    
    WriteCamera(f, m_camera0);
    WriteCamera(f, m_camera1);
    
    WriteVector(f, 9, m_C0);
    WriteVector(f, 9, m_C1);
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

void TwoFrameModel::WriteSparse(FILE *f)
{
    // WriteCameraPose(f, m_camera0);
    // WriteCameraPose(f, m_camera1);

    /* Compute the camera pose of camera1 relative to camera0 */
    double pos0[3], pos1[3];
    
    // matrix_transpose_product(3, 3, 3, 1, m_camera0.R, m_camera0.t, pos0);
    // matrix_transpose_product(3, 3, 3, 1, m_camera1.R, m_camera1.t, pos1);
    memcpy(pos0, m_camera0.t, 3 * sizeof(double));
    memcpy(pos1, m_camera1.t, 3 * sizeof(double));

    // matrix_scale(3, 1, pos0, -1.0, pos0);
    // matrix_scale(3, 1, pos1, -1.0, pos1);
    
    double diff[3];
    matrix_diff(3, 1, 3, 1, pos1, pos0, diff);
    
    double R1[9], tr[3];
    matrix_transpose_product2(3, 3, 3, 3, m_camera0.R, m_camera1.R, R1);
    // matrix_transpose_product(3, 3, 3, 3, m_camera1.R, m_camera0.R, R1);
    matrix_product(3, 3, 3, 1, m_camera0.R, diff, tr);

    double norm = matrix_norm(3, 1, tr);
    matrix_scale(3, 1, tr, 1.0 / norm, tr);

    double viewdir[3] = { -R1[2], -R1[5], -R1[8] };
    double twist_angle = GetTwist(R1);
    
    /* Compute the distance to the scene */
    double z_avg = 0.0;
    
    for (int p = 0; p < m_num_points; p++) {
        v3_t &pt = m_points[p];

        double diff1[3], diff2[3];
        matrix_diff(3, 1, 3, 1, pt.p, pos0, diff1);
        matrix_diff(3, 1, 3, 1, pt.p, pos1, diff2);

        double dist1 = matrix_norm(3, 1, diff1);
        double dist2 = matrix_norm(3, 1, diff2);

        z_avg += 0.5 * (dist1 + dist2) / norm;
    }

    z_avg /= m_num_points;
    
    WriteVector(f, 9, R1);
    /* Write the viewing direction */
    // WriteVector(f, 3, viewdir);
    /* Write the twist angle */
    // fprintf(f, "%0.8f\n", twist_angle);
    /* Write the translation */
    WriteVector(f, 3, tr);

    fprintf(f, "%0.6f\n", z_avg);
}

void TwoFrameModel::WriteBrief(FILE *f) const
{
    fprintf(f, "%d\n", m_num_points);

    fprintf(f, "%0.5f\n", m_angle);

    fprintf(f, "%0.5f\n", m_error);
    
    WriteVector(f, 9, m_C0);
    WriteVector(f, 9, m_C1);
}

double TwoFrameModel::ComputeTrace(bool flip) const
{
#if 0
    if (flip) {
        return m_C0[0] + m_C0[4] + m_C0[8];
    } else {
        return m_C1[0] + m_C1[4] + m_C1[8];
    }
#else
    if (flip) {
        return m_C1[0] + m_C1[4] + m_C1[8];
    } else {
        return m_C0[0] + m_C0[4] + m_C0[8];
    }
#endif
}

double TwoFrameModel::AverageDistanceToPoints() const {
    double dist_sum = 0.0;
    const double *pos1 = m_camera0.t;
    const double *pos2 = m_camera1.t;
    
    for (int i = 0; i < m_num_points; i++) {
        double diff1[3], diff2[3];
        matrix_diff(3, 1, 3, 1, m_points[i].p, (double *) pos1, diff1);
        matrix_diff(3, 1, 3, 1, m_points[i].p, (double *) pos2, diff2);
        
        double norm1 = matrix_norm(3, 1, diff1);
        double norm2 = matrix_norm(3, 1, diff2);

        dist_sum += (norm1 + norm2);
    }

    return dist_sum / (2 * m_num_points);
}

#if 0
/* Turn this edge from a scaffolding edge to a connector edge */
void TwoFrameModel::TurnOffScaffold() 
{
    m_scaffold = false;
    // matrix_ident(4, m_Sinit0);
    // matrix_ident(4, m_Sinit1);
    m_scale0 = m_scale1 = 1.0;

    matrix_ident(3, m_C0);
    matrix_ident(3, m_C1);

    m_C0[0] = m_C1[0] = 1.0e10;
    m_C0[4] = m_C1[4] = 1.0e10;
    m_C0[8] = m_C1[8] = 1.0e10;
    m_angle = 0.0;
}
#endif

void TwoFrameModel::
ComputeTransformedCovariance(bool flip, double *S, double *C) const
{
    double tmp[9];
    if (flip) {
        matrix_product(3, 3, 3, 3, S, m_C1, tmp);
        matrix_transpose_product2(3, 3, 3, 3, tmp, S, C);
    } else {
        matrix_product(3, 3, 3, 3, S, m_C0, tmp);
        matrix_transpose_product2(3, 3, 3, 3, tmp, S, C);        
    }    
}

void TwoFrameModel::SetPrevious(bool flip, int prev)
{
    if (flip) {
        m_prev1 = prev;
    } else {
        m_prev0 = prev;
    }
}

int TwoFrameModel::GetPrevious(bool flip)
{
    if (flip)
        return m_prev1;
    else
        return m_prev0;
}

void TwoFrameModel::SetPredecessor(bool flip, int pred)
{
    if (flip) {
        m_pred1 = pred;
    } else {
        m_pred0 = pred;
    }
}

int TwoFrameModel::GetPredecessor(bool flip)
{
    if (flip)
        return m_pred1;
    else
        return m_pred0;
}

void TwoFrameModel::SetFlag(bool flip, int flag)
{
    if (flip)
        m_flag1 = flag;
    else
        m_flag0 = flag;
}

int TwoFrameModel::GetFlag(bool flip)
{
    if (flip)
        return m_flag1;
    else
        return m_flag0;
}

ModelMap ReadModels(FILE *f, int *num_images_out) 
{
    char buf[256];

    int num_images;
    fscanf(f, "%d\n", &num_images);

    ModelMap models(num_images);

    // assert(num_images == num_images);

    while (fgets(buf, 256, f)) {
        int i1, i2;
        sscanf(buf, "%d %d", &i1, &i2);

        TwoFrameModel m;
        m.Read(f);
        
        if (m.ComputeTrace(true) < 0.0 || m.ComputeTrace(false) < 0.0) {
            printf("[ReadModels] Error! Trace(%d,%d) < 0!\n", i1, i2);
            continue;
        }

        if (m.m_num_points < 28 /*33*/) {
            // printf("[ReadModels] Error! Too few points [%d] for (%d,%d)\n",
            //        m.m_num_points, i1, i2);

            continue;
        }

        if (std::isnan(m.m_angle) || std::isnan(m.m_error)) {
            printf("[ReadModels] Error! NaNs in pair %d,%d!\n", i1, i2);
            continue;
        }

        assert(i1 < i2);
        models.AddModel(GetMatchIndex(i1, i2), m);
    }

    if (num_images_out != NULL)
        *num_images_out = num_images;

    return models;
}

PEdgeMap ReadPEdges(FILE *f, int num_images) 
{
    char buf[256];

    PEdgeMap p_edges;

    while (fgets(buf, 256, f)) {
        int i1, i2;
        sscanf(buf, "%d %d", &i1, &i2);

        p_edges[i1 * num_images + i2] = true;
    }

    return p_edges;
}

void WriteModels(ModelMap &models, int num_images, char *out_file)
{
    FILE *f = fopen(out_file, "w");
    if (f == NULL) {
        printf("[WriteModels] Error opening file %s for reading\n", out_file);
    } else {
        fprintf(f, "%d\n", num_images);

        // FIXME LOOP
        for (unsigned int i = 0; i < (unsigned int) num_images; i++) {
            // for (int j = i+1; j < num_images; j++) {
            ModelTable::iterator iter;
            for (iter = models.Begin(i); iter != models.End(i); iter++) {
                unsigned int j = iter->first; // iter->m_index;

                if (i >= j)
                    continue;

                MatchIndex idx = GetMatchIndex(i, j);
                if (models.Contains(idx)) {
                    fprintf(f, "%d %d\n", i, j);
                    models.GetModel(idx).Write(f);
                }
            }
        }
        
        fclose(f);
    }
}

void WriteModelsSparse(ModelMap &models, int num_images, char *out_file)
{
    FILE *f = fopen(out_file, "w");
    if (f == NULL) {
        printf("[WriteModels] Error opening file %s for reading\n", out_file);
    } else {
        unsigned int num_pairs = 0;

        for (unsigned int i = 0; i < (unsigned int) num_images; i++) {
            num_pairs += models.GetNeighbors(i).size();
        }

        num_pairs = num_pairs / 2;

        fprintf(f, "%d %u\n", num_images, num_pairs);

        // FIXME LOOP
        for (unsigned int i = 0; i < (unsigned int) num_images; i++) {
            // for (int j = i+1; j < num_images; j++) {
            ModelTable::iterator iter;
            for (iter = models.Begin(i); iter != models.End(i); iter++) {
                unsigned int j = iter->first; // iter->m_index;

                if (i >= j)
                    continue;

                MatchIndex idx = GetMatchIndex(i, j);
                if (models.Contains(idx)) {
                    fprintf(f, "%d %d\n", i, j);
                    models.GetModel(idx).WriteSparse(f);
                }
            }
        }
        
        fclose(f);
    }
}

void WritePEdges(PEdgeMap &p_edges, int num_images, char *out_file)
{
    FILE *f = fopen(out_file, "w");
    if (f == NULL) {
        printf("[WritePEdges] Error opening file %s for reading\n", out_file);
    } else {
        for (int i = 0; i < num_images; i++) {
            for (int j = i+1; j < num_images; j++) {
                if (p_edges.find(i * num_images + j) != p_edges.end()) {
                    fprintf(f, "%d %d\n", i, j);
                }
            }
        }
        
        fclose(f);
    }
}

void ThresholdTwists(int num_images, ModelMap &models, 
                     std::vector<ImageData> &image_data, bool panos_only) 
{
    int *num_large_twists = new int[num_images];
    int *degree = new int[num_images];

    for (int i = 0; i < num_images; i++) {
        num_large_twists[i] = 0;
        degree[i] = 0;
    }

    for (int i = 0; i < num_images; i++) {
        ModelTable::iterator iter;
        for (iter = models.Begin(i); iter != models.End(i); iter++) {
            unsigned int j = iter->first; // iter->m_index;
            
            if (i >= j)
                continue;
            
            MatchIndex idx = GetMatchIndex(i, j);
            if (models.Contains(idx)) {
                TwoFrameModel &m = models.GetModel(idx);
                
                /* Compute the twist */
                double Rp_i[9], Rp_j[9];
                matrix_transpose(3, 3, m.m_camera0.R, Rp_i);
                matrix_transpose(3, 3, m.m_camera1.R, Rp_j);
                
                double Rp_ij[9];
                matrix_transpose_product(3, 3, 3, 3, Rp_i, Rp_j, Rp_ij);
                
                double twist_angle = GetTwist(Rp_ij);
                
                if (fabs(RAD2DEG(twist_angle)) >= 12.0) { 
                    num_large_twists[i]++;
                    num_large_twists[j]++;
                }
                
                degree[i]++;
                degree[j]++;
            }
        }
    }
    
    for (int i = 0; i < num_images; i++) {
        if (degree[i] == 0)
            continue;

        double perc_large_twists = (double) num_large_twists[i] / degree[i];
        
        int w = image_data[i].GetWidth();
        int h = image_data[i].GetHeight();
        
        double ratio = (double) w / h;

        if ((panos_only || perc_large_twists < 0.4) && 
             ratio > 0.4 && ratio < 2.5) {
            continue;
        }

        printf("[ThresholdTwists] Removing image %d with score %0.3f, %0.3f\n",
               i, perc_large_twists, ratio);

        std::list<unsigned int> nbrs = models.GetNeighbors(i);
        std::list<unsigned int>::iterator iter;
        for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
            unsigned int j = *iter; // iter->m_index;
            
            if (i < j) {            
                MatchIndex idx = GetMatchIndex(i, j);
                if (models.Contains(idx)) {
                    models.RemoveModel(idx);
                }
            } else {
                MatchIndex idx = GetMatchIndex(j, i);
                if (models.Contains(idx)) {
                    models.RemoveModel(idx);
                }
            }
        }
    }
}
