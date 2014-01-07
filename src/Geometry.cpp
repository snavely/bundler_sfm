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

/* Geometry.cpp */
/* Geometric primitives */

#include <assert.h>
#include <float.h>
#include <math.h>


#include <algorithm>
#include <vector>

#include "Geometry.h"
#include "Camera.h"

#ifndef __BUNDLER__
#include "UtilGL.h"
#endif /* __BUNDLER__ */

#include "defines.h"
#include "fit.h"
#include "matrix.h"
#include "util.h"

#ifndef __DEMO__
#ifndef __BUNDLER__
#include "AHStroke/color.h"
#include "AHStroke/Stroke.h"
#endif
#endif



/*Write the Coordinates*/
void PointData::WriteCoordinates(FILE *f)
{
	fprintf(f,"%lf %lf %lf\n",m_pos[0],m_pos[1],m_pos[2]);
}


/* Write the point data in XML */
void PointData::WriteXML(FILE *f)
{
    static const char *spacer = "    ";

    /* Position and color */
    fprintf(f, "%s<point>\n", spacer);
    fprintf(f, "%s  <pos>\n", spacer);
    fprintf(f, "%s    <x> %0.8e </x>\n", spacer, m_pos[0]);
    fprintf(f, "%s    <y> %0.8e </y>\n", spacer, m_pos[1]);
    fprintf(f, "%s    <z> %0.8e </z>\n", spacer, m_pos[2]);
    fprintf(f, "%s  </pos>\n", spacer);
    fprintf(f, "%s  <col>\n", spacer);
    fprintf(f, "%s    <r> %d </r>\n", spacer, iround(m_color[0]));
    fprintf(f, "%s    <g> %d </g>\n", spacer, iround(m_color[1]));
    fprintf(f, "%s    <b> %d </b>\n", spacer, iround(m_color[2]));
    fprintf(f, "%s  </col>\n", spacer);

    /* Views */
    fprintf(f, "%s  <views>\n", spacer);

    int num_views = (int) m_views.size();
    for (int i = 0; i < num_views; i++) {
	fprintf(f, "%s    <view>\n", spacer);
	fprintf(f, "%s      <cam> %d </cam>\n", spacer, m_views[i].first);
	// fprintf(f, "%s      <key> %d </key>\n", spacer, m_views[i].second);
	fprintf(f, "%s    </view>\n", spacer);
    }

    fprintf(f, "%s  </views>\n", spacer);
    fprintf(f, "%s</point>\n", spacer);
}

/* Write the point data in XML */
void PointData::WriteGeoXML(FILE *f)
{

}

/* Create a planar patch for this point */
void PointData::CreatePlanarPatch(double size, PlaneData &plane)
{
    memcpy(plane.m_normal, m_norm, 3 * sizeof(double));

    double dist;
    matrix_product(1, 3, 3, 1, m_norm, m_pos, &dist);
    plane.m_dist = -dist;

    memcpy(plane.m_origin, m_pos, 3 * sizeof(double));

    // u-axis should be perp to the y-axis and to normal
    double y[3] = { 0.0, 1.0, 0.0 };
    matrix_cross(m_norm, y, plane.m_uaxis);

    // v-axis should be perp to normal and uaxis
    matrix_cross(m_norm, plane.m_uaxis, plane.m_vaxis);

    plane.m_u_extent = plane.m_v_extent = size;
}

/* Read/write routines for tracks */
void TrackData::Read(FILE *f) 
{
    int size;
    fscanf(f, "%d", &size);
    
    for (int i = 0; i < size; i++) {
	ImageKey ik;
	fscanf(f, "%d %d", &(ik.first), &(ik.second));
	m_views.push_back(ik);
    }
}

void TrackData::Write(FILE *f) 
{
    int size = (int) m_views.size();
    fprintf(f, "%d", size);
    
    for (int i = 0; i < size; i++) {
	fprintf(f, " %d %d", m_views[i].first, m_views[i].second);
    }

    fprintf(f, "\n");
}


static int cube_faces[24] = 
    { 0, 1, 3, 2,
      0, 1, 5, 4, 
      0, 4, 6, 2,
      1, 3, 7, 5, 
      3, 2, 6, 7, 
      4, 5, 7, 6 };

void Cube::Finalize() 
{
    int count = 0;
    for (int i = 0; i < 2; i++) {
	int x_sign = (i == 0) ? -1 : 1;

	for (int j = 0; j < 2; j++) {
	    int y_sign = (j == 0) ? -1 : 1;

	    for (int k = 0; k < 2; k++) {
		int z_sign = (k == 0) ? -1 : 1;
		
		double x[3], y[3], z[3];
		matrix_scale(3, 1, m_x_axis, x_sign * 0.5 * m_x_scale, x);
		matrix_scale(3, 1, m_y_axis, y_sign * 0.5 * m_y_scale, y);
		matrix_scale(3, 1, m_z_axis, z_sign * 0.5 * m_z_scale, z);
		
		m_vertices[3 * count + 0] = m_origin[0] + x[0] + y[0] + z[0];
		m_vertices[3 * count + 1] = m_origin[1] + x[1] + y[1] + z[1];
		m_vertices[3 * count + 2] = m_origin[2] + x[2] + y[2] + z[2];

		count++;
	    }
	}
    }
}

void Cube::Render()
{
#ifndef __BUNDLER__
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(4.0);

    for (int i = 0; i < 6; i++) {
	glBegin(GL_LINE_LOOP);
	int vidx0 = cube_faces[i * 4 + 0];
	int vidx1 = cube_faces[i * 4 + 1];
	int vidx2 = cube_faces[i * 4 + 2];
	int vidx3 = cube_faces[i * 4 + 3];

	glVertex3dv(m_vertices + 3 * vidx0);
	glVertex3dv(m_vertices + 3 * vidx1);
	glVertex3dv(m_vertices + 3 * vidx2);
	glVertex3dv(m_vertices + 3 * vidx3);

	glEnd();
    }
#endif /* __BUNDLER__ */
}

/* Project a point onto the plane */
void PlaneData::Project(double *p, double *p_proj) 
{
    /* Subtract the distance vector */
    double vec[3];    
    matrix_scale(3, 1, m_normal, m_dist, vec);
    
    double p_norm[3];
    matrix_diff(3, 1, 3, 1, p, vec, p_norm);

    double dot;
    matrix_product(1, 3, 3, 1, m_normal, p_norm, &dot);
    
    double p_par[3];
    matrix_scale(3, 1, m_normal, dot, p_par);

    double p_perp[3];
    matrix_diff(3, 1, 3, 1, p_norm, p_par, p_perp);
    
    matrix_sum(3, 1, 3, 1, p_perp, vec, p_proj);
}

double PlaneData::ProjectU(double *p, double *p_proj) 
{
    double dot;
    matrix_product(1, 3, 3, 1, p, m_uaxis, &dot);
    matrix_scale(3, 1, m_uaxis, dot, p_proj);

    return dot;
}

double PlaneData::ProjectV(double *p, double *p_proj) 
{
    double dot;
    matrix_product(1, 3, 3, 1, p, m_vaxis, &dot);
    matrix_scale(3, 1, m_vaxis, dot, p_proj);
    return dot;
}

void PlaneData::Transform(const double *M)
{
    double p[4] = { m_normal[0], m_normal[1], m_normal[2], m_dist };

    double Minv[16];
    matrix_invert(4, (double *)M, Minv);
    
    double pNew[4];
    matrix_transpose_product(4, 4, 4, 1, Minv, p, pNew);

    double len = matrix_norm(3, 1, pNew);
    
    m_normal[0] = pNew[0] / len;
    m_normal[1] = pNew[1] / len;
    m_normal[2] = pNew[2] / len;
    m_dist = pNew[3] / len;

#if 0
    double origin[4] = { m_origin[0], m_origin[1], m_origin[2], 1.0 };
    double Morigin[4];
    matrix_product(4, 4, 4, 1, (double *) M, origin, Morigin);

    double dot0, dot1;
    matrix_product(1, 4, 4, 1, p, origin, &dot0);
    matrix_product(1, 4, 4, 1, pNew, Morigin, &dot1);

    printf("dot0 = %0.3f\n", dot0);
    printf("dot1 = %0.3f\n", dot1);
#endif
}

/* Setup various planar aspect */
void PlaneData::Setup(std::vector<PointData> &point_data, 
		      double *origin, double *up)
{
    /* Compute the mean of the points on the plane */
    double mean[3] = { 0.0, 0.0, 0.0 };
    int num_plane_points = (int) m_points.size();

    for (int i = 0; i < num_plane_points; i++) {
        int pt_idx = m_points[i];
        mean[0] += point_data[pt_idx].m_pos[0];
        mean[1] += point_data[pt_idx].m_pos[1];
        mean[2] += point_data[pt_idx].m_pos[2];
    }

    matrix_scale(3, 1, mean, 1.0 / num_plane_points, mean);
    matrix_diff(3, 1, 3, 1, mean, origin, mean);

    /* Project the mean onto the plane */
    Project(mean, m_origin);

    matrix_sum(3, 1, 3, 1, m_origin, origin, m_origin);

    memcpy(m_vaxis, up, 3 * sizeof(double));
    matrix_cross(m_normal, m_vaxis, m_uaxis);

    /* Compute the variance in each direction */
    double u_variance = 0.0, v_variance = 0.0;
    for (int i = 0; i < num_plane_points; i++) {
        int pt_idx = m_points[i];

        /* Subtract out the origin */
        double diff[3];
        matrix_diff(3, 1, 3, 1, point_data[pt_idx].m_pos, m_origin, diff);

        /* Project onto u,v axes */
        double u_proj[3], v_proj[3];
        ProjectU(diff, u_proj);
        ProjectV(diff, v_proj);

        u_variance += matrix_normsq(3, 1, u_proj);
        v_variance += matrix_normsq(3, 1, v_proj);
    }

    u_variance = sqrt(u_variance / num_plane_points);
    v_variance = sqrt(v_variance / num_plane_points);

    m_u_extent = 2.0 * u_variance;
    m_v_extent = 2.0 * v_variance;


    /* Find all the views that see this plane */
    std::vector<int> views;
    for (int i = 0; i < num_plane_points; i++) {
        int pt_idx = m_points[i];

        for (int j = 0; j < (int) point_data[pt_idx].m_views.size(); j++) {
            int view = point_data[pt_idx].m_views[j].first;

            bool found = false;
            for (int k = 0; k < (int) views.size(); k++) {
                if (views[k] == view) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                views.push_back(view);
            }
        }
    }

    m_views = views;

}


PlaneData::PlaneData(double normal[3],double dist)
{
	for(int i=0;i<3;i++)
		m_normal[i]=normal[i];
	m_dist=dist;
}

/* check ifthe point lies inside the plane*/

bool PlaneData::CheckInside(double *point)
{
	double side_vec[3],poly_vec[3],point_vec[3],cross1[3],cross2[3];
	for(int i=0;i<4;i++)
	{
		int i1,i2,i3;
		i1=i;
		i2=(i1+1)%4;
		i3=(i2+1)%4;
		matrix_diff(3,1,3,1,m_corners+i3*3,m_corners+i2*3,side_vec);
		matrix_diff(3,1,3,1,m_corners+i1*3,m_corners+i2*3,poly_vec);
		matrix_diff(3,1,3,1,point,m_corners+i2*3,point_vec);
		matrix_cross(side_vec,poly_vec,cross1);
		matrix_cross(side_vec,point_vec,cross2);
		double dp=0;
		for(int j=0;j<3;dp+=cross1[j]*cross2[j],j++);
		if(dp<0) return false;
	}
	return true;
}

/* Intersect a ray with the plane */
double PlaneData::IntersectRay(double *pos, double *ray, double *pt) const
{
    double pos_dot, ray_dot;

    matrix_product(1, 3, 3, 1, (double *) m_normal, pos, &pos_dot);
    matrix_product(1, 3, 3, 1, (double *) m_normal, ray, &ray_dot);

    if (ray_dot==0.0) {
        return -DBL_MAX;
    }

    double t = (-m_dist - pos_dot) / ray_dot;

    pt[0] = pos[0] + t * ray[0];
    pt[1] = pos[1] + t * ray[1];
    pt[2] = pos[2] + t * ray[2];

    return t;
}

bool PlaneData::SetCorners(double *R, double *t, double f, int w, int h, int ymin, int ymax, int xmin, int xmax)
{
	if(xmax==-1)  xmax=w;
	if(ymax==-1)  ymax=h;
	if(xmax>w || ymax>h ||xmax<-1 || ymax<-1) printf("Error in bounding box data w=%d h=%d xmin=%d xmax=%d ymin=%d ymax=%d\n",w,h,xmin,xmax,ymin,ymax);
	//printf("xmin=%d xmax=%d ymin=%d ymax=%d\n",xmin,xmax,ymin,ymax);
    double Rt[9];
    matrix_transpose(3, 3, R, Rt);

    double center[3];
    matrix_product(3, 3, 3, 1, Rt, t, center);
    matrix_scale(3, 1, center, -1.0, center);

    /* Create rays for the four corners */
	//uncomment the follwing code to use the 4 image corners
    //double ray1[3] = { -0.5 * w, -0.5 * h, -f };
    //double ray2[3] = {  0.5 * w, -0.5 * h, -f };
    //double ray3[3] = {  0.5 * w,  0.5 * h, -f };
    //double ray4[3] = { -0.5 * w,  0.5 * h, -f };

	double ray1[3] = { xmin-0.5 * w, ymin-0.5 * h, -f };
	double ray2[3] = {  xmax-0.5 * w, ymin-0.5 * h, -f };
	double ray3[3] = {  xmax -0.5 * w,  ymax-0.5 * h, -f };
	double ray4[3] = {xmin -0.5 * w,  ymax - 0.5 * h, -f };

    double ray_world[18];
    matrix_product(3, 3, 3, 1, Rt, ray1, ray_world + 0);
    matrix_product(3, 3, 3, 1, Rt, ray2, ray_world + 3);
    matrix_product(3, 3, 3, 1, Rt, ray3, ray_world + 6);
    matrix_product(3, 3, 3, 1, Rt, ray4, ray_world + 9);

    double t0 = IntersectRay(center, ray_world + 0, m_corners + 0);
    double t1 = IntersectRay(center, ray_world + 3, m_corners + 3);
    double t2 = IntersectRay(center, ray_world + 6, m_corners + 6);
    double t3 = IntersectRay(center, ray_world + 9, m_corners + 9);

    if (t0 > 0.0 && t1 > 0.0 && t2 > 0.0 && t3 > 0.0)
        return true;
    else
        return false;
}

void PlaneData::Read(FILE *f) 
{
    /* Plane parameters */
    fscanf(f, "%lf %lf %lf %lf", 
        m_normal + 0, m_normal + 1, m_normal + 2, &m_dist);

    /* Plane normal */
    fscanf(f, "%lf %lf %lf", m_origin + 0, m_origin + 1, m_origin + 2);

    /* U, V axes */
    fscanf(f, "%lf %lf %lf", m_uaxis + 0, m_uaxis + 1, m_uaxis + 2);
    fscanf(f, "%lf %lf %lf", m_vaxis + 0, m_vaxis + 1, m_vaxis + 2);

    fscanf(f, "%lf %lf %lf %lf", 
        m_color + 0, m_color + 1, m_color + 2, m_color + 3);

    fscanf(f, "%lf %lf", &m_u_extent, &m_v_extent);

    int num_views;
    fscanf(f, "%d", &num_views);
    for (int i = 0; i < num_views; i++) {
        int view;
        fscanf(f, "%d", &view);
        m_views.push_back(view);
    }

    int num_points;
    fscanf(f, "%d", &num_points);
    for (int i = 0; i < num_points; i++) {
        int point;
        fscanf(f, "%d", &point);
        m_points.push_back(point);
    }
    fprintf(f, "\n");    
}



void PlaneData::Write(FILE *f) 
{
    /* Plane parameters */
    fprintf(f, "%0.6lf %0.6lf %0.6lf %0.6lf\n", 
	    m_normal[0], m_normal[1], m_normal[2], m_dist);

    /* Plane normal */
    fprintf(f, "%0.6lf %0.6lf %0.6lf\n", 
	    m_origin[0], m_origin[1], m_origin[2]);

    /* U, V axes */
    fprintf(f, "%0.6lf %0.6lf %0.6lf\n", m_uaxis[0], m_uaxis[1], m_uaxis[2]);
    fprintf(f, "%0.6lf %0.6lf %0.6lf\n", m_vaxis[0], m_vaxis[1], m_vaxis[2]);

    fprintf(f, "%0.6lf %0.6lf %0.6lf %0.6lf\n", 
	    m_color[0], m_color[1], m_color[2], m_color[3]);

    fprintf(f, "%0.6lf %0.6lf\n", m_u_extent, m_v_extent);

    fprintf(f, "%d\n", (int) m_views.size());
    for (int i = 0; i < (int) m_views.size(); i++) {
	fprintf(f, "%d ", m_views[i]);
    }
    fprintf(f, "\n");

    fprintf(f, "%d\n", (int) m_points.size());
    for (int i = 0; i < (int) m_points.size(); i++) {
	fprintf(f, "%d ", m_points[i]);
    }
    fprintf(f, "\n");
}

void PlaneData::Render() {
#ifndef __BUNDLER__
    GLCheckForError("PlaneData::Render[begin]");

    double p0[3], p1[3], p2[3], p3[3];
    double u_scale[3], v_scale[3];

    memcpy(u_scale, m_uaxis, 3 * sizeof(double));
    memcpy(v_scale, m_vaxis, 3 * sizeof(double));

    matrix_scale(3, 1, u_scale, m_u_extent, u_scale);
    matrix_scale(3, 1, v_scale, m_v_extent, v_scale);

    p0[0] = m_origin[0] + u_scale[0] + v_scale[0];
    p0[1] = m_origin[1] + u_scale[1] + v_scale[1];
    p0[2] = m_origin[2] + u_scale[2] + v_scale[2];

    p1[0] = m_origin[0] + u_scale[0] - v_scale[0];
    p1[1] = m_origin[1] + u_scale[1] - v_scale[1];
    p1[2] = m_origin[2] + u_scale[2] - v_scale[2];

    p2[0] = m_origin[0] - u_scale[0] - v_scale[0];
    p2[1] = m_origin[1] - u_scale[1] - v_scale[1];
    p2[2] = m_origin[2] - u_scale[2] - v_scale[2];

    p3[0] = m_origin[0] - u_scale[0] + v_scale[0];
    p3[1] = m_origin[1] - u_scale[1] + v_scale[1];
    p3[2] = m_origin[2] - u_scale[2] + v_scale[2];

    bool texture = (m_texture_index != -1);

    if (texture) {
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, m_texture_index);
    }

#if 0
    if (!texture)
	glColor4dv(m_color);
    else
	glColor4d(1.0, 1.0, 1.0, 0.3);
#else
    glColor4dv(m_color);
#endif

    double minx = m_bounds.m_min_x;
    double maxx = m_bounds.m_max_x;
    double miny = m_bounds.m_min_y;
    double maxy = m_bounds.m_max_y;

    glBegin(GL_QUADS);

    if (texture) 
	glTexCoord2d(maxx, maxy);
    glVertex3d(p0[0], p0[1], p0[2]); 

    if (texture) 
	glTexCoord2d(maxx, miny);
    glVertex3d(p1[0], p1[1], p1[2]); 

    if (texture) 
	glTexCoord2d(minx, miny);
    glVertex3d(p2[0], p2[1], p2[2]);

    if (texture) 
	glTexCoord2d(minx, maxy);
    glVertex3d(p3[0], p3[1], p3[2]);
    glEnd();

    if (texture)
	glDisable(GL_TEXTURE_2D);

    GLCheckForError("PlaneData::Render[end]");
#endif /* __BUNDLER__ */
}

/* Return the plane parameters */
void PlaneData::GetParams(double *params)
{
    params[0] = m_normal[0];
    params[1] = m_normal[1];
    params[2] = m_normal[2];
    params[3] = m_dist;    
}

/* Return the length of the segment */
double LineSegment2D::Length() {
    double dx = m_p2[0] - m_p1[0];
    double dy = m_p2[1] - m_p1[1];

    return sqrt(dx * dx + dy * dy);
}

/* Returns true (and stores the endpoints in t and u) if the
 * epipolar swath e1, e2 intersects this segment */
bool LineSegment2D::IntersectsEpipolarSwath(double e1[3], double e2[3], 
					    double &t, double &u)
{
    /* Find the line between the two endpoints */
    double p[3] = { m_p1[0], m_p1[1], 1.0 };
    double q[3] = { m_p2[0], m_p2[1], 1.0 };

    double l[3];
    matrix_cross(p, q, l);

    /* Intersect the line with the two epipolar lines */
    double i1[3], i2[3];

    matrix_cross(l, e1, i1);
    matrix_cross(l, e2, i2);

    double inv_i12 = 1.0 / i1[2];
    i1[0] *= inv_i12;
    i1[1] *= inv_i12;

    double inv_i22 = 1.0 / i2[2];
    i2[0] *= inv_i22;
    i2[1] *= inv_i22;

    double qp[3], i1p[3], i2p[3];
    matrix_diff(2, 1, 2, 1, q, p, qp);

    matrix_diff(2, 1, 2, 1, i1, p, i1p);
    matrix_diff(2, 1, 2, 1, i2, p, i2p);

    double dot1, dot2;
    matrix_product(1, 2, 2, 1, i1p, qp, &dot1);
    matrix_product(1, 2, 2, 1, i2p, qp, &dot2);

    int sign1 = SGN(dot1);
    int sign2 = SGN(dot2);

    double inv_norm = 1.0 / matrix_norm(2, 1, qp);
    double mag1 = matrix_norm(2, 1, i1p) * inv_norm; // matrix_norm(2, 1, qp);
    double mag2 = matrix_norm(2, 1, i2p) * inv_norm; // matrix_norm(2, 1, qp);

    t = sign1 * mag1;
    u = sign2 * mag2;

    /* Intersection check */
    if ((t < 0.0 && u < 0.0) || (t > 1.0 && u > 1.0))
	return false;
    else
	return true;
}

/* Return the sample at the given parameter value 
 * (l(0) = p1, l(1) = p2) */
void LineSegment2D::Sample(double t, double *p) 
{
    p[0] = (1.0 - t) * m_p1[0] + t * m_p2[0];
    p[1] = (1.0 - t) * m_p1[1] + t * m_p2[1];
}

/* Return the homogeneous form of the line */
void LineSegment2D::Homogeneous(double *l)
{
    double p1[3] = { m_p1[0], m_p1[1], 1.0 };
    double p2[3] = { m_p2[0], m_p2[1], 1.0 };    

    matrix_cross(p1, p2, l);
    
    double norm = matrix_norm(3, 1, l);
    matrix_scale(3, 1, l, 1.0 / norm, l);
}

void LineSegment2D::Read(FILE *f) 
{
    fscanf(f, "%lf %lf %lf %lf",
	   m_p1 + 0, m_p1 + 1, m_p2 + 0, m_p2 + 1);
}

void LineSegment2D::Write(FILE *f) 
{
    fprintf(f, "%0.8le %0.8le %0.8le %0.8le\n",
	    m_p1[0], m_p1[1], m_p2[0], m_p2[1]);
}

void LineSegment3D::Read(FILE *f) 
{
    fscanf(f, "%lf %lf %lf %lf %lf %lf",
	   m_p1 + 0, m_p1 + 1, m_p1 + 2, m_p2 + 0, m_p2 + 1, m_p2 + 2);

    /* Read visible views */
    int num_views;
    fscanf(f, "%d", &num_views);

    for (int i = 0; i < num_views; i++) {
	int v;
	fscanf(f, "%d", &v);
	m_views.push_back(v);
    }

    m_ignore = false;
}

void LineSegment3D::Write(FILE *f) 
{
    fprintf(f, "%0.8le %0.8le %0.8le %0.8le %0.8le %0.8le\n",
	    m_p1[0], m_p1[1], m_p1[2], m_p2[0], m_p2[1], m_p2[2]);

    /* Views */
    int num_views = (int) m_views.size();
    fprintf(f, "%d\n", num_views);

    for (int i = 0; i < num_views; i++)
	fprintf(f, "%d ", m_views[i]);

    fprintf(f, "\n");
}

void LineSegment3D::Render(const CameraInfo &camera, double max_width,
			   int stroke_texture, ParameterBound stroke_bounds)
{
#ifndef __BUNDLER__
#ifndef __DEMO__
#if 1
    /* Project the line into the camera */
    double p1[4] = { m_p1[0], m_p1[1], m_p1[2], 1.0 };
    double p2[4] = { m_p2[0], m_p2[1], m_p2[2], 1.0 };    
    double proj1[3], proj2[3];

    matrix_product(3, 4, 4, 1, (double *) camera.m_Pmatrix, p1, proj1);
    matrix_product(3, 4, 4, 1, (double *) camera.m_Pmatrix, p2, proj2);

    if (proj1[2] >= 0.0 || proj2[2] >= 0.0)
	return;

    double width = CLAMP(5.0 / (-proj1[2] - proj2[2]), 0.5, max_width);

    proj1[0] /= -proj1[2];
    proj1[1] /= -proj1[2];

    proj2[0] /= -proj2[2];
    proj2[1] /= -proj2[2];
#endif

#if 0
    bool in_front1 = camera.Project(m_p1, proj1);
    bool in_front2 = camera.Project(m_p2, proj2);

    if (!in_front1 || !in_front2)
	return;
#endif    

    /* Draw a line segment whose width is proportional to the distance
     * from the camera */

#if 0
    double pos[3];
    camera.GetPosition(pos);

    double disp[3];
    matrix_diff(3, 1, 3, 1, pos, m_p1, disp);
    
    double dist = matrix_norm(3, 1, disp);
#endif

    if (stroke_texture != -1) {
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, stroke_texture);
    }

#if 0
    #define LINE_WIDTH_SIGMA 1.0
    double width =
	max_width * exp(-dist * dist / (LINE_WIDTH_SIGMA * LINE_WIDTH_SIGMA));
#endif

    /* Create a stroke */
    Stroke s;
    s.radius() = 0.5 * width;
    s.cap() = false;
    s.depth() = 1.0;

    if (stroke_texture == -1) {
	s.useTexture() = false;
	s.color() = makeColor(0, 0, 0, 0.9f);
    } else {
	GLubyte b = 0x0;
	s.useTexture() = true;
	s.color() = makeColor(b, b, b, 0.9f);
    }

    s.addControlPoint(proj1[0], proj1[1]);
    s.addControlPoint(proj2[0], proj2[1]);

    s.render();

    if (stroke_texture != -1) {
		glDisable(GL_TEXTURE_2D);
    }
#endif /* __DEMO__ */
#endif /* __BUNDLER__ */
}

/* Return a unit vector pointing in the direction of the given line */
void LineToUnitVector(double *l, double *v) 
{
    v[0] = -l[1];
    v[1] = l[0];
    
    double mag = sqrt(v[0] * v[0] + v[1] * v[1]);
    
    v[0] /= mag;
    v[1] /= mag;

    v[2] = 0.0;
}

v3_t FindWeightedRobustMean(const std::vector<v3_t> &points,
			    const std::vector<double> &weights)
{
    int num_points = (int) points.size();
    double best_sum = DBL_MAX;
    int best_idx = -1;

    for (int i = 0; i < num_points; i++) {
	double sum = 0.0;

	for (int j = 0; j < num_points; j++) {
	    v3_t diff = v3_sub(points[i], points[j]);
	    
	    sum += weights[j] * 
		(fabs(Vx(diff)) + fabs(Vy(diff)) + fabs(Vz(diff)));
	}

	if (sum < best_sum) {
	    best_sum = sum;
	    best_idx = i;
	}
    }

    if (best_idx == -1)
	return v3_new(0.0, 0.0, 0.0);
    
    return points[best_idx];    
}

double FindRobustVariance(v3_t mean, const std::vector<v3_t> &points)
{
    std::vector<double> dists;
    int num_points = (int) points.size();
    
    for (int i = 0; i < num_points; i++) {
	v3_t disp = v3_sub(points[i], mean);
	double dist = v3_magsq(disp);
	
	dists.push_back(dist);
    }

    nth_element(dists.begin(), 
		dists.begin() + num_points / 2, 
		dists.end());

    return dists[num_points / 2];
}

v3_t FindRobustMean(const std::vector<v3_t> &points) 
{
    int num_points = (int) points.size();
    double best_sum = DBL_MAX;
    int best_idx = -1;

    for (int i = 0; i < num_points; i++) {
	double sum = 0.0;

	for (int j = 0; j < num_points; j++) {
	    v3_t diff = v3_sub(points[i], points[j]);
	    
	    sum += fabs(Vx(diff)) + fabs(Vy(diff)) + fabs(Vz(diff));
	}

	if (sum < best_sum) {
	    best_sum = sum;
	    best_idx = i;
	}
    }

    if (best_idx == -1)
	return v3_new(0.0, 0.0, 0.0);
    
    return points[best_idx];
}

/* Fit a plane to the points at the given indices */
std::vector<int> FitPlaneToPoints(const std::vector<PointData> &points,
                                  const std::vector<int> &indices,
                                  double *plane, 
                                  int ransac_rounds,
                                  double ransac_threshold,
                                  bool par_to_up, bool perp_to_up,
                                  double *up)
{
    if (par_to_up && perp_to_up) {
	printf("[FitPlaneToPoints] Error: cannot be both "
	       "parallel and perpendicular to the up vector!\n");
	perp_to_up = false;
    }

    std::vector<int> inliers;
    if (!par_to_up) {
        /* Marshall the points */
        int num_points = (int) indices.size();
        v3_t *pts = new v3_t[num_points];
        for (int i = 0; i < num_points; i++) {
            int pt_idx = indices[i];
            const PointData &pt = points[pt_idx];
            pts[i] = v3_new(pt.m_pos[0], pt.m_pos[1], pt.m_pos[2]);
        }

        /* Fit the plane */
        int num_inliers = 0;
        double error = fit_3D_plane_ortreg_ransac(num_points, pts, 
                                                  ransac_rounds, 
                                                  ransac_threshold,
                                                  &num_inliers, plane);

        printf("error = %0.3f\n", error);

        /* Gather the inliers */
        for (int i = 0; i < num_points; i++) {
            double dist = plane_point_distance(plane, pts[i]);
            if (dist < ransac_threshold) {
                inliers.push_back(indices[i]);
            }
        }

        if (perp_to_up) {
            /* Compute the mean of the inliers */
            int num_inliers = (int) inliers.size();
            v3_t *pts_inlier = new v3_t[num_inliers];

            for (int i = 0; i < num_inliers; i++) {
                int pt_idx = inliers[i];
                const PointData &pt = points[pt_idx];
                pts_inlier[i] = v3_new(pt.m_pos[0], pt.m_pos[1], pt.m_pos[2]);
            }
            
            v3_t mean = v3_mean(num_inliers, pts_inlier);

            double dot;
            matrix_product(1, 3, 3, 1, up, mean.p, &dot);

            plane[0] = up[0];
            plane[1] = up[1];
            plane[2] = up[2];
            plane[3] = -dot;

            delete [] pts_inlier;
        }

        delete [] pts;
    } else {
        assert(fabs(up[1] - 1.0) < 1.0e-5);

        /* Marshall the points */
        int num_points = (int) indices.size(); // points.size();
        v2_t *pts = new v2_t[num_points];
        for (int i = 0; i < num_points; i++) {
            int pt_idx = indices[i];
            const PointData &pt = points[pt_idx];

            pts[i] = v2_new(pt.m_pos[0], pt.m_pos[2]);
        }

        /* Fit the plane */
        double line[3];
        int num_inliers = 0;
        double error = fit_2D_line_ortreg_ransac(num_points, pts, 
                                                 ransac_rounds, 
                                                 ransac_threshold,
                                                 &num_inliers, line);

        plane[0] = line[0];
        plane[1] = 0.0;
        plane[2] = line[1];
        plane[3] = line[2];

        printf("error = %0.3f\n", error);
        printf("num_inliers = %d\n", num_inliers);

        /* Gather the inliers */
        for (int i = 0; i < num_points; i++) {
            int pt_idx = indices[i];
            const PointData &p = points[pt_idx];
            v3_t pt = v3_new(p.m_pos[0], p.m_pos[1], p.m_pos[2]);
            double dist = plane_point_distance(plane, pt);

            if (dist < ransac_threshold) {
                inliers.push_back(indices[i]);
            }
        }
    }

    return inliers;
}

std::vector<v2_t> GetPointProjections(const CameraInfo &cam, 
                                      const std::vector<PointData> &points,
                                      const std::vector<int> &indices,
                                      bool inside_only, 
                                      int &num_inside)
{
    int num_points = (int) indices.size();
    std::vector<v2_t> projs;
    BoundingBox bbox = cam.GetBoundingBox();

    num_inside = 0;

    for (int i = 0; i < num_points; i++) {
	int pidx = indices[i];
	const PointData &p = points[pidx];
	
	double proj[2];

#if 1
	bool in_front = cam.Project(p.m_pos, proj);

	if (!in_front)
	    continue;

	bool inside = bbox.Contains(proj[0], proj[1]);
#else
	bool inside = cam.Project(p.m_pos, proj);
#endif

	if (inside)
	    num_inside++;

	if (inside_only && inside) 
	    projs.push_back(v2_new(proj[0], proj[1]));
	else if (!inside_only)
	    projs.push_back(v2_new(proj[0], proj[1]));
    }

    return projs;
}
