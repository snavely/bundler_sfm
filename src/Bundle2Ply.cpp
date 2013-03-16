/* Bundle2Ply.cpp */
/* Convert a bundle file to a ply file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "matrix.h"
#include "sfm.h"

typedef struct 
{
    double pos[3];
    double color[3];
} point_t;

void ReadBundleFile(char *bundle_file, 
                    std::vector<camera_params_t> &cameras,
                    std::vector<point_t> &points, double &bundle_version)
{
    FILE *f = fopen(bundle_file, "r");
    if (f == NULL) {
	printf("Error opening file %s for reading\n", bundle_file);
	return;
    }

    int num_images, num_points;
    bool coalesced;

    char first_line[256];
    fgets(first_line, 256, f);
    if (first_line[0] == '#') {
        double version;
        sscanf(first_line, "# Bundle file v%lf", &version);

        bundle_version = version;
        printf("[ReadBundleFile] Bundle version: %0.3f\n", version);

        fscanf(f, "%d %d\n", &num_images, &num_points);
    } else if (first_line[0] == 'v') {
        double version;
        sscanf(first_line, "v%lf", &version);
        bundle_version = version;
        printf("[ReadBundleFile] Bundle version: %0.3f\n", version);

        fscanf(f, "%d %d\n", &num_images, &num_points);
    } else {
        bundle_version = 0.1;
        sscanf(first_line, "%d %d\n", &num_images, &num_points);
    }

    printf("[ReadBundleFile] Reading %d images and %d points...\n",
	   num_images, num_points);

    /* Read cameras */
    for (int i = 0; i < num_images; i++) {
	double focal_length, k0, k1;
	double R[9];
	double t[3];
        
        if (bundle_version < 0.2) {
            /* Focal length */
            fscanf(f, "%lf\n", &focal_length);
        } else {
            fscanf(f, "%lf %lf %lf\n", &focal_length, &k0, &k1);
        }

	/* Rotation */
	fscanf(f, "%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n", 
	       R+0, R+1, R+2, R+3, R+4, R+5, R+6, R+7, R+8);
	/* Translation */
	fscanf(f, "%lf %lf %lf\n", t+0, t+1, t+2);

        if (focal_length == 0.0)
            continue;

        camera_params_t cam;

        cam.f = focal_length;
        memcpy(cam.R, R, sizeof(double) * 9);
        memcpy(cam.t, t, sizeof(double) * 3);

        cameras.push_back(cam);
    }
    
    /* Read points */
    for (int i = 0; i < num_points; i++) {
	point_t pt;

	/* Position */
	fscanf(f, "%lf %lf %lf\n", 
	       pt.pos + 0, pt.pos + 1, pt.pos + 2);

	/* Color */
	fscanf(f, "%lf %lf %lf\n", 
	       pt.color + 0, pt.color + 1, pt.color + 2);

	int num_visible;
	fscanf(f, "%d", &num_visible);

	for (int j = 0; j < num_visible; j++) {
	    int view, key;
	    fscanf(f, "%d %d", &view, &key);

            double x, y;
            if (bundle_version >= 0.3)
                fscanf(f, "%lf %lf", &x, &y);
	}

        // if (num_visible > 0) {
        points.push_back(pt);
        // }
    }

    fclose(f);
}

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
void WritePlyFile(char *ply_file,
                  std::vector<camera_params_t> &cameras,
                  std::vector<point_t> &points, double bundle_version,
                  int decimate_step, int write_cameras)
{
    FILE *f = fopen(ply_file, "w");

    if (f == NULL) {
	printf("Error opening file %s for writing\n", ply_file);
	return;
    }

    int num_cameras = (int) cameras.size();
    int num_points = (int) points.size();

    int num_points_decimate = 
        num_points / decimate_step; // + ((num_points % decimate_step) == 0);
    
    int num_points_out = 0;
    if (write_cameras)
        num_points_out = num_points_decimate + 2 * num_cameras;
    else
        num_points_out = num_points_decimate;

    /* Print the ply header */
    fprintf(f, ply_header, num_points_out);

    bool reflect = bundle_version < 0.3;

    /* Now triangulate all the correspondences */
    for (int i = 0; i < num_points; i += decimate_step) {
	/* Output the vertex */
	fprintf(f, "%0.16e %0.16e %0.16e %d %d %d\n", 
		points[i].pos[0], points[i].pos[1], 
                (reflect ? -1.0 : 1.0) * points[i].pos[2],
                (int) points[i].color[0], 
                (int) points[i].color[1], 
                (int) points[i].color[2]);
    }

    if (write_cameras) {
        for (int i = 0; i < num_cameras; i++) {
            double c[3];

            double Rinv[9];
            matrix_transpose(3, 3, cameras[i].R, Rinv);

            matrix_product(3, 3, 3, 1, Rinv, cameras[i].t, c);
            matrix_scale(3, 1, c, -1.0, c);
	
            if (reflect)
                c[2] = -c[2];

            if ((i % 2) == 0)
                fprintf(f, "%0.16e %0.16e %0.16e 0 255 0\n", c[0], c[1], c[2]);
            else
                fprintf(f, "%0.16e %0.16e %0.16e 255 0 0\n", c[0], c[1], c[2]);

            double p_cam[3] = { 0.0, 0.0, -0.3 };
            double p[3];

            if (reflect)
                p_cam[2] = -p_cam[2];

            matrix_product(3, 3, 3, 1, Rinv, p_cam, p);

            p[0] += c[0];
            p[1] += c[1];
            p[2] += c[2];

            fprintf(f, "%0.16e %0.16e %0.16e 255 255 0\n",
                    p[0], p[1], p[2]);
        }
    }

    fclose(f);
}


int main(int argc, char **argv) 
{
    if (argc != 3 && argc != 4 && argc != 5) {
        printf("Usage: %s <bundle.out> <out.ply> "
               "[decimate:1] [write_cameras:1]\n", argv[0]);
        return 1;
    }
    
    char *bundle_file = argv[1];
    char *ply_file = argv[2];

    int decimate_step = 1;
    if (argc >= 4) {
        decimate_step = atoi(argv[3]);
    }

    int write_cameras = 1;
    if (argc >= 5) {
        write_cameras = atoi(argv[4]);
    }
    
    /* Read the bundle file */
    std::vector<camera_params_t> cameras;
    std::vector<point_t> points;
    double bundle_version;
    ReadBundleFile(bundle_file, cameras, points, bundle_version);
    WritePlyFile(ply_file, cameras, points, 
                 bundle_version, decimate_step, write_cameras);

    return 0;
}
