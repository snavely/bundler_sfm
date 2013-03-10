/* Bundle2Vis.cpp */
/* Convert a bundle file to a PMVS vis.dat file */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "sfm.h"

typedef std::pair<int,int> ImageKey;

typedef struct 
{
    double pos[3];
    double color[3];
    std::vector<ImageKey> views;
} point_t;

void ReadBundleFile(const char *bundle_file, 
                    std::vector<camera_params_t> &cameras,
                    std::vector<point_t> &points, double &bundle_version)
{
    FILE *f = fopen(bundle_file, "r");
    if (f == NULL) {
	printf("Error opening file %s for reading\n", bundle_file);
	return;
    }

    int num_images, num_points;

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

    // int *map = new int[num_images];

    /* Read cameras */
    // int count = 0;
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

        // if (focal_length == 0.0)
        //     continue;

        camera_params_t cam;

        cam.f = focal_length;
        memcpy(cam.R, R, sizeof(double) * 9);
        memcpy(cam.t, t, sizeof(double) * 3);

        // map[i] = count;
        // count++;

        cameras.push_back(cam);
    }
    
    /* Read points */
    int total_num_visible = 0;
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
        total_num_visible += num_visible;

	for (int j = 0; j < num_visible; j++) {
	    int view, key;
	    fscanf(f, "%d %d", &view, &key);

            // pt.views.push_back(ImageKey(map[view],key));
            assert(view >= 0 && view < num_images);
            pt.views.push_back(ImageKey(view, key));

            double x, y;
            if (bundle_version >= 0.3)
                fscanf(f, "%lf %lf", &x, &y);
	}

        if (num_visible > 0) {
            points.push_back(pt);
        }
    }

    printf("Num visible: %d\n", total_num_visible);

    fclose(f);
}

void WriteVisFile(const char *vis_file, 
                  std::vector<camera_params_t> &cameras,
                  std::vector<point_t> &points)
{
    int nCameras = (int) cameras.size();
    int nPoints = (int) points.size();

    printf("Num cameras: %d\n", nCameras);

    FILE *f = fopen(vis_file, "w");
    if (!f) {
        printf("[WriteVisFile] Error opening "
               "file %s for writing\n", vis_file);
        return;
    }

    /* Fill in the matches matrix */
    unsigned int *matches = new unsigned int[nCameras * nCameras];
    for (int i = 0; i < nCameras; i++) {
        for (int j = 0; j < nCameras; j++) {
            matches[i * nCameras + j] = 0;
        }
    }
    
    for (int i = 0; i < nPoints; i++) {
        // bool seen = false;
        int nViews = (int) points[i].views.size();
        for (int j = 0; j < nViews; j++) {
            int i1 = points[i].views[j].first;
            for (int k = j+1; k < nViews; k++) {
                if (j == k) continue;
                int i2 = points[i].views[k].first;

                matches[i1 * nCameras + i2]++;
                matches[i2 * nCameras + i1]++;
            }
        }
    }

    fprintf(f, "VISDATA\n");
    fprintf(f, "%d\n", nCameras);    

    // write camera rows
    const unsigned int MATCH_THRESHOLD = 32;
    for (int i = 0; i < nCameras; i++) {
        std::vector<int> vis;
        for (int j = 0; j < nCameras; j++) {
            if (matches[i * nCameras + j] >= MATCH_THRESHOLD)
                vis.push_back(j);
        }
         
        int nVis = (int) vis.size();
        fprintf(f, "%d %d", i, (int) nVis);

        for (int j = 0; j < nVis; j++) {
            fprintf(f, " %d", vis[j]);
        }

        fprintf(f, "\n");
    }

    fclose(f);
}

int main(int argc, char **argv) 
{
    if (argc != 3) {
        printf("Usage: %s <bundle.out> <vis.dat>\n", argv[0]);
        return 1;
    }
    
    char *bundle_file = argv[1];
    char *vis_file = argv[2];

    /* Read the bundle file */
    std::vector<camera_params_t> cameras;
    std::vector<point_t> points;
    double bundle_version;
    ReadBundleFile(bundle_file, cameras, points, bundle_version);
    WriteVisFile(vis_file, cameras, points);

    return 0;
}
