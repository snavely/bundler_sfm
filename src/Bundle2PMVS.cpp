/* Bundle2PMVS.cpp */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>

/* For mkdir */
#include <sys/stat.h>
#include <sys/types.h>

#include "image.h"
#include "matrix.h"
#include "sfm.h"
#include "LoadJPEG.h"

typedef struct 
{
    double pos[3];
    double color[3];
    
} point_t;

void ReadListFile(const char *list_file, std::vector<std::string> &files)
{
    FILE *f = fopen(list_file, "r");
    
    char buf[256];
    while (fgets(buf, 256, f)) {
        if (buf[strlen(buf)-1] == '\n')
            buf[strlen(buf)-1] = 0;

        char *space = strchr(buf, ' ');
        if (space) *space = 0;

        files.push_back(std::string(buf));
    }

    fclose(f);
}

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

        // if (focal_length == 0.0)
        //     continue;

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

        if (num_visible > 0) {
            points.push_back(pt);
        }
    }

    fclose(f);
}

void WritePMVS(const char *output_path, 
               const char *list_file, const char *bundle_file,
               std::vector<std::string> images, 
               std::vector<camera_params_t> &cameras)
{
    int num_cameras = (int) cameras.size();

    /* Make sure output_path exists */
    mkdir(output_path, 0770);

    char buf[2048];
    sprintf(buf, "%s/prep_pmvs.sh", output_path);

    FILE *f_scr = fopen(buf, "w");

    fprintf(f_scr, "# Script for preparing images and calibration data \n"
            "#   for Yasutaka Furukawa's PMVS system\n\n");
    fprintf(f_scr, "BUNDLER_BIN_PATH="" # Edit this line before running\n");
    fprintf(f_scr, "if [ \"$BUNDLER_BIN_PATH\" == \"\" ] ; then "
            "echo Please edit prep_pmvs.sh to specify the path to the "
            " bundler binaries.; exit; fi\n");

    fprintf(f_scr, "# Apply radial undistortion to the images\n");
    fprintf(f_scr, "$BUNDLER_BIN_PATH/RadialUndistort %s %s %s\n", 
            list_file, bundle_file, output_path);
    fprintf(f_scr, "\n# Create directory structure\n");
    // fprintf(f_scr, "mkdir -p %s/\n", output_path);
    fprintf(f_scr, "mkdir -p %s/txt/\n", output_path);
    fprintf(f_scr, "mkdir -p %s/visualize/\n", output_path);
    fprintf(f_scr, "mkdir -p %s/models/\n", output_path);
    fprintf(f_scr, "\n# Copy and rename files\n");

    int count = 0;
    for (int i = 0; i < num_cameras; i++) {
        if (cameras[i].f == 0.0)
            continue;

        char buf[256];
        sprintf(buf, "%s/%08d.txt", output_path, count);
        FILE *f = fopen(buf, "w");
        assert(f);

        /* Compute the projection matrix */
        double focal = cameras[i].f;
        double *R = cameras[i].R;
        double *t = cameras[i].t;

        int w, h;
        GetJPEGDimensions(images[i].c_str(), w, h);

        double K[9] = 
            { -focal, 0.0, 0.5 * w - 0.5,
              0.0, focal, 0.5 * h - 0.5,
              0.0, 0.0, 1.0 };

        double Ptmp[12] = 
            { R[0], R[1], R[2], t[0],
              R[3], R[4], R[5], t[1],
              R[6], R[7], R[8], t[2] };
        
        double P[12];
        matrix_product(3, 3, 3, 4, K, Ptmp, P);
        matrix_scale(3, 4, P, -1.0, P);

        fprintf(f, "CONTOUR\n");
        fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[0], P[1], P[2],  P[3]);
        fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[4], P[5], P[6],  P[7]);
        fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[8], P[9], P[10], P[11]);

        fclose(f);

        int last_dot = images[i].rfind('.', images[i].length()-1);
        int last_slash = images[i].rfind('/', last_dot-1);
        std::string basename = images[i].substr(last_slash + 1, last_dot - last_slash - 1);

        fprintf(f_scr, "mv pmvs/%s.rd.jpg %s/visualize/%08d.jpg\n", 
                basename.c_str(), output_path, count);
        fprintf(f_scr, "mv %s %s/txt/\n", buf, output_path);

        count++;
    }

    /* Write the options file */
    sprintf(buf, "%s/pmvs_options.txt", output_path);
    FILE *f_opt = fopen(buf, "w");

    fprintf(f_opt, "level 1\n");
    fprintf(f_opt, "csize 2\n");
    fprintf(f_opt, "threshold 0.7\n");
    fprintf(f_opt, "wsize 7\n");
    fprintf(f_opt, "minImageNum 3\n");
    fprintf(f_opt, "CPU 8\n");
    fprintf(f_opt, "setEdge 0\n");
    fprintf(f_opt, "useBound 0\n");
    fprintf(f_opt, "useVisData 1\n");
    fprintf(f_opt, "sequence -1\n");
    fprintf(f_opt, "timages -1 0 %d\n", count);
    fprintf(f_opt, "oimages -3\n");

    fclose(f_opt);

    fprintf(f_scr, "\necho \"Running Bundle2Vis to generate vis.dat\n\"\n");
    fprintf(f_scr, "$BUNDLER_BIN_PATH/Bundle2Vis %s/bundle.rd.out %s/vis.dat\n", 
            output_path, output_path);

    fprintf(f_scr, "\n\n\necho @@ Sample command for running pmvs:\n");
    fprintf(f_scr, "echo \"   pmvs2 %s/ pmvs_options.txt\"\n", output_path);
    fprintf(f_scr, "echo \"    - or - \"\n");
    fprintf(f_scr, "echo \"   use Dr. Yasutaka Furukawa's view clustering algorithm to generate a set of options files.\"\n");
    fprintf(f_scr, "echo \"       The clustering software is available at http://grail.cs.washington.edu/software/cmvs\"\n");

    fclose(f_scr);
}

int main(int argc, char **argv) 
{
    if (argc != 3 && argc != 4) {
        printf("Usage: %s <list.txt> <bundle.out> [pmvs_output_path (default: pmvs)]\n", 
               argv[0]);
        return 1;
    }
    
    const char *list_file = argv[1];
    const char *bundle_file = argv[2];
    const char *output_path = "pmvs";

    if (argc == 4)
        output_path = argv[3];

    /* Read the list file */
    FILE *f = fopen(list_file, "r");
    if (f == NULL) {
        printf("Error opening file %s for reading\n", list_file);
        return 1;
    }

    std::vector<std::string> images;
    ReadListFile(list_file, images);

    /* Read the bundle file */
    std::vector<camera_params_t> cameras;
    std::vector<point_t> points;
    double bundle_version;
    ReadBundleFile(bundle_file, cameras, points, bundle_version);

    /* Write camera geometry in the PMVS file format */
    WritePMVS(output_path, list_file, bundle_file, images, cameras);

    printf("\n\n");
    printf("@@ Conversion complete, execute \"sh %s/prep_pmvs.sh\" to finalize\n", output_path);
    printf("@@ (you will first need to edit prep_pmvs.sh to specify your bundler path, \n");
    printf("@@  so that the script knows where to find your\n"
           "@@  RadialUndistort and Bundle2Vis binaries)\n");

    return 0;
}
