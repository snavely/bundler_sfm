/* FisheyeUndistort.cpp */

#include <string>
#include <vector>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "image.h"
#include "matrix.h"
#include "resample.h"
#include "util.h"

#include "LoadJPEG.h"

#include "BundleUtil.h"

typedef struct fisheye_params_t {
    fisheye_params_t() : 
        fCx(0.0), fCy(0.0), fRad(0.0), fAngle(0.0), fFocal(0.0)
    {}
        
    double fCx;
    double fCy;
    double fRad;
    double fAngle;
    double fFocal;
} fisheye_params_t;

static void ReadFisheyeParameters(char *filename, 
                                  fisheye_params_t &fisheye_params)
{
    fisheye_params.fCx = 0.0;
    fisheye_params.fCy = 0.0;

    /* Read fisheye params */
    FILE *f = fopen(filename, "r");

    if (f == NULL) {
	printf("Error opening fisheye parameters file %s for reading\n", 
	       filename);
	exit(1);
    }

    char buf[256];
    while (fgets(buf, 256, f) != NULL) {
	/* Split the command into tokens */
        std::string str(buf);    
	std::vector<std::string> toks;
	
        Tokenize(str, toks, " ");

	if (toks[0] == "FisheyeCenter:") {
	    if (toks.size() != 3) {
		printf("FisheyeCenter needs two arguments!\n");
		exit(1);
	    } else {
                fisheye_params.fCx = atof(toks[1].c_str());
                fisheye_params.fCy = atof(toks[2].c_str());
	    }
	} else if (toks[0] == "FisheyeRadius:") {
	    if (toks.size() != 2) {
		printf("FisheyeRadius needs one argument!\n");
		exit(1);
	    } else {
                fisheye_params.fRad = atof(toks[1].c_str());
	    }
	} else if (toks[0] == "FisheyeAngle:") {
	    if (toks.size() != 2) {
		printf("FisheyeAngle needs one argument!\n");
		exit(1);
	    } else {
                fisheye_params.fAngle = atof(toks[1].c_str());
	    }	    
	} else if (toks[0] == "FisheyeFocal:") {
	    if (toks.size() != 2) {
		printf("FisheyeFocal needs one argument!\n");
		exit(1);
	    } else {
                fisheye_params.fFocal = atof(toks[1].c_str());
	    }
	} else {
	    printf("Unrecognized fisheye field %s\n", toks[0].c_str());
	}
    }

    fclose(f);
}

void DistortPoint(const double x, const double y, 
		  const fisheye_params_t &p, double *R,
                  double &x_out, double &y_out)
{
    double xn = x; // - m_fCx;
    double yn = y; // - m_fCy;

    double ray[3] = { xn, yn, p.fFocal }, ray_rot[3];
    matrix_product(3, 3, 3, 1, R, ray, ray_rot);

    if (ray_rot[2] <= 0.0) {
	xn = -DBL_MAX;
	yn = -DBL_MAX;
	return;
    } else {
	xn = ray_rot[0] * p.fFocal / ray_rot[2];
	yn = ray_rot[1] * p.fFocal / ray_rot[2];
    }
    
    double r = sqrt(xn * xn + yn * yn);
    double angle = RAD2DEG(atan(r / p.fFocal));
    double rnew = p.fRad * angle / (0.5 * p.fAngle);
    
    x_out = xn * (rnew / r) + p.fCx;
    y_out = yn * (rnew / r) + p.fCy;
}

img_t *UndistortImage(img_t *img, 
		      const fisheye_params_t &fisheye_params) {
    double R[9];
    matrix_ident(3, R);

    int w = img->w;
    int h = img->h;

    img_t *img_out = img_new(w, h);

    for (int y = 0; y < h; y++) {
	for (int x = 0; x < w; x++) {
	    double xn = x - 0.5 * w;
	    double yn = y - 0.5 * h;
	    double x_new = 0.0, y_new = 0.0;
	    
	    DistortPoint(xn, yn, fisheye_params, R, x_new, y_new);

	    x_new += 0.5 * w;
	    y_new += 0.5 * h;

	    if (x_new < 0 || x_new >= w || y_new < 0 || y_new >= h)
		continue;

	    fcolor_t c = pixel_lerp(img, x_new, y_new);

	    img_set_pixel(img_out, x, y, 
			  iround(c.r), iround(c.g), iround(c.b));
	}
    }

    return img_out;
}


void UndistortImages(const std::vector<image_t> &images, 
                     const fisheye_params_t &fisheye_params)
{
    int num_images = (int) images.size();

    for (int i = 0; i < num_images; i++) {
        img_t *img = LoadJPEG(images[i].name.c_str());
        img_t *img_u;

        if (images[i].is_fisheye) {
	  printf("Undistorting image %s\n", images[i].name.c_str());
	  img_u = UndistortImage(img, fisheye_params);
	} else {
	  printf("Skipping image %s (not marked as fisheye).\n", 
		 images[i].name.c_str());
	  img_u = img_copy(img);
	}

        const std::string out = 
            images[i].name.substr(0, 
                    images[i].name.length() - 3).append("fd.jpg");
        WriteJPEG(img_u, (const char *) out.c_str());

        img_free(img);
        img_free(img_u);
    }
}

int main(int argc, char **argv) 
{
    if (argc != 3) {
        printf("Usage: %s <list.txt> <fisheye.txt>\n", argv[0]);
        return 1;
    }
    
    char *list_file = argv[1];
    char *fisheye_file = argv[2];

    std::vector<image_t> images;

    fisheye_params_t fisheye_params;
    ReadFisheyeParameters(fisheye_file, fisheye_params);
    ReadListFile(list_file, images);
    UndistortImages(images, fisheye_params);

    return 0;
}
