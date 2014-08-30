/* BundleCeres.cpp */
/* Bundle adjust using Ceres */

#include <algorithm>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>

#include "matrix.h"
#include "qsort.h"
#include "util.h"

#include "Bundle.h"
#include "BundlerApp.h"
#include "Distortion.h"

#include <ceres/ceres.h>
#include "snavely_reprojection_error.h"

#include <omp.h>

// #define USE_L2_NORM

static int compare_doubles(const void *d1, const void *d2)
{
    double a = *(double *) d1;
    double b = *(double *) d2;
    
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

static void aa2rot(double const * x, double * R)
{
    // go from a vector representing a rotation in 
    // axis-angle format to a 3x3 rotation matrix 
    // in column major form
  
    const double epsilon = 1e-18;
    double theta = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    // dblas::norm2(3,const_cast<double *>(x));
    double wx = x[0]/(theta+epsilon);
    double wy = x[1]/(theta+epsilon);
    double wz = x[2]/(theta+epsilon);
    
    double costheta = cos(theta);
    double sintheta = sin(theta);
  
    R[0] = costheta + wx*wx*(1-costheta);
    R[1] = wz*sintheta + wx*wy*(1-costheta);
    R[2] = -wy*sintheta + wx*wz*(1-costheta);
    R[3] = wx*wy*(1-costheta) - wz*sintheta;
    R[4] = costheta + wy*wy*(1-costheta);
    R[5] = wx*sintheta + wy*wz*(1-costheta);
    R[6] = wy*sintheta + wx*wz*(1-costheta);
    R[7] = -wx*sintheta + wy*wz*(1-costheta);
    R[8] = costheta + wz*wz*(1-costheta);
};

// rotation to axis angle vector format
static void rot2aa(double *R, double *aa)
{
    double tr = R[0] + R[4] + R[8];
    double costheta = CLAMP(0.5 * (tr - 1.0), -1, 1);

    double RT[9], RRT[9];
    matrix_transpose(3, 3, R, RT);
    matrix_diff(3, 3, 3, 3, R, RT, RRT);
    double sintheta = matrix_norm(3, 3, RRT)/sqrt(8.0);

    double theta = atan2(sintheta,costheta);
    double factor = theta / (2.0 * sintheta + 1.0e-10);

    aa[0] = factor * (R[7]-R[5]);
    aa[1] = factor * (R[2]-R[6]);
    aa[2] = factor * (R[3]-R[1]);
}

/* Penalize a camera variable for deviating from a given prior value */
struct PriorError {
    PriorError(int prior_index, double prior_value, double prior_scale) : 
        prior_index(prior_index), prior_value(prior_value), prior_scale(prior_scale) {}

    template <typename T>
    bool operator()(const T * const x, T * residual) const {
        residual[0] = prior_scale * (prior_value - x[prior_index]);
        return true;
    }

    int prior_index;
    double prior_value;
    double prior_scale;
};

double BundlerApp::RunSFM_Ceres(int num_pts, int num_cameras, 
                                int start_camera,
                                bool fix_points, 
                                camera_params_t *init_camera_params,
                                v3_t *init_pts, int *added_order, 
                                v3_t *colors,
                                std::vector<ImageKeyVector> &pt_views, 
                                int max_iter, int max_iter2, 
                                int verbosity, double eps2, 
                                double *S, double *U, double *V, double *W,
                                bool remove_outliers, bool final_bundle, 
                                bool write_intermediate)
{
    CheckPointKeyConsistency(pt_views, added_order);

#define MIN_POINTS 20
    int num_outliers = 0;
    int total_outliers = 0;
    double dist_total = 0.0;
    int num_dists = 0;

    int *remap = new int [num_pts];
    v3_t *nz_pts = new v3_t[num_pts];

    int round = 0;
    const int MIN_OUTLIERS = 40; /* was 10 */
    const double HUBER_PARAM = 25.0; /*16.0;*/ /* was 4.0 */

    timeval start_all, stop_all;
    gettimeofday(&start_all, NULL);

    assert(m_estimate_distortion);

    bool use_dense = true;
    if (num_cameras > 200)
        use_dense = false;

    do {
        ceres::Problem problem;

        timeval start_iter, stop_iter;
        gettimeofday(&start_iter, NULL);

        if (num_pts - total_outliers < MIN_POINTS) {
            printf("[RunSFM] Too few points remaining, exiting!\n");
            fflush(stdout);

            dist_total = DBL_MAX;
            break;
        }

	/* Set up the vmask and projections */
	// char *vmask = NULL;
	double *projections = NULL;
        int *pidx = NULL;
        int *cidx = NULL;
        unsigned int *num_vis = NULL;
	
	int num_projections = 0;
	for (int i = 0; i < num_pts; i++) {
	    num_projections += (int) pt_views[i].size();
	}

	// vmask = new char[num_pts * num_cameras];
	projections = new double[2 * num_projections];
        pidx = new int[num_projections];
        cidx = new int[num_projections];
        num_vis = new unsigned int[num_cameras];

	// for (int i = 0; i < num_pts * num_cameras; i++)
        //    vmask[i] = 0;

        for (int i = 0; i < num_cameras; i++)
            num_vis[i] = 0;

	int arr_idx = 0;
        int nz_count = 0;
	for (int i = 0; i < num_pts; i++) {
	    int num_views = (int) pt_views[i].size();

            if (num_views > 0) {
                for (int j = 0; j < num_views; j++) {
                    int c = pt_views[i][j].first;
                    int v = added_order[c];
                    int k = pt_views[i][j].second;

                    // vmask[nz_count * num_cameras + c] = 1;

                    projections[2 * arr_idx + 0] = GetKey(v,k).m_x;
                    projections[2 * arr_idx + 1] = GetKey(v,k).m_y;

                    pidx[arr_idx] = nz_count;
                    cidx[arr_idx] = c;

                    num_vis[c]++;

                    arr_idx++;
                }

                remap[i] = nz_count;
                nz_pts[nz_count] = init_pts[i];
                nz_count++;
            } else {
                remap[i] = -1;
            }
	}

        /* Set up initial parameters */
        int num_nz_points = nz_count;
        double *init_x = NULL;

        int cnp = 0;
        if (!m_estimate_distortion)
            cnp = 7;
        else 
            cnp = 9;

        int pnp = 3;

        unsigned int num_parameters = 
            pnp * num_nz_points + cnp * num_cameras;

        init_x = new double[num_parameters];

        double *cameras = init_x;
        double *points = init_x + cnp * num_cameras;

        /* Insert camera parameters */
        double *ptr = init_x;

        int idx = 0;
        for (int i = 0; i < num_cameras; i++) {
            /* Get the rotation */
            double axis[3]; //, angle;
            // double RT[9];
            // matrix_transpose(3, 3, init_camera_params[i].R, RT);
            rot2aa(init_camera_params[i].R, axis); // , &angle);
            // matrix_scale(3, 1, axis, angle, axis);

            double *c = init_camera_params[i].t;
            double t[3];
            matrix_product331(init_camera_params[i].R, c, t);
            matrix_scale(3, 1, t, -1.0, t);

            double f = init_camera_params[i].f;
            double k1 = init_camera_params[i].k[0];
            double k2 = init_camera_params[i].k[1];

            ptr[idx] = axis[0]; idx++;
            ptr[idx] = axis[1]; idx++;
            ptr[idx] = axis[2]; idx++;
            ptr[idx] = t[0]; idx++;
            ptr[idx] = t[1]; idx++;
            ptr[idx] = t[2]; idx++;
            ptr[idx] = f; idx++;
            ptr[idx] = k1; idx++; // / (f * f); idx++;
            ptr[idx] = k2; idx++; // / (f * f * f * f); idx++;
        }
        
        /* Insert point parameters */
	for (int i = 0; i < num_pts; i++) {
	    int num_views = (int) pt_views[i].size();

            if (num_views > 0) {
                ptr[idx] = Vx(init_pts[i]); idx++;
                ptr[idx] = Vy(init_pts[i]); idx++;
                ptr[idx] = Vz(init_pts[i]); idx++;
            }
        }
        
        assert(idx == (int) num_parameters);

        for (int i = 0; i < num_projections; ++i) {
            ceres::CostFunction *cost_function;

            /* Each Residual block takes a point and a camera as input and
             * outputs a 2 dimensional residual. */
            cost_function =
                new ceres::AutoDiffCostFunction<ceres::examples::SnavelyReprojectionError, 2, 9, 3>(new ceres::examples::SnavelyReprojectionError(projections[2 * i + 0],
                                                                                                                              projections[2 * i + 1]));

            // If enabled use Huber's loss function.
            ceres::LossFunction *loss_function = NULL;

#ifndef USE_L2_NORM
            loss_function = new ceres::HuberLoss(HUBER_PARAM);
#endif

            // Each observation correponds to a pair of a camera and a point
            // which are identified by camera_index()[i] and point_index()[i]
            // respectively.
            double *camera = cameras + cnp * cidx[i];
            double *point = points + pnp * pidx[i];

            problem.
                AddResidualBlock(cost_function, loss_function, camera, point);
        }

#if 1
        /* Now add the priors */
        for (int i = 0; i < num_cameras; i++) {
            ceres::CostFunction *prior_cost_function;
            if (init_camera_params[i].constrained[6]) {
                prior_cost_function =
                    new ceres::AutoDiffCostFunction<PriorError, 1, 9>(
                        new PriorError(6, init_camera_params[i].constraints[6],
                                       init_camera_params[i].weights[6] * num_vis[i]));

                problem.AddResidualBlock(prior_cost_function, NULL,
                                         cameras + cnp * i);
            }

            /* Take care of priors on distortion parameters */
            prior_cost_function = 
                new ceres::AutoDiffCostFunction<PriorError, 1, 9>(
                    new PriorError(7, 0.0, 0.0001 * m_distortion_weight * num_vis[i]));

            problem.AddResidualBlock(prior_cost_function, NULL,
                                     cameras + cnp * i);

            prior_cost_function = 
                new ceres::AutoDiffCostFunction<PriorError, 1, 9>(
                    new PriorError(8, 0.0, 0.0001 * m_distortion_weight * num_vis[i]));

            problem.AddResidualBlock(prior_cost_function, NULL,
                                     cameras + cnp * i);
        }
#endif

	dist_total = 0.0;
	num_dists = 0;

        clock_t start = clock();

        /* Make call to Ceres */
        printf("[Ceres] num_points = %d\n", num_nz_points);
        printf("[Ceres] num_cameras = %d\n", num_cameras);
        printf("[Ceres] num_params = %d\n", num_parameters);
        printf("[Ceres] num_projections = %d\n", num_projections);
        fflush(stdout);
        fflush(stderr);

#if 0
        if (true /*round < 2 == 0*/) {
            /* Write to file */
            char problem_file[256];
            sprintf(problem_file, "problem-%d-%d-pre%d.txt", 
                    num_cameras, num_nz_points, round);
            FILE *f = fopen(problem_file, "w");
            fprintf(f, "%d %d %d\n", 
                    num_cameras, num_nz_points, num_projections);

            /* write the projections */
            for (int i = 0; i < num_projections; i++) {
                fprintf(f, "%d %d     %0.6e %0.6e\n", 
                        cidx[i], pidx[i], 
                        projections[2*i], projections[2*i + 1]);
            }

            /* write init_x */
            for (unsigned int i = 0; i < num_parameters; i++) {
                fprintf(f, "%0.16e\n", init_x[i]);
            }

            fclose(f);
        }
#endif
        
        double *final_x = NULL;

        ceres::Solver::Options options;

        /* Select solver */
        bool iterative = false;
        if (use_dense && !final_bundle) {
            printf("[Ceres] Using CHOLMOD\n");
            options.linear_solver_type = ceres::SPARSE_SCHUR;
        } else {
            printf("[Ceres] Using CG\n");
            options.linear_solver_type = ceres::ITERATIVE_SCHUR;
            options.preconditioner_type = ceres::SCHUR_JACOBI;
            iterative = true;
        }

        options.parameter_tolerance = 1e-8;
        options.function_tolerance = 1e-4; // 1e1;

        if (max_iter == 0) {
            options.max_num_iterations = ((final_bundle) ? 60 : 25);

            if (round > 0) {
                if (!iterative) 
                    options.max_num_iterations = 10;
                else
                    options.max_num_iterations = 20;
            }
        } else {
            if (round == 0 || max_iter2 == 0)
                options.max_num_iterations = max_iter;
            else
                options.max_num_iterations = max_iter2;
        }
             
        // options.tau = (round > 0) ? 1e-8 : 1e-3; /* was 1e-5 */
        options.minimizer_progress_to_stdout = true;

        if (iterative) {
            /* seabad.precon_type = 1; 0: implicit-jacobi
                                       1: implicit-ssor
                                       2: implicit-none */

            // KNS: trying defaults...
            // options.linear_solver_min_num_iterations = 10;
            // options.linear_solver_max_num_iterations = 100;
            // options.eta = 0.1;
        }

        options.num_threads = omp_get_max_threads();

        options.num_linear_solver_threads = omp_get_max_threads();

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.FullReport() << "\n";
        fflush(stdout);
        fflush(stderr);
                
        if (summary.preprocessor_time_in_seconds >= 20) {
            use_dense = false;
        }

#if 0
        if (round == 0) {
            /* Write to file */
            char problem_file[256];
            sprintf(problem_file, "problem-%d-%d-post%d.txt", 
                    num_cameras, num_nz_points, round);
            FILE *f = fopen(problem_file, "w");
            fprintf(f, "%d %d %d\n", 
                    num_cameras, num_nz_points, num_projections);

            /* write the projections */
            for (int i = 0; i < num_projections; i++) {
                fprintf(f, "%d %d     %0.6e %0.6e\n", 
                        cidx[i], pidx[i], 
                        projections[2*i], projections[2*i + 1]);
            }

            /* write init_x */
            for (int i = 0; i < num_parameters; i++) {
                fprintf(f, "%0.16e\n", init_x[i]);
            }

            fclose(f);
        }
#endif
        
        final_x = init_x;
        ptr = final_x;
        for (int i = 0; i < num_cameras; i++) {
            /* Get the rotation */
            double axis[3], angle;

            axis[0] = *ptr; ptr++;
            axis[1] = *ptr; ptr++;
            axis[2] = *ptr; ptr++;

            // angle = matrix_norm(3, 1, axis);
            // matrix_scale(3, 1, axis, 1.0 / angle, axis);

            // axis_angle_to_matrix(axis, angle, init_camera_params[i].R);
            double RT[9];
            aa2rot(axis, RT);
            matrix_transpose(3, 3, RT, init_camera_params[i].R);

            double t[3];
            double *c = init_camera_params[i].t;

            t[0] = *ptr; ptr++;
            t[1] = *ptr; ptr++;
            t[2] = *ptr; ptr++;

            matrix_transpose_product(3, 3, 3, 1, 
                                     init_camera_params[i].R, t, c);
            matrix_scale(3, 1, c, -1.0, c);

            double f = init_camera_params[i].f = *ptr; ptr++;

            if (m_estimate_distortion) {
                init_camera_params[i].k[0] = *ptr; ptr++;
                init_camera_params[i].k[1] = *ptr; ptr++;
                // init_camera_params[i].k[0] = *ptr * (f * f); ptr++;
                // init_camera_params[i].k[1] = *ptr * (f * f * f * f); ptr++;
            }
        }

        /* Insert point parameters */
	for (int i = 0; i < num_nz_points; i++) {
            Vx(nz_pts[i]) = *ptr; ptr++;
            Vy(nz_pts[i]) = *ptr; ptr++;
            Vz(nz_pts[i]) = *ptr; ptr++;
        }
        
        assert(ptr == final_x + num_parameters);


	clock_t end = clock();
	
	printf("[RunSFM] run_sfm took %0.3fs\n",
	       (double) (end - start) / (double) CLOCKS_PER_SEC);

	/* Check for outliers */

	start = clock();

        double global_reprojection_error = 0;
        int global_num_observations = 0;

	std::vector<int> outliers;
        std::vector<int> outlier_views;
	std::vector<double> reproj_errors;
	    
	for (int i = 0; i < num_cameras; i++) {
            ImageData &data = m_image_data[added_order[i]];

	    double K[9] = { init_camera_params[i].f, 0.0, 0.0, 
			    0.0, init_camera_params[i].f, 0.0,
			    0.0, 0.0, 1.0 };
	    // double w[3] = { 0.0, 0.0, 0.0 };
	    double dt[3] = { init_camera_params[i].t[0],
			     init_camera_params[i].t[1],
			     init_camera_params[i].t[2] };		    

            /* Compute inverse distortion parameters */
            if (m_estimate_distortion) {
                double *k = init_camera_params[i].k;
                double k_dist[6] = { 0.0, 1.0, 0.0, k[0], 0.0, k[1] };
                double w_2 = 0.5 * data.GetWidth();
                double h_2 = 0.5 * data.GetHeight();
                double max_radius = 
                    sqrt(w_2 * w_2 + h_2 * h_2) / init_camera_params[i].f;

                InvertDistortion(6, 6, 0.0, max_radius, 
                                 k_dist, init_camera_params[i].k_inv);
            }

#ifdef __ENABLE_INTRINSICS__
            if (data.m_known_intrinsics) {
                double *k = init_camera_params[i].k_known;
                double k_dist[8] = 
                    { 0.0, 1.0, 0.0, k[0], 0.0, k[1], 0.0, k[4]};
                double w_2 = 0.5 * data.GetWidth();
                double h_2 = 0.5 * data.GetHeight();
                double max_radius = 
                    sqrt(w_2 * w_2 + h_2 * h_2) /
                        init_camera_params[i].K_known[0];

                InvertDistortion(8, 6, 0.0, max_radius, k_dist, 
                                 init_camera_params[i].k_inv);
            }
#endif /* __ENABLE_INTRINSICS__ */

	    int num_keys = GetNumKeys(added_order[i]);
	    
	    int num_pts_proj = 0;
	    for (int j = 0; j < num_keys; j++) {
		if (GetKey(added_order[i], j).m_extra >= 0) {
		    num_pts_proj++;
		}
	    }

	    double *dists = new double[num_pts_proj];
            std::vector<bool> in_back;
            in_back.resize(num_pts_proj);
	    int pt_count = 0;

            std::vector<Keypoint>::iterator iter;
	    // for (int j = 0; j < num_keys; j++) {
            for (iter = m_image_data[added_order[i]].m_keys.begin();
                 iter != m_image_data[added_order[i]].m_keys.end();
                 iter++) {

                const Keypoint &key = *iter;

		if (key.m_extra >= 0) {
		    double b[3], pr[2];
		    double dx, dy, dist;
		    int pt_idx = key.m_extra;
		    
		    b[0] = Vx(nz_pts[remap[pt_idx]]);
		    b[1] = Vy(nz_pts[remap[pt_idx]]);
		    b[2] = Vz(nz_pts[remap[pt_idx]]);

		    int in_front = sfm_project_rd(&(init_camera_params[i]), K, 
                                                  init_camera_params[i].k,
                                                  init_camera_params[i].R, 
                                                  dt, b, pr, 
                                                  m_estimate_distortion, 
                                                  true);

                    if (!in_front) {
                        printf("[RunSFM] Cheirality violation...\n");
                        in_back[pt_count] = true;
                    }

                    if (m_optimize_for_fisheye) {
                        /* Distort the points */
                        double x = pr[0], y = pr[1];
                        m_image_data[added_order[i]].
                            DistortPoint(x, y, pr[0], pr[1]);
                    }

		    dx = pr[0] - key.m_x;
		    dy = pr[1] - key.m_y;

		    dist = sqrt(dx * dx + dy * dy);
		    dist_total += dist;
		    num_dists++;

		    dists[pt_count] = dist;

		    pt_count++;
		}
	    }

            if (num_pts_proj == 0)
                continue;

	    /* Estimate the median of the distances */
	    double med = kth_element_copy(num_pts_proj, 
					  iround(0.8 /* 0.9 */ * num_pts_proj),
					  dists);
	    
	    median_copy(num_pts_proj, dists);

#define NUM_STDDEV 2.0 // 3.0 // 6.0
	    double thresh = 1.2 * NUM_STDDEV * med; /* k * stddev */
	    thresh = CLAMP(thresh, m_min_proj_error_threshold, 
                           m_max_proj_error_threshold);  

            /* Compute the average reprojection error for this
             * camera */

            double sum = 0.0;
            for (int j = 0; j < num_pts_proj; j++) {
                sum += dists[j];
            }

            double avg = sum / num_pts_proj;
            printf("[RunSFM] Mean error cam %d[%d] [%d pts]: %0.3e "
                   "[med: %0.3e, %0.3e]\n",
                   i, added_order[i], num_pts_proj, avg, 
                   kth_element_copy(num_pts_proj, 
                                    iround(0.5 * num_pts_proj), dists), 
                   thresh);

            global_reprojection_error += sum;
            global_num_observations += num_pts_proj;

	    pt_count = 0;
	    for (int j = 0; j < num_keys; j++) {
		int pt_idx = GetKey(added_order[i],j).m_extra;

		if (pt_idx < 0)
		    continue;

		/* Don't remove constrained points */
		if (m_use_point_constraints && 
		    Vx(m_point_constraints[pt_idx]) != 0.0) {
		    pt_count++;
		    continue;
		}
		
		if (dists[pt_count] > thresh || in_back[pt_count]) {
#define REMOVE_OUTLIER_VIEWS
#ifndef REMOVE_OUTLIER_VIEWS
		    /* Remove this point from consideration */
		    bool found = false;
		    for (int k = 0; k < (int) outliers.size(); k++) {
			if (outliers[k] == pt_idx) {
			    found = true;
			    break;
			}
		    }
#else
                    bool found = false;
#endif
		    
		    if (!found) {
			outliers.push_back(pt_idx);
                        outlier_views.push_back(i);
			reproj_errors.push_back(dists[pt_count]);
		    }
		}
		pt_count++;
	    }

#define OUTPUT_VERBOSE_STATS
#ifdef OUTPUT_VERBOSE_STATS
#define NUM_ERROR_BINS 10
            qsort(dists, num_pts_proj, sizeof(double), compare_doubles);

            double pr_min = dists[0];
            double pr_max = dists[num_pts_proj-1];
            double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

            /* Break histogram into 10 bins */
            int idx_count = 0;
            for (int i = 0; i < NUM_ERROR_BINS; i++) {
                double max = pr_min + (i+1) * pr_step;
                int start = idx_count;
                
                while (idx_count < num_pts_proj && dists[idx_count] <= max)
                    idx_count++;

                int bin_size = idx_count - start;
                printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n", 
                       max - pr_step, max, bin_size, 
                       bin_size / (double) num_pts_proj);
            }
#endif

	    delete [] dists;
	}

        printf("[RunSFM] Global mean reprojection error: %0.3e "
               "(%d observations)\n",
               global_reprojection_error / global_num_observations,
               global_num_observations);

	/* Remove outlying points */
        if ((!final_bundle || round > 0) && remove_outliers) {
            int num_dead = 0;
            for (int i = 0; i < (int) outliers.size(); i++) {
                int idx = outliers[i];
                int img_idx = outlier_views[i];

#ifndef REMOVE_OUTLIER_VIEWS
                printf("[RunSFM] Removing outlier %d "
                       "(reproj error: %0.3f)\n", idx, reproj_errors[i]);
#else
                printf("[RunSFM] Removing outlier view ( %d , %d ) "
                       "(reproj error: %0.3f)\n", idx, img_idx, 
                       reproj_errors[i]);
#endif
                
                int num_views = (int) pt_views[idx].size();
                int rem_idx = -1;

                for (int j = 0; j < num_views; j++) {
                    int v = pt_views[idx][j].first;
                    int k = pt_views[idx][j].second;
                    
                    // vmask[idx * num_cameras + v] = 0;
		
                    /* Sanity check */
                    if (GetKey(added_order[v], k).m_extra != idx)
                        printf("Error!  Entry for (%d,%d) "
                               "should be %d, but is %d\n",
                               added_order[v], k,
                               idx, GetKey(added_order[v], k).m_extra);

#ifndef REMOVE_OUTLIER_VIEWS
                    GetKey(added_order[v], k).m_extra = -2;
#else
                    if (v == img_idx) {
                        GetKey(added_order[v], k).m_extra = -2;
                        rem_idx = j;
                        break;
                    }
#endif
                }
                
#ifndef REMOVE_OUTLIER_VIEWS
                pt_views[idx].clear();
#else
                // assert(rem_idx != -1);
                if (rem_idx != -1) {
                    pt_views[idx].erase(pt_views[idx].begin() + rem_idx);
                    if (pt_views[idx].size() == 1) {
                        int v = pt_views[idx][0].first;
                        int k = pt_views[idx][0].second;
                        GetKey(added_order[v], k).m_extra = -2;
                        pt_views[idx].clear();

                        if (colors != NULL) {
                            Vx(colors[idx]) = 0x0;
                            Vy(colors[idx]) = 0x0;
                            Vz(colors[idx]) = 0xff;
                        }

                        num_dead++;
                    }
                }
#endif
            }
            
#ifndef REMOVE_OUTLIER_VIEWS
            num_outliers = outliers.size();
            total_outliers += num_outliers;
#else
            num_outliers = num_dead;
            total_outliers += num_outliers;
#endif

            end = clock();
            printf("[RunSFM] outlier removal took %0.3fs\n",
                   (double) (end - start) / (double) CLOCKS_PER_SEC);

            printf("[RunSFM] Removing %d outlier views\n", (int) outliers.size());
            printf("[RunSFM] %d points completely removed\n", num_dead);

            RemoveBadPointsAndCameras(num_pts, num_cameras, 
                                      added_order, init_camera_params, 
                                      init_pts, colors, pt_views);
        } else if (round == 0) {
            RemoveBadPointsAndCameras(num_pts, num_cameras, 
                                      added_order, init_camera_params, 
                                      init_pts, colors, pt_views);

            num_outliers = MIN_OUTLIERS + 1; /* so we don't stop */
        }
        
	// delete [] vmask;
	delete [] projections;
        delete [] init_x;
        delete [] cidx;
        delete [] pidx;
        delete [] num_vis;

        for (int i = 0; i < num_pts; i++) {
            if (remap[i] != -1) {
                init_pts[i] = nz_pts[remap[i]];
            }
        }

        if (!remove_outliers) break;

#if 1
        if (write_intermediate) {
            /* Write out the intermediate file */
            char buf[256];
            sprintf(buf, "%sBA.%d.out", m_bundle_output_base, round);
            DumpOutputFile(m_output_directory, buf, 
                           GetNumImages(), num_cameras, num_pts,
                           added_order, init_camera_params, init_pts, 
                           colors, pt_views);
        }
#endif

        gettimeofday(&stop_iter, NULL); 
        time_t time_iter = stop_iter.tv_sec - start_iter.tv_sec;
        printf("[RunSFM] Round %d took %d s\n", round, (int) time_iter);

        round++;

        fflush(stdout);
        fflush(stderr);

    } while (num_outliers > MIN_OUTLIERS);

    gettimeofday(&stop_all, NULL); 
    time_t time_all = stop_all.tv_sec - start_all.tv_sec;
    printf("[RunSFM] Structure from motion with outlier removal took %d s\n",
           (int) time_all);

    delete [] remap;
    delete [] nz_pts;

    return dist_total / num_dists;
}
