// Standard C headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Mathlib headers
#include <xyz.h>
#include <spatial_rotations.h>

// Our headers
#include "sdl_draw.h"
#include "aero_model.h"
#include "aero_table.h"

#define UU __attribute__((unused))

#define MEGARES 0.3
#define DRAW

struct indirect_element {
	int index;
	double val;
};

int compare_indirect_elements(const void *a, const void *b);

int project_o_tron(aero_table_output_t * out, aero_model_t * model,
		   quat_t ref_to_body, int w, int h, int *pix_buf);

int att_of_params(quat_t *att, double a, double b, int c);

void params_of_att(double *a, double *b, int *c, const quat_t *att);

int compare_indirect_elements(const void *a, const void *b)
{
	double diff =
	    ((struct indirect_element *)a)->val -
	    ((struct indirect_element *)b)->val;
	if (diff > 0)
		return 1;
	if (diff < 0)
		return -1;
	return 0;
}

int project_o_tron(aero_table_output_t * out,
		   aero_model_t * model,
		   quat_t ref_to_body, int w, int h, int *pix_buf_panel)
{

	int n = model->n_panels;

	// *** Input parameters
	xyz_t UU e_v_ref = { 1, 0, 0 };	// velocity unit vector in ref frame
	// todo: check that is normalized
	double res = 0.001 / MEGARES;	// Resolution
	double sig_n = 0.8;	// normal momentum exchange coeff
	double sig_t = 0.8;	// tangential momentum exchange coeff

	// *** Let's choose a new "freestream" frame so that the velocity vector is just +X
	quat_t body_to_fs;
	xyz_t e_v_fs = { 1, 0, 0 };
	// Don't need to rotate at all, the frames are already aligned
	// Now form the quaternion rotating the body frame to the freestream frame
	quat_inv(&body_to_fs, &ref_to_body);

	// *** Rotate "everything" into the new frame
	//

	//  Rotate the center of mass
	xyz_t cg_fs;
	rot_vec_by_quat_a2b(&cg_fs, &body_to_fs, &model->cg);

	//  Rotate the panels' coordinates 
	panel_t panels_fs[model->n_panels];
	for (int i = 0; i < n; i++) {
		rot_vec_by_quat_a2b(&panels_fs[i].origin, &body_to_fs,
				    &model->panels[i].origin);
		rot_vec_by_quat_a2b(&panels_fs[i].side_a, &body_to_fs,
				    &model->panels[i].side_a);
		rot_vec_by_quat_a2b(&panels_fs[i].side_b, &body_to_fs,
				    &model->panels[i].side_b);
	}

	//  Find the centroids
	xyz_t centroid_fs[n];
	for (int i = 0; i < n; i++) {
		centroid_fs[i] = panels_fs[i].origin;
		//  This handy function does (a <= a*A + b*B + c*C)
		xyz_sum_scale(&centroid_fs[i], 1, &panels_fs[i].side_a, 0.5,
			      &panels_fs[i].side_b, 0.5);
	}

	// Sort by the X coordinate of the centroids
	struct indirect_element sort_list[n];
	for (int i = 0; i < n; i++) {
		sort_list[i].index = i;
		sort_list[i].val = centroid_fs[i].x;
	}
	qsort(sort_list, n, sizeof(sort_list[0]), compare_indirect_elements);

	double *pix_buf_x = malloc(w * h * sizeof(double));
	// Holds the fs frame X coordinate
	// of the panel element at each pixel
	if (!pix_buf_x) {
		fprintf(stderr, "Unable to allocate %zu bytes for pix_buf_x\n",
			w * h * sizeof(double));
		return 1;
	}
	// Iterate over the panels, rearmost centroid first
	for (int i = 0; i < n; i++) {	// if(sort_list[i].index==0){
		panel_t *p = &panels_fs[sort_list[i].index];
		//printf("Panel %d:::\n",sort_list[i].index);

		double res_a = res / xyz_norm(&p->side_a) / sqrt(2);
		double res_b = res / xyz_norm(&p->side_b) / sqrt(2);

		// Iterate over side a
		for (double a = 0; a <= 1; a += res_a) {
			//  Iterate over side b
			for (double b = 0; b <= 1; b += res_b) {
				xyz_t panel_cell = p->origin;
				xyz_sum_scale(&panel_cell, 1, &p->side_a, a,
					      &p->side_b, b);
				int pix_y =
				    (int)round(panel_cell.y / res + w / 2);
				int pix_z =
				    (int)round(panel_cell.z / res + h / 2);
				pix_buf_panel[pix_z + pix_y * w] =
				    sort_list[i].index + 1;
				pix_buf_x[pix_z + pix_y * w] = panel_cell.x;
				//   printf("\t%+7.4f\t%+7.4f\t%d\t%d\n",panel_cell.y,panel_cell.z,pix_y,pix_z);
			}
		}
	}

	// Find the unit normal vector in fs frame, and e_v.e_n for all panels
	xyz_t e_n_fs[n];
	double e_v_dot_e_n[n];
	for (int p = 0; p < n; p++) {
		xyz_cross(&e_n_fs[p], &panels_fs[p].side_a,
			  &panels_fs[p].side_b);
		xyz_normalize_self(&e_n_fs[p]);
		e_v_dot_e_n[p] = xyz_dot(&e_n_fs[p], &e_v_fs);
//    if (e_v_dot_e_n[p] < 0)
//      printf("e_v dot e_n negative for panel %d\n",p);
	}

//  double A = 0;
	int n_pix[n];
	memset(n_pix, 0, sizeof(n_pix));
	xyz_t panel_sum_pos[n];
	memset(panel_sum_pos, 0, sizeof(panel_sum_pos));

	//printf("sizeof(n_pix)=%zu\n",sizeof(n_pix));

  xyz_t F_fs = {0, 0, 0};
  xyz_t T_fs = {0, 0, 0};

	// Now iterate over the pixels
	for (int pix_y = 0; pix_y < w; pix_y++)
		for (int pix_z = 0; pix_z < h; pix_z++) {
			int p = pix_buf_panel[pix_z + pix_y * w];
			if (p) {
				// There is a panel at this pixel.
        p--;  // Actual panel index number is one less

				if (e_v_dot_e_n[p] < 0) {
					// Skip panels facing the wrong way
					continue;
				}

				n_pix[p]++;

				// Find the coordinates in the fs frame
				xyz_t fs_p;
				fs_p.x = pix_buf_x[pix_z * w + pix_y];
				fs_p.y = (pix_y - w / 2) * res;
				fs_p.z = (pix_z - h / 2) * res;

				// Where is it wrt the center of mass?
				xyz_t r_from_cg_fs;
				xyz_diff(&r_from_cg_fs, &fs_p, &cg_fs);

				// Let's compute the forces on this element
				// ref: NASA SP-8058 eq 2-2 (page 5).
				// Dynamic pressure (1/2 rho v^2) factored out, hence leading 2 *
				xyz_t dF = { 0, 0, 0 };
				xyz_sum_scale(&dF, 0,
					      &e_n_fs[p],
					      2 * -1 * res * res * (2 - sig_n -
							       sig_t) *
					      xyz_dot(&e_v_fs, &e_n_fs[p]),
					      &e_v_fs, 2 * -1 * res * res * sig_t);

				// Torque contribution from this element
				xyz_t dT;
				xyz_cross(&dT, &r_from_cg_fs, &dF);

				xyz_sum(&F_fs, &F_fs, &dF);
				xyz_sum(&T_fs, &T_fs, &dT);

			}
		}


  rot_vec_by_quat_b2a(&out->force,  &body_to_fs, &F_fs);
  rot_vec_by_quat_b2a(&out->torque, &body_to_fs, &T_fs);
  
  printf("F = %f %f %f (norm=%f)\n",out->force.x, out->force.y, out->force.z, xyz_norm(&out->force));
  printf("T = %f %f %f (norm=%f)\n",out->torque.x, out->torque.y, out->torque.z, xyz_norm(&out->torque));
//	 printf("%d pixels dropped due to e_v dot e_n negative\n",n_pix);

	free(pix_buf_x);
	return 0;
}


void params_of_att(double *a, double *b, int *c, const quat_t *att) {
  const xyz_t e_x={1,0,0};
  xyz_t v;
  rot_vec_by_quat_a2b(&v, att, &e_x);
  *a = v.y;
  *b = v.z;
  *c = (v.x > 0);
}


int att_of_params(quat_t *att, double a, double b, int c) {
  
  xyz_t v, cross;
  quat_t q;
  int flag=0;

  v.y = a;
  v.z = b;
  v.x = sqrt(1 - a*a - b*b);
  
  if (isnan(v.x)) {
    fprintf(stderr,"att_of_params: %.4f^2 + %.4f^2 > 1.  Normalizing and ploughing on regardless.\n",a,b);
//    return 1;       //invalid parameters
    v.x = 0;
    xyz_normalize_self(&v);
  }
  
  if (!c) {
    v.x = -v.x;
    if (v.x < -.9999) {  // Fix singularity
      printf(" singularity\n");
      v.y = 1;
      flag=1;
    }
  }



  cross.x = 0;
  cross.y = -v.z;
  cross.z = v.y;

  q.q1 = cross.x;
  q.q2 = cross.y;
  q.q3 = cross.z;

  q.q0 = 1 + v.x;
  if (quat_norm(&q) < 0.1)
  printf(" q_of_params: norm was %.4f for %.3f, %.3f, %d (v.x = %.3f)\n",quat_norm(&q),a,b,c,v.x);

  if (flag)
    printf(" %.3f, %.3f, %d -> %.3f, %.3f, %.3f, %.3f\n",a,b,c,q.q0,q.q1,q.q2,q.q3);

  quat_normalize(&q);
 
  if (flag)
    printf(" and after normalizing, that's %.3f, %.3f, %.3f, %.3f\n",q.q0,q.q1,q.q2,q.q3);
  
  quat_inv(att,&q);


  return 0;

}




int main(int argc, char *argv[])
{
	char *filename_base;
	char filename_table[100];

  aero_model_t model;
	aero_env_t env;
	aero_table_output_t result;

	load_aero_model(&model, "");
	load_aero_env(&env, "");

  int w = 2 * 600 * MEGARES;
	int h = w;
	int *pix_buf = malloc(w * h * sizeof(int));
	if (!pix_buf) {
		fprintf(stderr, "Unable to allocate %zu bytes for pix_buf\n",
			w * h * sizeof(int));
		return 1;
	};

	if (argc > 2) {
		if (strlen(argv[2]) + 10 > sizeof(filename_table)) {
			printf("Specify a shorter filename, please (max %zu)\n",
			       sizeof(filename_table) - 10);
      return 1;
    }
		else
			filename_base = argv[2];
	}
  else {
    printf("syntax: %s {w|r} basename\n",argv[0]);
    return 1;
  }


	snprintf(filename_table, sizeof(filename_table), "%s.table",
		 filename_base);

  if (argv[1][0]=='r') {
    aero_table_t table;
    load_aero_table(&table,filename_table);
    
    aero_table_output_t result;
   
		sdl_open(w, h);

    for (int c=0; c<=1; c++)
    for (double a=-1; a<=1; a+=0.01)
      for (double b=-1; b<=1; b+=0.01) {
        if (a*a + b*b > 1.0001)
          continue;
        if (lookup_aero_table(&result,&table,a,b,c)) {
          printf("lookup fail for %.3f, %.3f, %d\n",a,b,c);
          continue;
        }
        double torque = xyz_norm(&result.torque);
        if (torque < 0.00015) {
          printf("Low torque position found: %.3f, %.3f, %d  (%.4f)\n",a,b,c, torque);
          quat_t att;

          att_of_params(&att,a,b,c);

			    memset(pix_buf, 0, w * h * sizeof(*pix_buf));
    			project_o_tron(&result, &model, att, w, h, pix_buf);

			    sdl_plot(w, h, pix_buf);

			    if (sdl_waitkey()) return 1;

        }

      }

    
    close_aero_table(&table);
    return 0;
  }

 if (argv[1][0]=='a') {
    aero_table_t table;
    load_aero_table(&table,filename_table);
    
    aero_table_output_t result;
   
		sdl_open(w, h);

    double a,b;
    int c;

    sscanf(argv[3],"%lf",&a);
    sscanf(argv[4],"%lf",&b);
    sscanf(argv[5],"%d",&c);


        if (lookup_aero_table(&result,&table,a,b,c)) {
          printf("lookup fail for %.3f, %.3f, %d\n",a,b,c);
          return 1;
        }
        double torque = xyz_norm(&result.torque);
          printf("Params: %.3f, %.3f, %d  (%.4f)\n",a,b,c, torque);
          quat_t att;

          att_of_params(&att,a,b,c);

			    memset(pix_buf, 0, w * h * sizeof(*pix_buf));
    			project_o_tron(&result, &model, att, w, h, pix_buf);

			    sdl_plot(w, h, pix_buf);

			    if (sdl_waitkey()) return 1;

        
    close_aero_table(&table);
    return 0;
  }



	FILE UU *table_file = fopen(filename_table, "wb");


  const double param_step = 0.02;
	

	int UU o = 0;


	double max_torque = 0;

  int UU col=1;

  int coverage[600][600];
  memset(coverage, 0, sizeof(coverage));

	memset(pix_buf, 0, w * h * sizeof(*pix_buf));


  int param_dim = round(2/param_step + 1);

  
  for (int c=0; c<2; c++)
  for (int ia = 0; ia < param_dim; ia ++) 
    {
  for (int ib = 0; ib < param_dim; ib ++) 
		{

    double a,b;
    a = 2.0 * ia / (param_dim-1) - 1;
    b = 2.0 * ib / (param_dim-1) - 1;

    if (a*a + b*b > 1.0001) {
      printf("skipping %f, %f\n",a,b);
      continue;
    }
    else   printf("@ %f, %f \n",a,b);


    quat_t att;

    att_of_params(&att,a,b,c);

/*
    if (c) 
      pix_buf[(int)(100*a)+120 + w * ((int)(100*b)+120)] = col;
    else
      pix_buf[(int)(100*a)+420 + w * ((int)(100*b)+120)] = col;
  
    col%=7;
    col++;
*/

			// printf("norm = %f\n",quat_norm(&att));

			memset(pix_buf, 0, w * h * sizeof(*pix_buf));

  

			project_o_tron(&result, &model, att, w, h, pix_buf);

			write_aero_table_entry(table_file, a, b,c, &result);

#ifdef DRAW
      /*
			// Plot some color squares at the top, as a crude key / legend
			for (int i = 0; i < model.n_panels; i++)
				for (int y = 0; y < 22; y++)
					for (int x = 0; x < 22; x++)
						pix_buf[x + 22 * 2 * i +
							y * w] = i + 1;
      */

			if (!o) {
				sdl_open(w, h);
				o = 1;
			}
			sdl_plot(w, h, pix_buf);
			if (sdl_poll_key()) goto end;
    //  if (sdl_waitkey()) goto end;
#endif

			double torque =
			    .5 * env.rho * env.v * env.v *
			    xyz_norm(&result.torque);
			if (torque > max_torque) {
				max_torque = torque;
				printf
				    ("New max torque, %.2f microNewton meters\n",
				     max_torque / 1E-6);
    //   if (torque > 4e-6 &&  sdl_waitkey()) goto end;
			}
		}
	}
#ifdef DRAW
 end:
  sdl_waitkey();
	sdl_close();
#endif
  fclose(table_file);
	free(pix_buf);
	return 0;
}
