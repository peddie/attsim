#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>

#include <xyz.h>
#include <spatial_rotations.h>

#define MAX_ITER 22
#define NUM_TOL 1e-9
#define SYS_SIZE 9

const double I_b[9] = {1, 0, 0,
                       0, 1, 0,
                       0, 0, 1};
double I_b_inv[9];

/* State: (q0 q1 q2 q3 Lx Ly Lz W) */

static inline void matrix_multiply(size_t n, size_t m, size_t p,
                                   const double *a, const double *b,
                                   double *c) {
  unsigned int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < p; j++) {
      c[p*i + j] = 0;
      for (k = 0; k < m; k++)
        c[p*i + j] += a[m*i+k] * b[p*k + j];
    }
}


static inline int inv3(const double *a, double *b) {
  double det = ((a[3*1 + 0]*-(a[3*0 + 1]*a[3*2 + 2]-a[3*0 + 2]*a[3*2 + 1])
                 +a[3*1 + 1]*(a[3*0 + 0]*a[3*2 + 2]-a[3*0 + 2]*a[3*2 + 0]))
                +a[3*1 + 2]*-(a[3*0 + 0]*a[3*2 + 1]-a[3*0 + 1]*a[3*2 + 0]));

  if (det < NUM_TOL)
    return -1;

  b[3*0 + 0] = (a[3*1 + 1]*a[3*2 + 2]-a[3*1 + 2]*a[3*2 + 1])/det;
  b[3*1 + 0] = -(a[3*1 + 0]*a[3*2 + 2]-a[3*1 + 2]*a[3*2 + 0])/det;
  b[3*2 + 0] = (a[3*1 + 0]*a[3*2 + 1]-a[3*1 + 1]*a[3*2 + 0])/det;

  b[3*0 + 1] = -(a[3*0 + 1]*a[3*2 + 2]-a[3*0 + 2]*a[3*2 + 1])/det;
  b[3*1 + 1] = (a[3*0 + 0]*a[3*2 + 2]-a[3*0 + 2]*a[3*2 + 0])/det;
  b[3*2 + 1] = -(a[3*0 + 0]*a[3*2 + 1]-a[3*0 + 1]*a[3*2 + 0])/det;

  b[3*0 + 2] = (a[3*0 + 1]*a[3*1 + 2]-a[3*0 + 2]*a[3*1 + 1])/det;
  b[3*1 + 2] = -(a[3*0 + 0]*a[3*1 + 2]-a[3*0 + 2]*a[3*1 + 0])/det;
  b[3*2 + 2] = (a[3*0 + 0]*a[3*1 + 1]-a[3*0 + 1]*a[3*1 + 0])/det;

  return 0;
}

static int
attitude(double t __attribute__((unused)),
         const double y[], double f[], void *params __attribute__((unused)))
{
  /* Normalize attitude quat q_b2i */
  quat_t q;
  quat_memcpy(&q, (quat_t *) y);
  quat_normalize(&q);

  double I_i_inv[9];
  double dcm_b2i[9];
  dcm_of_quat_a2b(dcm_b2i, &q);

  /* compute torques */
  const xyz_t torque_b = {0.5, 0, 0};
  xyz_t torque_i = {0, 0, 0};
  xyz_t body_torques_tmp;

  /* Convert body torques into inertial frame */
  rot_vec_by_quat_a2b(&body_torques_tmp, &q, &torque_b);
/*   matrix_multiply(3, 3, 1, dcm_b2i, (double *) &torque_b, */
/*                   (double *) &body_torques_tmp); */
  xyz_sum((xyz_t *) &f[4], &body_torques_tmp, &torque_i);

  /* Rotate body-frame inertia tensor into inertial frame. */
  matrix_multiply(3, 3, 3, dcm_b2i, I_b_inv, I_i_inv);
  
  /* Compute current angular velocity */
  xyz_t w_i, w;
  matrix_multiply(3, 3, 1, I_i_inv, &y[4], (double *) &w_i);
  rot_vec_by_quat_b2a(&w, &q, &w_i);

  /* Compute differential attitude quaternion */
  f[0] = 0.5*(              - (w.x)*(q.q1) - (w.y)*(q.q2) - (w.z)*(q.q3)  );
  f[1] = 0.5*(  (w.x)*(q.q0)               + (w.z)*(q.q2) - (w.y)*(q.q3)  );
  f[2] = 0.5*(  (w.y)*(q.q0) - (w.z)*(q.q1)               + (w.x)*(q.q3)  );
  f[3] = 0.5*(  (w.z)*(q.q0) + (w.y)*(q.q1)   - (w.x)*(q.q2)              );

  /* Compute power and rotational impulse to track work and
   * momentum as a check on the solver */
  f[7] = xyz_dot(&w_i, (xyz_t *) &f[4]);
  f[8] = xyz_norm((xyz_t *) &f[4]);
  
  return GSL_SUCCESS;
}

int
main (void)
{
  int i, status;

  /* Pre-compute inverse of body-frame inertia tensor */
  inv3(I_b, I_b_inv);

  /* Initial conditions */
  quat_t q0 = {1, 0, 0, 0};
  xyz_t w0 = {0, 1, 0};

  quat_normalize(&q0);

  double ke0;
  double I_tmp[3];
  matrix_multiply(1, 3, 3, (double *) &w0, I_b, I_tmp);
  matrix_multiply(1, 3, 1, I_tmp, (double *)&w0, &ke0);
  ke0 /= 2;
  
  xyz_t L0;
  xyz_t w_i;
  double I_i[9], I_i_inv[9], dcm_b2i[9];
  dcm_of_quat_b2a(dcm_b2i, &q0);
  xyz_mult_3x3_by_xyz(&w_i, dcm_b2i, &w0);
  /* Rotate body-frame inertia tensor into inertial frame. */
  matrix_multiply(3, 3, 3, dcm_b2i, I_b, I_i);
  xyz_mult_3x3_by_xyz(&L0, I_i, &w_i);

  /* Integration parameters */
  double h = 1e-3, t = 0.0, t1 = 10;

  /* GSL setup */
  const gsl_odeiv2_step_type * T
      = gsl_odeiv2_step_rk8pd;
     
  gsl_odeiv2_step * s
      = gsl_odeiv2_step_alloc(T, SYS_SIZE);
  gsl_odeiv2_control * c
      = gsl_odeiv2_control_standard_new(1e-10, 1e-10, 1, 1);
  gsl_odeiv2_evolve * e
      = gsl_odeiv2_evolve_alloc(SYS_SIZE);

  if (!s || !c || !e) fprintf(stderr, "Eit!\n");

  gsl_odeiv2_system sys = {attitude, NULL, SYS_SIZE, NULL};

  double y[SYS_SIZE] = {q0.q0, q0.q1, q0.q2, q0.q3,
                        L0.x, L0.y, L0.z,
                        ke0, xyz_norm(&L0)};

  /* Step using GSL integrator */
  while (t < t1) {
    if ((status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, y))
        != GSL_SUCCESS)
      break;

    /* Spew a lot of stuff */
    printf("%2.5lf\t", t);
    for (i = 0; i < 4; i++) printf("%.5e ", y[i]);
    printf("\t");

    xyz_t w;
    dcm_of_quat_a2b(dcm_b2i, (quat_t *) y);
      
    matrix_multiply(3, 3, 3, dcm_b2i, I_b_inv, I_i_inv);
    matrix_multiply(3, 3, 1, I_i_inv, &y[4], (double *) &w);

    printf("%.5e %.5e %.5e\t", w.x, w.y, w.z);

    euler_t eulers;
    euler321_of_quat(&eulers, (quat_t *) y);
    printf("%.5e %.5e %.5e\t", eulers.roll, eulers.pitch, eulers.yaw);

    /* Check quaternion norm */
    double qnorm = y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
    if (fabs(qnorm - 1) > NUM_TOL)
      fprintf(stderr, "Quaternion norm differs from 1: %lf\n", qnorm);

    /* Check work and energy */
    double ke;
    matrix_multiply(1, 3, 3, (double *) &w, I_b, I_tmp);
    matrix_multiply(1, 3, 1, I_tmp, (double *)&w, &ke);
    ke /= 2;

    printf("%.5e %.5e %.5e\t", y[7], ke, y[7] - ke);

    /* Check momentum */
    printf("%.5e %.5e\t", xyz_norm((xyz_t *) &y[4]), y[8]);

    printf("%.5e %.5e %.5e\t", y[4], y[5], y[6]);
    printf("\n");
  }

  /* GSL cleanup */
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);
  return 0;
}
