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

/* State: (q_i2b L_bi_i W P) */

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
  /* Normalize attitude quat q_i2b */
  quat_t q_i2b;
  quat_memcpy(&q_i2b, (quat_t *) y);
  quat_normalize(&q_i2b);

  double I_i_inv[9];
  double dcm_b2i[9];
  dcm_of_quat_b2a(dcm_b2i, &q_i2b);

  /* compute torques */
  const xyz_t torque_b = {0, 0, 0.5};
  xyz_t torque_i = {0.5, 0, 0};
  xyz_t body_torques_tmp;

  /* Convert body torques into inertial frame */
  rot_vec_by_quat_b2a(&body_torques_tmp, &q_i2b, &torque_b);
  fprintf(stderr, "T_bi_b (body subset): (%lf %lf %lf)\n",
          torque_b.x, torque_b.y, torque_b.z);
  fprintf(stderr, "T_bi_i (inertial subset): (%lf %lf %lf)\n",
          torque_i.x, torque_i.y, torque_i.z);
  /* dot(L_bi_i) = torque_i */
  xyz_sum((xyz_t *) &f[4], &body_torques_tmp, &torque_i);

  /* Rotate body-frame inertia tensor into inertial frame. */
  matrix_multiply(3, 3, 3, dcm_b2i, I_b_inv, I_i_inv);
  
  /* Compute current angular velocity */
  xyz_t w_bi_i, w;                      /* w = w_bi_b */
  matrix_multiply(3, 3, 1, I_i_inv, &y[4], (double *) &w_bi_i); /* y[4] = L_bi_i */
  rot_vec_by_quat_a2b(&w, &q_i2b, &w_bi_i);

  /* Compute differential attitude quaternion */
  f[0] = 0.5*(            - (w.x)*(q_i2b.q1) - (w.y)*(q_i2b.q2) - (w.z)*(q_i2b.q3));
  f[1] = 0.5*((w.x)*(q_i2b.q0)               + (w.z)*(q_i2b.q2) - (w.y)*(q_i2b.q3));
  f[2] = 0.5*((w.y)*(q_i2b.q0) - (w.z)*(q_i2b.q1)               + (w.x)*(q_i2b.q3));
  f[3] = 0.5*((w.z)*(q_i2b.q0) + (w.y)*(q_i2b.q1)   - (w.x)*(q_i2b.q2)            );

  /* Compute power and rotational impulse to track work and
   * momentum as a check on the solver */
  fprintf(stderr, "\t\tw_bi_i: (%lf %lf %lf)\n\t\tw_bi_b: (%lf %lf %lf)\n"
          "\t\tT_bi_i: (%lf %lf %lf)\n",
          w_bi_i.x, w_bi_i.y, w_bi_i.z,
          w.x, w.y, w.z,
          f[4], f[5], f[6]);
  f[7] = xyz_dot(&w, (xyz_t *) &f[4]);
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
  quat_t q_i2b_0 = {1, 0, 1, 0};
  xyz_t w_bi_b_0 = {0, 1, 0};

  quat_normalize(&q_i2b_0);

  double ke0;
  double I_tmp[3];
  matrix_multiply(1, 3, 3, (double *) &w_bi_b_0, I_b, I_tmp);
  matrix_multiply(1, 3, 1, I_tmp, (double *)&w_bi_b_0, &ke0);
  ke0 /= 2;
  
  xyz_t L_bi_i_0;
  xyz_t w_bi_i_0;
  double I_i[9], I_i_inv[9], dcm_b2i[9];
  dcm_of_quat_b2a(dcm_b2i, &q_i2b_0);
  /* xyz_mult_3x3_by_xyz(&w_bi_i_0, dcm_b2i, &w_bi_b_0); */
  rot_vec_by_quat_b2a(&w_bi_i_0, &q_i2b_0, &w_bi_b_0);
  /* Rotate body-frame inertia tensor into inertial frame. */
  matrix_multiply(3, 3, 3, dcm_b2i, I_b, I_i);
  xyz_mult_3x3_by_xyz(&L_bi_i_0, I_i, &w_bi_i_0);

  /* Integration parameters */
  double h = 1e-6, t = 0.0, t1 = 10;

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

  double y[SYS_SIZE] = {q_i2b_0.q0, q_i2b_0.q1, q_i2b_0.q2, q_i2b_0.q3,
                        L_bi_i_0.x, L_bi_i_0.y, L_bi_i_0.z,
                        ke0, xyz_norm(&L_bi_i_0)};

  /* Step using GSL integrator */
  while (t < t1) {
    if ((status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, y))
        != GSL_SUCCESS)
      break;

    /* Spew a lot of stuff */
    printf("%2.5lf\t", t);
    for (i = 0; i < 4; i++) printf("%.5e ", y[i]);
    printf("\t");

    xyz_t w_bi_i, w_bi_b;
    dcm_of_quat_b2a(dcm_b2i, (quat_t *) y);
      
    matrix_multiply(3, 3, 3, dcm_b2i, I_b_inv, I_i_inv);
    matrix_multiply(3, 3, 1, I_i_inv, &y[4], (double *) &w_bi_i);
    rot_vec_by_quat_a2b(&w_bi_b, (quat_t *) y, &w_bi_i);

    printf("%.5e %.5e %.5e\t", w_bi_b.x, w_bi_b.y, w_bi_b.z);

    euler_t eulers;
    euler321_of_quat(&eulers, (quat_t *) y);
    printf("%.5e %.5e %.5e\t", eulers.roll, eulers.pitch, eulers.yaw);

    /* Check quaternion norm */
    double qnorm = y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
    if (fabs(qnorm - 1) > NUM_TOL)
      fprintf(stderr, "Quaternion norm differs from 1: %lf\n", qnorm);

    /* Check work and energy */
    double ke;
    matrix_multiply(1, 3, 3, (double *) &w_bi_b, I_b, I_tmp);
    matrix_multiply(1, 3, 1, I_tmp, (double *) &w_bi_b, &ke);
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
