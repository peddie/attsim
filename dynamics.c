/***********************************************************/
/* dynamics.c                                              */
/* Copyright (C) 2012 Matthew Peddie <peddie@alum.mit.edu> */
/* Attitude dynamics                                       */
/***********************************************************/

#include <math.h>
#include <string.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>

#include <xyz.h>
#include <spatial_rotations.h>

#include "dynamics.h"


/* State: (q_i2b L_bi_i W P) */

/* math functions */
inline void matrix_multiply(size_t n, size_t m, size_t p,
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


inline int inv3(const double *a, double *b) {
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

double
compute_kinetic_energy(const full_state *s, const dynamics_params *dp)
{
  double T;
  xyz_t wJ;
  matrix_multiply(1, 3, 3, (double *) &s->w_bi_b, dp->I_b, (double *) &wJ);
  matrix_multiply(1, 3, 1, (double *) &wJ, (double *) &s->w_bi_b, &T);
  T /= 2;

  return T;
}

void
compute_momentum(const full_state *s, const dynamics_params *dp,
                 xyz_t *L_bi_b)
{
  xyz_mult_3x3_by_xyz(L_bi_b, dp->I_b, &s->w_bi_b);
}

static int
compute_torques(const full_state *state __attribute__((unused)),
                xyz_t *torque_b, xyz_t *torque_i)
{
  torque_b->x = 0.5;
  torque_b->y = 0.0;
  torque_b->z = 0.0;

  torque_i->x = 0.0;
  torque_i->y = 0.0;
  torque_i->z = 0.5;

  return 0;                             /* computation went OK */
}

int
attitude(double t __attribute__((unused)), const double y[],
         double f[], void *params)
{
  int fail = 0;
  /* Unmarshal states */
  xyz_t w;                              /* w_bi_b */
  xyz_memcpy(&w, (xyz_t *) &y[4]);
  
  quat_t q;                             /* q_i2b */
  quat_memcpy(&q, (quat_t *) y);

  dynamics_params *p = (dynamics_params *) params;

  /* Normalize attitude quat */
  quat_normalize(&q);

  /* get torques, computed in different frames */
  xyz_t torque_b, torque_i;
  full_state s;
  xyz_memcpy(&s.w_bi_b, &w);
  quat_memcpy(&s.q_i2b, &q);
  if ((fail = compute_torques(&s, &torque_b, &torque_i)) != 0)
    return 22;
  
  /* Convert inertial-frame torques into the body frame */
  xyz_t inertial_torques_tmp, m_total;
  rot_vec_by_quat_a2b(&inertial_torques_tmp, &q, &torque_i);
  xyz_sum(&m_total, &inertial_torques_tmp, &torque_b);

  /* \dot{\omega} = J^{-1} (M - \omega \times (J \omega)) */
  xyz_t Jw, wJw, MwJw;
  xyz_mult_3x3_by_xyz(&Jw, p->I_b, &w);
  xyz_cross(&wJw, &w, &Jw);
  xyz_subtract(&MwJw, &m_total, &wJw);
  xyz_mult_3x3_by_xyz((xyz_t *) &f[4], p->I_b_inv, &MwJw);

  /* Compute differential attitude quaternion */
  f[0] = 0.5*(            - (w.x)*(q.q1) - (w.y)*(q.q2) - (w.z)*(q.q3));
  f[1] = 0.5*((w.x)*(q.q0)               + (w.z)*(q.q2) - (w.y)*(q.q3));
  f[2] = 0.5*((w.y)*(q.q0) - (w.z)*(q.q1)               + (w.x)*(q.q3));
  f[3] = 0.5*((w.z)*(q.q0) + (w.y)*(q.q1)   - (w.x)*(q.q2)            );

  /* Compute a correction quaternion from the GSL state and the
   * quaternion we re-normalized; use it to make corrections to the
   * GSL state updates. */
  f[0] += q.q0 - y[0];
  f[1] += q.q1 - y[1];
  f[2] += q.q2 - y[2];
  f[3] += q.q3 - y[3];

  /* Compute power and rotational impulse to track work and
   * momentum as a check on the solver */
  f[7] = xyz_dot(&w, &m_total);

  xyz_t L_bi_b;
  xyz_mult_3x3_by_xyz(&L_bi_b, p->I_b, (xyz_t *) &f[4]);
  f[8] = xyz_norm(&L_bi_b);

  return GSL_SUCCESS;
}

int
dynamics_init(double y0[], dynamics_params *dp)
{
    /* Pre-compute inverse of body-frame inertia tensor */
  double J[9] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1};
  memcpy(dp->I_b, J, sizeof(J));
  inv3(dp->I_b, dp->I_b_inv);
  
  /* Initial conditions */
  quat_t q_i2b_0 = {1, 0, 0, 0};
  xyz_t w_bi_b_0 = {0, 1, 0};
  

  quat_normalize(&q_i2b_0);

  /* Initial conditions for conserved parameters */
  double ke0;
  double I_tmp[3];
  matrix_multiply(1, 3, 3, (double *) &w_bi_b_0, dp->I_b, I_tmp);
  matrix_multiply(1, 3, 1, I_tmp, (double *)&w_bi_b_0, &ke0);
  ke0 /= 2;

  xyz_t L_bi_b_0;
  xyz_mult_3x3_by_xyz(&L_bi_b_0, dp->I_b, &w_bi_b_0);
  
  quat_memcpy((quat_t *) y0, &q_i2b_0);
  xyz_memcpy((xyz_t *) &y0[4], &w_bi_b_0);
  y0[7] = ke0;
  y0[8] = xyz_norm(&L_bi_b_0);

  return 0;
}
