/***********************************************************/
/* dynamics.h                                              */
/* Copyright (C) 2012 Matthew Peddie <peddie@alum.mit.edu> */
/* Attitude dynamics                                       */
/***********************************************************/

#include <stddef.h>

#include <xyz.h>
#include <quat.h>

#include "controller.h"

#ifndef __ATTSIM_DYNAMICS_H__
#define __ATTSIM_DYNAMICS_H__

#define NUM_TOL 1e-12
#define SYS_SIZE 9

/* Callback type to compute magnetic field */
typedef int(*mag_lvlh_f)(double, xyz_t *);

typedef struct actuators {
  xyz_t mag_coil;
} actuators;

typedef struct dynamics_params {
  double I_b[9];
  double I_b_inv[9];
  actuators act;
  mag_lvlh_f compute_mag;
} dynamics_params;

typedef struct full_state {
  double t;
  quat_t q_i2b;
  xyz_t w_bi_b;
  double P, T;
} full_state;

/* Math helpers */
void matrix_multiply(size_t n, size_t m, size_t p, const double *a,
                       const double *b, double *c);
int inv3(const double *a, double *b);

/* Dynamics helpers */
double compute_kinetic_energy(const full_state *s, const dynamics_params *dp);
void compute_momentum(const full_state *s, const dynamics_params *dp,
                      xyz_t *L_bi_b);


/* Attitude dynamics */
int attitude(double t, const double y[], double f[], void *params);
int dynamics_init(double y0[SYS_SIZE], dynamics_params *dp);

#endif  /* __ATTSIM_DYNAMICS_H__ */
