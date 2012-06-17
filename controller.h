#ifndef __ATTSIM_CONTROLLER_H__
#define __ATTSIM_CONTROLLER_H__

#include <xyz.h>
#include <quat.h>

#include "sensors.h"

typedef struct estimator_params {
  double kp_pqr_bias;
  double kp_mag_bias;
} estimator_params;

typedef struct estimator_state {
  quat_t q_lvlh2b;
  xyz_t w_blvlh_b;
  xyz_t pqr_bias;
  xyz_t mag_bias;
} estimator_state;

typedef struct bdot_params {
  const double diff_tau;
  const double kp;
} bdot_params;

typedef struct bdot_state {
  xyz_t mag_b_1, mag_b_dot;
  xyz_t mag_b_dot_f;
} bdot_state;

typedef struct controller_params {
  const bdot_params bdot;
} controller_params;

typedef struct controller_state {
  bdot_state bdot;
  estimator_state est;
} controller_state;

int run_controller(double ts, const measurements *meas, actuators *act);
int controller_init(void);

#endif  // __ATTSIM_CONTROLLER_H__

