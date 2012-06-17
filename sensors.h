/***********************************************************/
/* dynamics.h                                              */
/* Copyright (C) 2012 Matthew Peddie <peddie@alum.mit.edu> */
/* Attitude dynamics sensors                               */
/***********************************************************/

#include <xyz.h>

#include "dynamics.h"

#ifndef __ATTSIM_SENSORS_H__
#define __ATTSIM_SENSORS_H__

typedef struct sensors_params {
  xyz_t mag_noise_b;
} sensors_params;

typedef struct measurements {
  xyz_t mag_b;
  xyz_t pqr;
} measurements;

int sensor_outputs(const full_state *state, const dynamics_params *dp,
                   measurements *meas);
int sensors_init(void);

#endif  /* __ATTSIM_SENSORS_H__ */
