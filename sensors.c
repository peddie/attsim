#include <quat.h>
#include <xyz.h>
#include <spatial_rotations.h>

#include "dynamics.h"

static const sensors_params p = {
  .mag_noise_b = {
    .x = 0.0,
    .y = 0.0,
    .z = 0.0
  }
};

static int
sense_mat_b(const full_state *state, const dynamics_params *dp, xyz_t *mag_b)
{
  xyz_t mag_lvlh;
  dp->compute_mag(state->t, &mag_lvlh);
  rot_vec_by_quat_a2b(mag_b, &state->q_i2b, &mag_lvlh);

  return 0;
}

int
sensor_outputs(const full_state *state,
               const dynamics_params *dp,
               measurements *meas)
{
  return sense_mat_b(state, dp, &meas->mag_b);
}

int
sensors_init(void)
{
  return 0;
}
