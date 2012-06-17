#include <xyz.h>
#include <filters.h>
#include <spatial_rotations.h>

#include "sensors.h"
#include "controller.h"

#define BDOT_GAIN 0.2

static controller_state st;

static const controller_params p = {
  .bdot = {
    .diff_tau = 5,
    .kp = 1
  }
};

static int
bdot_init(void) {
  xyz_set_all(&st.bdot.mag_b_1, 0);
  xyz_set_all(&st.bdot.mag_b_dot, 0);
  xyz_set_all(&st.bdot.mag_b_dot_f, 0);
  return 0;
}

static int
bdot(double ts, const measurements *meas, const bdot_params *bp,
     bdot_state *bst, actuators *act)
{
  /* differentiate magnetic field */
  xyz_diff(&bst->mag_b_dot, &meas->mag_b, &bst->mag_b_1);
  xyz_memcpy(&bst->mag_b_1, &meas->mag_b);
  
  /* 1st-order LPF on bdot */
  xyz_simple_lowpass(ts, bp->diff_tau, &bst->mag_b_dot_f, &bst->mag_b_dot);

  /* Set coils to damp the rotation */
  xyz_scale(&act->mag_coil, bp->kp, &bst->mag_b_dot_f);

  return 0;
}

int
controller_init(void) {
  /* initialize bdot controller */
  return bdot_init();
}

int
run_controller(double ts, const measurements *meas, actuators *act)
{
  /* bdot only for now */
  return bdot(ts, meas, &p.bdot, &st.bdot, act);
}
