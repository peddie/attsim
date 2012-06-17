#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>

#include <xyz.h>
#include <spatial_rotations.h>

#include "dynamics.h"
#include "sensors.h"
#include "controller.h"

#define NUM_TOL 1e-9

#define PERIODIC_POINTS_PER_SECOND 5

static void
print_sim_log(const double t, const double y[],
              const dynamics_params *dp, FILE *out)
{
  full_state s;
  /* Emit time */
  fprintf(out, "%2.8lf\t", t);

  /* Emit basic states */
  quat_memcpy(&s.q_i2b, (quat_t *) y);
  fprintf(out, "%.8e %.8e %.8e %.8e\t",
         s.q_i2b.q0, s.q_i2b.q1, s.q_i2b.q2, s.q_i2b.q3);

  xyz_memcpy(&s.w_bi_b, (xyz_t *) &y[4]);
  fprintf(out, "%.8e %.8e %.8e\t", s.w_bi_b.x, s.w_bi_b.y, s.w_bi_b.z);

  /* Emit euler angles */
  euler_t eulers;
  euler321_of_quat(&eulers, &s.q_i2b);
  fprintf(out, "%.8e %.8e %.8e\t", eulers.roll, eulers.pitch, eulers.yaw);

  /* Check work and energy */
  double T = compute_kinetic_energy(&s, dp);
  fprintf(out, "%.8e %.8e %.8e\t", y[7], T, y[7] - T);

  /* Check momentum */
  xyz_t L_bi_b;
  compute_momentum(&s, dp, &L_bi_b);
  fprintf(out, "%.8e %.8e\t", xyz_norm(&L_bi_b), y[8]);
  fprintf(out, "%.8e %.8e %.8e\t", L_bi_b.x, L_bi_b.y, L_bi_b.z);

  fprintf(out, "\n");
}

static void
print_ctrl_log(const double t, const actuators *act, FILE *out)
{
  fprintf(out, "%2.8lf\t", t);
  fprintf(out, "%.8e %.8e %.8e\t",
          act->mag_coil.x, act->mag_coil.y, act->mag_coil.z);

  fprintf(out, "\n");
}

static int
integrate_periodic_output(double tmax, dynamics_params *dp,
                          FILE *simlog, FILE *ctrllog,
                          gsl_odeiv2_system *sys, double reltol, double abstol,
                          const gsl_odeiv2_step_type *step, double y[SYS_SIZE])
{
  int i;
  double qerr = 1.0;
  double t = 0.0;
  gsl_odeiv2_driver *d
      = gsl_odeiv2_driver_alloc_standard_new(sys, step, 1e-6,
                                             abstol, reltol, 1, 1);
  if (!d) {
    fprintf(stderr, __FILE__ ":%d: Error in GSL initialization!\n", __LINE__);
    return 1;
  }

  fprintf(stderr, "SIM:  Using periodic output driver due to long simulation\n"
          "      time.  Configure the output time resolution with\n"
          "      PERIODIC_POINTS_PER_SECOND (currently %d).\n",
          PERIODIC_POINTS_PER_SECOND);

  /* Step using GSL integrator */
  for (i = 1; i <= tmax * PERIODIC_POINTS_PER_SECOND; i++) {
    double ti = (double) i / PERIODIC_POINTS_PER_SECOND;
    int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
     
    if (status != GSL_SUCCESS) {
      fprintf(stderr, "SIM:  GSL integration error"
              ", return value = %d\n", status);
      break;
    }

    /* Check quaternion norm */
    double qnorm = y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
    if ((fabs(qnorm - 1) > 1e-5) && (fabs(qnorm - qerr) > 1e-4)) {
      fprintf(stderr, "SIM:  Warning: Quaternion norm differs from 1: %lf\n",
              qnorm);
      qerr = qnorm;
    }

    /* Create measurements (passthrough for now) */
    measurements meas;
    full_state state;
    state.t = t;
    quat_memcpy(&state.q_i2b, (quat_t *) y);
    xyz_memcpy(&state.w_bi_b, (xyz_t *) &y[4]);
    state.T = y[7];
    state.P = y[8];
    
    if ((status = sensor_outputs(&state, dp, &meas)) < 0) {
      fprintf(stderr, "SIM:  Error in sensor output computation; exiting.\n");
      return 1;
    }
    
/*     xyz_memcpy(&meas.pqr, (xyz_t *) &y[4]); */
/*     xyz_set_all(&meas.mag_b, 0); */
/*     meas.panel_voltage[0] = meas.panel_voltage[1] = 0; */
    
    /* Run controller */
    if ((status = run_controller(((double) i+1) / PERIODIC_POINTS_PER_SECOND,
                                 &meas, &dp->act)) < 0) {
      fprintf(stderr, "SIM:  Error in controller computation; exiting.\n");
      return 1;
    }

    /* output to log file */
    print_sim_log(t, y, dp, simlog);
    print_ctrl_log(t, &dp->act, ctrllog);
  }
     
  gsl_odeiv2_driver_free(d);
  return 0;
}

static int
integrate_solver_output(double tmax, dynamics_params *dp,
                        FILE *simlog, FILE *ctrllog __attribute__((unused)),
                        gsl_odeiv2_system *sys, double abstol, double reltol,
                        const gsl_odeiv2_step_type *step, double y[SYS_SIZE])
{
  int status;
  double t = 0.0;
  double h = 1e-6;
  double qerr = 1.0;
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(step, SYS_SIZE);
  gsl_odeiv2_control * c =
      gsl_odeiv2_control_standard_new(abstol, reltol, 1, 1);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(SYS_SIZE);

  if (!s || !c || !e) {
    fprintf(stderr, __FILE__ ":%d: Error in GSL initialization!\n", __LINE__);
    return 1;
  }

  fprintf(stderr, "SIM:  Using direct solver output driver due to short "
          "simulation time.\n");

  /* Step using GSL integrator */
  while (t < tmax) {
    if ((status = gsl_odeiv2_evolve_apply(e, c, s, sys, &t, tmax, &h, y))
        != GSL_SUCCESS)
      break;

    /* Check quaternion norm */
    double qnorm = y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
    if ((fabs(qnorm - 1) > 1e-6) && (fabs(qnorm - qerr) > 1e-4)) {
      fprintf(stderr, "SIM:  Warning: Quaternion norm differs from 1: %lf\n",
              qnorm);
      qerr = qnorm;
    }
    
    /* output to log file */
    print_sim_log(t, y, dp, simlog);
  }

  /* GSL cleanup */
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  return 0;
}

int
main(int argc, char **argv)
{
  int status;
  double y[SYS_SIZE];
  dynamics_params dp;

  /* Integration parameters */
  double tmax = 10;
  double tolabs = 1e-12, tolrel = 1e-12;

  const char *simlogfile = "sim.log";
  const char *ctrllogfile = "controller.log";
  
  FILE *simlog = fopen(simlogfile, "w");
  FILE *ctrllog = fopen(ctrllogfile, "w");
  
  if (argc > 1)
    tmax = atof(argv[1]);

  /* set up dynamics */
  if ((status = dynamics_init(y, &dp)) < 0) {
    fprintf(stderr, __FILE__ ":%d: ERROR: Couldn't set up initial conditions"
            " for dynamics simulation!\n", __LINE__);
    return 1;
  }

  /* set up controller */
  if ((status = controller_init()) < 0) {
    fprintf(stderr, __FILE__":%d: ERROR: Couldn't initialize controller!\n",
            __LINE__);
    return 1;
  }
  
  /* GSL setup */
  /* I tried msadams, but it's around 8 times slower. */
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_system sys = {attitude, NULL, SYS_SIZE, &dp};

  printf("SIM:  Simulating attitude dynamics for %3.3lf seconds.\n", tmax);
  printf("SIM:  Using (absolute, relative) tolerance of (%2.1e, %2.1e)\n",
         tolabs, tolrel);
  printf("SIM:  Initial state is as follows:\n"
         "\t\t<|q_i2b | %lf %lf %lf %lf>\n"
         "\t\t<|w_bi_b| %lf %lf %lf>\n"
         "\t\tmomentum: %lf\n"
         "\t\tkinetic energy: %lf\n",
         y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]);

  double runtime;
  clock_t before = clock();

  if (tmax > 1)
    integrate_periodic_output(tmax, &dp, simlog, ctrllog,
                              &sys, tolabs, tolrel, T, y);
  else
    integrate_solver_output(tmax, &dp, simlog, ctrllog,
                            &sys, tolabs, tolrel, T, y);

  runtime = (double) clock() - (double) before;
  runtime /= CLOCKS_PER_SEC;

  printf("SIM:  Completed %3.3lf seconds of simulation.\n", tmax);
  printf("SIM:  Simulation took %3.3lf seconds (%3.2lf times realtime).\n",
         runtime, tmax / runtime);
  printf("SIM:  State output was written to file ``%s''.\n", simlogfile);

  fclose(simlog);
  printf("SIM:  Goodbye.\n");
  return 0;
}
