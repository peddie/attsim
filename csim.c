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

#define NUM_TOL 1e-9
#define SYS_SIZE 9

static void
print_gnuplot_log(const double t, const double y[],
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

int
main(int argc, char **argv)
{
  int status;
  double y[SYS_SIZE];
  dynamics_params dp;

  /* Integration parameters */
  double h = 1e-6, t = 0.0, t1 = 10;
  double tolabs = 1e-10, tolrel = 1e-10;

  const char *logfilename = "sim.log";
  
  FILE *simlog = fopen(logfilename, "w");
  
  if (argc > 1)
    t1 = atof(argv[1]);

  /* set up dynamics */
  if ((status = dynamics_init(y, &dp)) < 0) {
    fprintf(stderr, __FILE__ ":%d: ERROR: Couldn't set up initial conditions"
            " for dynamics simulation!\n", __LINE__);
    return 1;
  }

  /* GSL setup */
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
     
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, SYS_SIZE);
  gsl_odeiv2_control * c =
      gsl_odeiv2_control_standard_new(tolabs, tolrel, 1, 1);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(SYS_SIZE);

  if (!s || !c || !e) {
    fprintf(stderr, __FILE__ ":%d: Error in GSL initialization!\n", __LINE__);
    return 1;
  }

  gsl_odeiv2_system sys = {attitude, NULL, SYS_SIZE, &dp};

  printf("SIM:  Simulating attitude dynamics for %3.3lf seconds.\n", t1);
  printf("SIM:  Using (absolute, relative) tolerance of (%2.1e, %2.1e)\n",
         tolabs, tolrel);
  printf("SIM:  Initial state is as follows:\n"
         "\t\t<|q_i2b | %lf %lf %lf %lf>\n"
         "\t\t<|w_bi_b| %lf %lf %lf>\n"
         "\t\tmomentum: %lf\n"
         "\t\tkinetic energy: %lf\n",
         y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]);

  double time;
  clock_t before = clock();

  /* Step using GSL integrator */
  while (t < t1) {
    if ((status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, y))
        != GSL_SUCCESS)
      break;

    /* Check quaternion norm */
    double qnorm = y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
    if (fabs(qnorm - 1) > 1e-5)
      fprintf(stderr, "SIM:  Warning: Quaternion norm differs from 1: %lf\n",
              qnorm);

    /* output to stdout for now */
    print_gnuplot_log(t, y, &dp, simlog);
  }

  time = (double) clock() - (double) before;
  time /= CLOCKS_PER_SEC;

  printf("SIM:  Completed %3.3lf seconds of simulation.\n", t);
  printf("SIM:  Simulation took %3.3lf seconds (%3.2lf times realtime).\n",
         time, t / time);
  printf("SIM:  State output was written to file ``%s''.\n", logfilename);

  /* GSL cleanup */
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  fclose(simlog);
  printf("SIM:  Goodbye.\n");
  return 0;
}
