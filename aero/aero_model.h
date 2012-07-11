#ifndef AERO_MODEL_H
#define AERO_MODEL_H

#include <xyz.h>


typedef struct {
  xyz_t origin;
  xyz_t side_a;   // side_a cross side_b should point "out"
  xyz_t side_b;
} panel_t;

typedef struct {
  int n_panels;
  panel_t *panels;
  xyz_t cg;
} aero_model_t;

typedef struct {
  double rho;     // density
  double v;    // scalar spacecraft velocity
  xyz_t e_v_ref;  // spacecraft velocity unit vector in ref frame (e.g. {1,0,0})
} aero_env_t;


int load_aero_model(aero_model_t *model, const char *filename);
int load_aero_env(aero_env_t *env, const char *filename);


#endif
