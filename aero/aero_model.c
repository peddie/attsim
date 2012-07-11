#include "aero_model.h"

#define UU __attribute__((unused))

int load_aero_model(aero_model_t *model, const char *filename UU) {

  // Geometry
  
#include "geom.c"

  model->n_panels = sizeof(panels_body)/sizeof(panel_t);
  model->panels = panels_body;

  return 0;
}

int load_aero_env(aero_env_t *env, const char *filename UU) {

  env->rho = 1.27E-11;
  env->v = 7700;

  return 0;
}

