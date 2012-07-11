#ifndef AERO_TABLE_H
#define AERO_TABLE_H

#include "xyz.h"




typedef struct {
  xyz_t force;          // in body frame, Newtons, normalized to unit dynamic pressure
  xyz_t torque;         // likewise (Newton meters)
} aero_table_output_t;

typedef struct {
  double step;
  int dim;
  aero_table_output_t *table;
} aero_table_t;

// attsim should just care about these two functions:
int open_aero_table(const char *filename);
int get_aero_torque(xyz_t *torque, xyz_t *velocity_body, double atmosph_density);



  
int write_aero_table_entry(FILE *f, double a, double b, int c, aero_table_output_t *result); 

int load_aero_table(aero_table_t *t, const char *filename);
void close_aero_table(aero_table_t *t);

int lookup_aero_table(aero_table_output_t *out, aero_table_t *t, double a, double b, int c);

#endif
