#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "aero_table.h"

typedef struct {
  double a;
  double b;
  int32_t c;
  aero_table_output_t result;
} __attribute__((__packed__)) aero_table_entry_t;


void result_mac(aero_table_output_t *sum, aero_table_output_t *val, double weight);

int check_presence(aero_table_t *t, unsigned int check_a, unsigned int check_b, int c);
  
int write_aero_table_entry(FILE *f, double a, double b, int c, aero_table_output_t *result) {
  aero_table_entry_t entry;
  entry.a = a;
  entry.b = b;
  entry.c = c;
  entry.result = *result;  

  return (1 != fwrite(&entry, sizeof(entry), 1, f));
}

#define INDEX(a,b,c,dim)((unsigned int)(a) + (unsigned int)(b) * (dim) + (unsigned int)(c) * (dim) * (dim))

int load_aero_table(aero_table_t *t, char *filename) {
  FILE *f = fopen(filename, "rb");
  if (!f) {
    fprintf(stderr, "Error opening %s\n",filename);
    return 1;
  }

  aero_table_entry_t r[2];

  if (fread(r,sizeof(aero_table_entry_t),2,f) != 2) {
    fprintf(stderr,"Error reading first 2 entries from %s\n",filename);
    return 1;
  }
  

  printf("r[0].a,b = %.3f, %.2f\n",r[0].a,r[0].b);
  printf("r[1].a,b = %.3f, %.2f\n",r[1].a,r[1].b);


  // Measure the interval between entries
  if (r[1].a == r[0].a)
    t->step = r[1].b-r[0].b;
  else
    t->step = r[1].a - r[0].a;


  printf("Step size is %.3f\n",t->step);
  
  // Calculate how many entries per row
  t->dim = (round(2/t->step)+1);
  
  size_t s = t->dim * t->dim * 2 * sizeof(aero_table_output_t);

  printf("Dim is %u @ %zu bytes, about to allocate %.3f MB\n",t->dim, sizeof(aero_table_output_t), s/1E6);
  
  t->table = malloc(s);
  if (!t->table)
    fprintf(stderr, "Unable to allocate %.3f MB for table\n",s/1E6);

  memset(t->table, 0, s);

  rewind(f);

  unsigned int n=0;

  while (fread(r,sizeof(aero_table_entry_t),1,f) == 1) {
    n++;
    
    if (r->a < -1 || r->a > 1 || r->b < -1 || r->b > 1 || r->c < 0 || r->c > 1) {
      fprintf(stderr, "Rejecting invalid entry %u (%.3f, %.3f, %d)",n, r->a, r->b, r->c);
      continue;
    }

//    printf("Got entry for %.3f, %.3f, %d\n", r->a, r->b, r->c);

    unsigned int ia, ib, index;

    ia = round(r->a / t->step) + (t->dim - 1)/2;
    ib = round(r->b / t->step) + (t->dim - 1)/2;

    index = INDEX(ia, ib, r->c, t->dim);
//    printf("Putting it in %u, %u, %u (=%u)\n", ia, ib, r->c, index);
    
    t->table[index] = r->result;


  }

  fclose(f);

  printf("Read %u entries from %s.\n",n,filename);

  return 0;
}

void result_mac(aero_table_output_t *sum, aero_table_output_t *val, double weight) {
  xyz_t temp;
  
  xyz_scale(&temp, weight, &val->force);
  xyz_sum(&sum->force, &sum->force, &temp);
  
  xyz_scale(&temp, weight, &val->torque);
  xyz_sum(&sum->torque, &sum->torque, &temp);
}

int check_presence(aero_table_t *t, unsigned int check_a, unsigned int check_b, int c) {

  xyz_t r = t->table[INDEX(check_a, check_b, c, t->dim)].force;
  if (r.x == 0 && r.y == 0 && r.z == 0) return 0;
  return 1;
}


int lookup_aero_table(aero_table_output_t *out, aero_table_t *t, double a, double b, int c) {
  // Bilinear interpolation

  double ia, ib;
 
  // Convert coordinates to index form
  ia = a / t->step + (t->dim - 1)/2;
  ib = b / t->step + (t->dim - 1)/2;
 
  if (c<0 || c > 1) {
    fprintf(stderr,"lookup_aero_table called with bad c (%.3f, %.3f, %d)\n",a,b,c);
    return 1;
  }



  // Check we have enough neighboring points
  int n=0;
  for (unsigned int check_a = floor(ia); check_a <= ceil(ia); check_a++)
    for (unsigned int check_b = floor(ib); check_b <= ceil(ib); check_b++) {
      if (check_presence(t, check_a, check_b, c)) {
        *out = t->table[INDEX(check_a, check_b, c, t->dim)];
        n++;
      }
    }
  
  if (n==0) {
    fprintf(stderr,"lookup_aero_table couldn't find any neighbors for (%.3f, %.3f, %d): %f, %f. \n",a,b,c, ia, ib);
    return 1;
  }
  
//  if (n==1)
//    return; //with the only one we found

  if (n<4) {
  /*  printf("lookup_aero_table found only %d neighbors: ",n);
    for (unsigned int check_a = floor(ia); check_a <= ceil(ia); check_a++)
      for (unsigned int check_b = floor(ib); check_b <= ceil(ib); check_b++) {
        printf("(%d, %d)=%d  ", check_a, check_b, check_presence(t, check_a, check_b, c));
      }
    printf("\n");
    */
    return 0; //with some random one
  }

  // Interpolate along a
  aero_table_output_t r = {{0,0,0},{0,0,0}};
  aero_table_output_t r0,r1;
  r0=r1=r;

  if (floor(ia) == ceil(ia)) {
    r0 = t->table[INDEX(floor(ia),floor(ib),c,t->dim)];
    r1 = t->table[INDEX(floor(ia),ceil(ib),c,t->dim)];
  } else {   
    result_mac(&r0, &t->table[INDEX(floor(ia),floor(ib),c,t->dim)], ceil(ia)-ia);
    result_mac(&r0, &t->table[INDEX(ceil(ia), floor(ib),c,t->dim)], ia-floor(ia));
    
    result_mac(&r1, &t->table[INDEX(floor(ia),ceil(ib),c,t->dim)], ceil(ia)-ia);
    result_mac(&r1, &t->table[INDEX(ceil(ia), ceil(ib),c,t->dim)], ia-floor(ia));
  }
  
  // Interpolate along b
  if (floor(ib) == ceil(ib)) {
    r = r0;
  } else {   
    result_mac(&r, &r0, ceil(ib)-ib);
    result_mac(&r, &r1, ib-floor(ib));
  }
  
  *out = r;

//  printf(" looked up (%.3f, %.3f, %d) and found: F= %8.5f %8.5f %8.5f    T= %8.5f %8.5f %8.5f\n",
//      a,b,c, r.force.x, r.force.y, r.force.z, r.torque.x, r.torque.y, r.torque.z);
  return 0;
}



void close_aero_table(aero_table_t *t) {
  free(t->table);
}
