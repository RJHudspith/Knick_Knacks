#ifndef ISING_H
#define ISING_H

// some squirreled includes
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// Number of dimensions
#define ND 2

// where did the 1 go?
typedef enum { 
  RAWDATA = 0 ,
  JACKDATA = 2 ,
  BOOTDATA = 3 } resample_type ;

// struct containing our statistics
struct resampled {
  double *resampled ;
  double avg ;
  double err_hi ;
  double err_lo ;
  double err ;
  int NSAMPLES ;
  resample_type restype ;
} ;

// general struct 
struct general {
  int dims[ ND ] ;
  int vol ;
  int MAX_ITERS ; //10000
  int NMEAS ; // 100
  int THERM ; // 2000
} ;

// everyone sees this
extern struct general Latt ;

// handy little definition
#define LVOLUME Latt.vol

#endif
