/**
   Common header stuff
 */
#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <math.h>

// a point in space
typedef struct {
  double x , y ;
} point ;

// struct to keep everyting clean
typedef struct {
  double (*step_fwd)( double (*f)( const point P ) ,
		      const double h , 
		      const point P ) ;
  double adaptive_growth , adaptive_shrink , 
    adaptive_errcon , adaptive_safety ;
  double error ;
  double tolerance ;
  size_t nsteps , notoksteps , maxsteps ;
} Integrator ;

// enum that lists all of the integration schemes we support
typedef enum {
  EULER , BACKWARD_EULER , MIDPOINT , IMPLICIT_MIDPOINT , HEUN , RALSTON , 
  RK3 , RK4 , RK4_38 , GAUSS4 , GAUSS5 , RADAU3 , RADAU5 
} integration_scheme ;

#endif
