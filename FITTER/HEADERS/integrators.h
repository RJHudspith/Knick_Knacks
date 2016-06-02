/**
   @file integrators.h
   @brief numerical integrators
 */
#ifndef INTEGRATORS_H
#define INTEGRATORS_H

double
adaptive_simpsons( double (*f)( const double ) ,
 		   const double low , 
		   const double upp ,
		   const double eps ) ;

double
simpsons_arr5( const double *y , 
	       const double *x ,
	       const int N ) ;

#endif
