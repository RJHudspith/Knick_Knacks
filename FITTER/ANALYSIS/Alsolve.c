/**
   @file Alsolve.c
   @brief solve for alpha
 */
#include "fitfunc.h"
#include "equivalents.h"
#include "cylinder.h"
#include "momentum_stuff.h"

// solves PERT - DATA == 0 for alpha using Newton-Raphson
struct resampled
alpha_solve( const struct resampled data ,
	     const double t1 ,
	     const double t2 )
{
  const double term1 = 1.6398212 - 1.125 * ( t1 + t2 ) ;
  const double term2 = 49.0757 
    + t1 * ( -0.8382151 + t1 * ( 0.40025907 - 0.07213198 * t1 ) ) 
    + t2 * ( -0.8382151 + t2 * ( 0.40025907 - 0.07213198 * t2 ) )
    + t1 * t2 * ( 0.40025907 - 0.07213198 * ( t1 + t2 ) ) ;
		  
  struct resampled alpha = init_dist( NULL , data.NSAMPLES , data.restype ) ;
  size_t j ;
  for( j = 0 ; j < data.NSAMPLES ; j++ ) {
    double a = 0.1 , eps = 1 ;
    int iters = 0 ;
    while( eps > 1E-14 && iters < 100 ) {

      const double fx = a + a*a*( term1 ) + a*a*a*( term2 ) -data.resampled[j] ;
      const double fxp = 1 + 2*a*term1 + 3*a*a*term2 ;
      const double fxpp = 2*term1 + 6*term2 ;

      const double up = a - 2*fx*fxp / ( 2*fxp*fxp - fx*fxpp ) ;
      eps = fabs( up - a ) ;
      a = up ;
      iters++ ;
    }
    if( iters >= 100 || eps > 1E-14 ) {
      printf( "Bad iteration %zu \n" , j ) ;
    }
    alpha.resampled[j] = a * M_PI ;
  }
  compute_err( &alpha ) ;
  //printf( "alpha :: %f +/- %f \n" , alpha.avg , alpha.err ) ;
  return alpha ;
}

// Alpha_s computation
void
alphas_solve_eval( double **xavg ,
		   struct resampled **bootavg ,
		   struct mom_info **mominfo ,
		   struct input_params *INPARAMS ,
		   const int NSLICES ,
		   const int LT ,
		   const bool renormalise )
{
  printf( "\n--> Recylinder <--\n" ) ;

  int j , i ;
  const double widths[3] = { 0.25 , 0.25 , 0.25 } ;
  recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , LT , 4 , widths ) ;

  // BOOT and X are freed in momavg
  struct mom_info **mavg ;
  struct resampled **BAVG = momentum_average( &mavg , INPARAMS , 
					      (const struct mom_info**)mominfo , 
					      (const struct resampled**)bootavg , 
					      NSLICES , PSQ_AVERAGE ) ;
  // set up x-avg
  double **x = malloc( NSLICES * sizeof( double* ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    x[ i ] = ( double* )malloc( INPARAMS -> NDATA[i] * sizeof( double ) ) ;
    int k ;
    for( k = 0 ; k < INPARAMS -> NDATA[ i ] ; k++ ) {
      x[ i ][ k ] = mavg[ i ][ k ].p2 ;
    }
  }

  printf( "\n--> Converting to physical momenta <--\n" ) ;

  convert_to_physmom( mavg , x , INPARAMS , NSLICES ) ; 

  printf( "\n--> solving <--\n" ) ;

  // solve goes here
  for( j = 0 ; j < NSLICES ; j++ ) {

    // renormalise 
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      mult_constant( &BAVG[j][i] , INPARAMS -> quarks[j].ZV ) ;
    }
  }

  return ;
}
