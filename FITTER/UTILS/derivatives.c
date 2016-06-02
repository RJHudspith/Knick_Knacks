/**
   @file derivatives.c
   @brief dirty numerical derivative
 */
#include "fitfunc.h"

static double
symmetric_der( const fitfunc fit ,  
	       const double *data , 
	       const double extrap_point , 
	       const struct chiral quarks , 
	       const struct mom_info mominfo ,
	       const int LT ,
	       const int SLICE )
{
  const double eps = 1E-12 ;
  size_t iters = 0 ;
  double dh = 0.01 , der1 = 1.0 , der2 = 0.0 ;
  while( fabs( der1 - der2 ) > eps || iters < 5 ) {

    // first one
    struct x_descriptor XX1 = { extrap_point + dh , quarks , mominfo , LT } ;
    der1 = fit.f( data , XX1 , SLICE ) ;
    struct x_descriptor XX2 = { extrap_point - dh , quarks , mominfo , LT } ;
    der1 -= fit.f( data , XX2 , SLICE ) ;
    der1 /= ( 2 * dh ) ;
    
    // second one
    struct x_descriptor XX3 = { extrap_point + dh/2. , quarks , mominfo , LT } ;
    der2 = fit.f( data , XX3 , SLICE ) ;
    struct x_descriptor XX4 = { extrap_point - dh/2. , quarks , mominfo , LT } ;
    der2 -= fit.f( data , XX4 , SLICE ) ;
    der2 /= ( dh ) ;
    
    // make h smaller
    dh *= 2.0 * sqrt( 2 * fabs( der1 - der2 ) ) ;
    if( dh < eps ) break ;
    iters++ ;
  }
  return der2 ;
}

// evaluate derivative at point "x" using a stencil
struct resampled
fit_der( const struct resampled *fitparams ,
	 const double extrap_point ,
	 const int NPARAMS ,
	 const struct chiral quarks ,
	 const struct mom_info mominfo ,
	 const fitfunc fit ,
	 const int LT ,
	 const int SLICE )
{
  struct resampled y ;
  y.resampled = malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
  y.NSAMPLES  = fitparams[0].NSAMPLES ;
  y.restype   = fitparams[0].restype ;

  size_t j ;
  #pragma omp parallel for private(j)
  for( j = 0 ; j < y.NSAMPLES ; j++ ) {

    double data[ NPARAMS ] ;
    // switch to individual axis
    int k ;
    for( k = 0 ; k < NPARAMS ; k++ ) {
      data[ k ] = fitparams[ k ].resampled[ j ] ;
    }
    
    y.resampled[j] = symmetric_der( fit , data , extrap_point , 
				    quarks , mominfo , LT , SLICE ) ;
  }

  double data[ NPARAMS ] ;
  // switch to individual axis
  int k ;
  for( k = 0 ; k < NPARAMS ; k++ ) {
    data[ k ] = fitparams[ k ].avg ;
  }
  y.avg = symmetric_der( fit , data , extrap_point , 
			 quarks , mominfo , LT , SLICE ) ;

  compute_err( &y ) ;
  printf( "%f %f %f \n" , extrap_point , y.avg , y.err ) ;
  return y ;
}
