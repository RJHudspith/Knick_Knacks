/**
   @file conformal_amu.c
   @brief conformal polynomial definition of a_\mu
 */
#include "fitfunc.h"

#include "equivalents.h"
#include "fit_and_plot.h"
#include "padefit.h"
#include "integrators.h"
#include "tmoment_routines.h"
#include "Utils.h"
#include "polyfit.h"

#define FOURMPI2 0.55

// pass this into the integrator
static int NPARAMS = 0 ;
static double *poly_params = NULL ;

static inline double
conformal_w( const double Q2 )
{
  return ( 1.0 - sqrt( 1 + Q2 / FOURMPI2 ) ) / ( 1.0 + sqrt( 1 + Q2 / FOURMPI2 ) ) ;
}

// evaluate f(q^2)
static double
evaluate_fq_squared( const double Q2 )
{
  if( Q2 < 1E-32 ) return 0.0 ;
  const double musq = 0.1056583715 * 0.1056583715 ; 
  const double ZQ = -( 1.0 - sqrt( 1.0 + 4. * musq / Q2 ) ) / ( 2.0 * musq ) ;
  return musq * Q2 * pow( ZQ , 3.0 ) * ( 1.0 - Q2 * ZQ ) / ( 1.0 + musq * Q2 * ZQ * ZQ ) ;
}

// compute the integrand, x is supposed to be momentum^2
static double
evaluate_integrand_HVP( const double x )
{
  struct x_descriptor X ; X.X = conformal_w( x ) ;
  return poly_eval( poly_params , X , NPARAMS ) * evaluate_fq_squared( x ) ;
}

// numerically integrate pade
static void
integrate_poly( struct resampled *INT ,
		const struct resampled *POLYS ,
		const double **x ,
		const struct input_params *INPARAMS ,
		const int NSLICES )
{
  // numerically integrate pade using the fit expression
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    // globalise the scale
    equate_constant( &INT[i] , 0.0 , POLYS[0].NSAMPLES , POLYS[0].restype ) ;

    // find the matchpoint
    const double Q2 = -FOURMPI2 * 4.0 * INPARAMS -> fit_hi / pow( INPARAMS -> fit_hi + 1. , 2 ) ;

    // x is physical
    const int lopos = find_idx( Q2 , x[i] , INPARAMS->NDATA[i] , 0 ) ;
    printf( "[PADE AMU] Head LOPOS %d :: %f %f \n" , lopos , x[i][lopos] , conformal_w( x[i][lopos] ) ) ;

    // and integrate the boots
    size_t k , j ;
    for( j = 0 ; j < POLYS[ 0 ].NSAMPLES ; j++ ) {
      for( k = 0 ; k < NPARAMS ; k++ ) {
	poly_params[ k ] = POLYS[ k ].resampled[ j ] ;
      }
      INT[ i ].resampled[ j ] = adaptive_simpsons( evaluate_integrand_HVP ,
						   0.0 , x[i][lopos] , 1E-12 ) ;
    }
    // do the average
    for( k = 0 ; k < NPARAMS ; k++ ) {
      poly_params[ k ] = POLYS[ k ].avg ;
    }
    INT[ i ].avg = adaptive_simpsons( evaluate_integrand_HVP ,
    				      0.0 , x[i][lopos] , 1E-12 ) ;
    compute_err( &INT[ i ] ) ;
  }
  return ;
}

// integrate the tail
static void
integrate_tail( struct resampled *INT ,
		const struct resampled **PIP ,
		const double **x ,
		const struct mom_info **momenta ,
		struct input_params *INPARAMS ,
		const int NSLICES )
{
  // compute the integration
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    equate_constant( &INT[i] , 0.0 , PIP[i][0].NSAMPLES , PIP[i][0].restype ) ;

    const double Q2 = -FOURMPI2 * 4.0 * INPARAMS -> fit_hi / pow( INPARAMS -> fit_hi + 1. , 2 ) ;
    const int lopos = find_idx( Q2 , x[i] , INPARAMS->NDATA[i] , 0 ) ;

    // tell us where the matching point is
    printf( "[CONPOLY] Tail LOPOS %d %f \n" , lopos , x[ i ][lopos] ) ;
    const size_t newmax = INPARAMS -> NDATA[i] - lopos ;

    // set up the x-axis
    double xloc[ newmax ] ;
    size_t k ;
    for( k = 0 ; k < newmax ; k++ ) {
      xloc[ k ] = x[ i ][ lopos + k ] ;
    }
    // loop boots
    #pragma omp parallel for private(k)
    for( k = 0 ; k < INT[i].NSAMPLES ; k++ ) {
      double yloc[ newmax ] ;
      size_t j ;
      for( j = 0 ; j < newmax ; j++ ) {
	yloc[ j ] = PIP[ i ][ lopos + j ].resampled[ k ] * evaluate_fq_squared( xloc[ j ] ) ;
      }
      INT[i].resampled[k] = simpsons_arr5( yloc , xloc , newmax ) ;
    }
    double yloc[ newmax ] ;
    for( k = 0 ; k < newmax ; k++ ) {
      yloc[ k ] = PIP[ i ][ lopos + k ].avg * evaluate_fq_squared( xloc[ k ] ) ;
    }
    INT[i].avg = simpsons_arr5( yloc , xloc , newmax ) ;
    compute_err( &INT[i] ) ;
  }
  return ;
}

// time moment eval
void
conformal_amu_eval( double **xavg ,
		    struct resampled **bootavg ,
		    struct mom_info **mominfo ,
		    struct input_params *INPARAMS ,
		    const int NSLICES ,
		    const int LT ,
		    const bool renormalise )
{
 switch( INPARAMS->fittype ) {
  case POLY0 : NPARAMS = 1 ; break ;
  case POLY1 : NPARAMS = 1 ; break ;
  case POLY2 : NPARAMS = 2 ; break ;
  case POLY3 : NPARAMS = 3 ; break ;
  case POLY4 : NPARAMS = 4 ; break ;
  case POLY5 : NPARAMS = 5 ; break ;
  case POLY6 : NPARAMS = 6 ; break ;
  case POLY7 : NPARAMS = 7 ; break ;
  case POLY8 : NPARAMS = 8 ; break ;
  case POLY9 : NPARAMS = 9 ; break ;
  default :
    printf( "[CONF_AMU] poly not recognised, leaving \n" ) ;
    return ;
  }

  // renormalise tcorr
  printf( "\n--> Renormalising <--\n" ) ;

  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    const double ainv2 = INPARAMS -> quarks[i].ainverse *
      INPARAMS -> quarks[i].ainverse ;
    size_t j ;
    for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
      mult_constant( &bootavg[i][j] , INPARAMS -> quarks[i].ZV ) ;
      xavg[i][j] = conformal_w( ainv2 * mominfo[i][j].p2 ) ;
    }
  }

  // fit a polynomial in this
  struct resampled *POLY = fit_data_plot_data( (const struct resampled**)bootavg , 
					       (const double**)xavg , 
					       (const struct mom_info**)mominfo ,
					       *INPARAMS , NSLICES , LT ) ;

  // subtract the fitted 0 and make momenta physical
  for( i = 0 ; i < NSLICES ; i++ ) {
    const double ainv2 = INPARAMS -> quarks[i].ainverse *
      INPARAMS -> quarks[i].ainverse ;
    size_t j ;
    for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
      subtract( &bootavg[i][j] , POLY[0] ) ;
      xavg[i][j] = mominfo[i][j].p2 * ainv2 ;
    }
    // set pade[0] to zero
    equate_constant( &POLY[ 0 ] , 0.0 , POLY[ 0 ].NSAMPLES , POLY[ 0 ].restype ) ;
  }

  // allocate "poly params" array
  poly_params = malloc( ( NPARAMS ) * sizeof( double ) ) ;

  // numerically integrate the top part of the data
  struct resampled *HEAD = malloc( NSLICES * sizeof( struct resampled ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    HEAD[i].resampled = malloc( POLY[0].NSAMPLES * sizeof( double ) ) ;
  }
  integrate_poly( HEAD , POLY , (const double**)xavg , INPARAMS , NSLICES ) ;

  // numerically integrate the smooth tail
  struct resampled *TAIL = malloc( NSLICES * sizeof( struct resampled ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    TAIL[i].resampled = malloc( bootavg[i][0].NSAMPLES * sizeof( double ) ) ;
  }
  integrate_tail( TAIL , (const struct resampled **)bootavg ,
		  (const double **)xavg , 
		  (const struct mom_info **)mominfo ,
		  INPARAMS , NSLICES ) ;

  for( i = 0 ; i < NSLICES ; i++ ) {
    printf( "[PADE_AMU] HEAD :: %e %e\n" , HEAD[ i ].avg , HEAD[i].err ) ;
    printf( "[PADE_AMU] TAIL :: %e %e\n" , TAIL[ i ].avg , TAIL[i].err ) ;
    add( &HEAD[i] , TAIL[i] ) ;
    mult_constant( &HEAD[ i ] , 4.0 / ( 137.03599878 * 137.03599878 ) ) ; 
    printf( "[PADE_AMU] amu :: %e %e\n" , HEAD[ i ].avg , HEAD[i].err ) ;
  }

  free( poly_params ) ;
  free_resampled( HEAD , NSLICES ) ;
  free_resampled( TAIL , NSLICES ) ;

  return ;
}
