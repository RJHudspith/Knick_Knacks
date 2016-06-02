#include "fitfunc.h"

#include "dft.h"
#include "fit_chooser.h"
#include "equivalents.h"
#include "graph_data.h"
#include "integrators.h"
#include "padefit.h"
#include "tmoment_routines.h"
#include "shuffalgo.h"
#include "Utils.h"

// pass this into the integrator
static double *pade_params = NULL , scalea = 1.0 ;

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
evaluate_integrand_HVP( const double Q2 )
{
  int n , m ;
  pade_get_nm( &n , &m ) ;
  struct x_descriptor X ; X.X = sqrt( Q2 ) / scalea ; // rescale these
  return pade_eval( pade_params , X , n + m + 1 ) * evaluate_fq_squared( Q2 ) ;
}

// numerically integrate pade
void
integrate_pade( struct resampled *INT ,
		const struct resampled **PADES ,
		const double **xavg ,
		const struct input_params *INPARAMS ,
		const int idx ,
		const int n , 
		const int m ,
		const size_t NSLICES , 
		const int matchpoint[ NSLICES ] )
{
  // numerically integrate pade using the fit expression
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    const double ainv2 = pow( INPARAMS -> quarks[ i ].ainverse , 2.0 )  ;
    // globalise the scale
    scalea = INPARAMS -> quarks[i].ainverse ;
    equate_constant( &INT[i] , 0.0 , PADES[i][0].NSAMPLES , PADES[i][0].restype ) ;
    // x is physical
    printf( "[TMOMENTS] Head LOPOS %d :: %f \n" , matchpoint[i] , xavg[idx+i][matchpoint[i]] ) ;
    // and integrate the boots
    size_t k , j ;
    // do not parallelise!
    for( j = 0 ; j < PADES[ i ][ 0 ].NSAMPLES ; j++ ) {
      for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
	pade_params[ k ] = PADES[ i ][ k ].resampled[ j ] ;
      }
      INT[ i ].resampled[ j ] = adaptive_simpsons( evaluate_integrand_HVP ,
						   0.0 , xavg[idx+i][matchpoint[i]] * ainv2 , 1E-12 ) ;
    }
    // do the average
    for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
      pade_params[ k ] = PADES[ i ][ k ].avg ;
    }
    INT[ i ].avg = adaptive_simpsons( evaluate_integrand_HVP ,
    				      0.0 , xavg[idx+i][matchpoint[i]] * ainv2 , 1E-12 ) ;
    compute_err( &INT[ i ] ) ;
  }
  return ;
}

// integrate the tail
void
integrate_tail( struct resampled *INT ,
		const struct resampled **bootavg ,
		const double **xavg ,
		const struct mom_info **mominfo ,
		struct input_params *INPARAMS ,
		const int idx ,
		const size_t i ,
		const int matchpoint )
{
  // compute the integration
  equate_constant( &INT[i] , 0.0 , bootavg[idx][0].NSAMPLES , bootavg[idx][0].restype ) ;

  // a^-2
  const double ainv2 = pow( INPARAMS -> quarks[ i ].ainverse , 2.0 )  ;

  // tell us where the matching point is
  printf( "[TMOMENTS] Tail LOPOS %d %f \n" , matchpoint , xavg[ idx ][matchpoint] ) ;
  const size_t newmax = INPARAMS -> NDATA[idx] - matchpoint ;
  printf( "newmax %zu \n" , newmax ) ;

  // set up the x-axis
  double xloc[ newmax ] ;
  size_t k ;
  for( k = 0 ; k < newmax ; k++ ) {
    xloc[ k ] = xavg[ idx ][ matchpoint + k ] * ainv2 ;
  }
  // loop boots
#pragma omp parallel for private(k)
  for( k = 0 ; k < INT[0].NSAMPLES ; k++ ) {
    double yloc[ newmax ] ;
    size_t j ;
    for( j = 0 ; j < newmax ; j++ ) {
      yloc[ j ] = bootavg[ idx ][ matchpoint + j ].resampled[ k ] * evaluate_fq_squared( xloc[ j ] ) ;
    }
    INT[i].resampled[k] = simpsons_arr5( yloc , xloc , newmax ) ;
  }
  double yloc[ newmax ] ;
  for( k = 0 ; k < newmax ; k++ ) {
    yloc[ k ] = bootavg[ idx ][ matchpoint + k ].avg * evaluate_fq_squared( xloc[ k ] ) ;
  }
  INT[i].avg = simpsons_arr5( yloc , xloc , newmax ) ;
  compute_err( &INT[i] ) ;

  return ;
}

// print to stdout the integrand
void
inspect_integrand( const struct resampled **PADES ,
		   const double **x ,
		   const struct input_params *INPARAMS ,
		   const size_t NSLICES )
{
  printf( "--> Integrand <-- \n\n" ) ;
  int n , m ;
  pade_get_nm( &n , &m ) ;
  pade_params = malloc( ( n+m+1 ) * sizeof( double ) ) ;
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    scalea = INPARAMS -> quarks[i].ainverse ;
    size_t j ;
    for( j = 0 ; j < (n+m+1) ; j++ ) {
      pade_params[ j ] = PADES[ i ][ j ].avg ;
      printf( "[TMOMENTS] pade :: %zu %1.12f\n" , j , pade_params[ j ] ) ;
    }
    for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
      printf( "[TMOMENTS] integrand %e %e \n" , x[i][j] , evaluate_integrand_HVP( x[i][j] ) ) ;
    }
  }
  free( pade_params ) ;
  return ;
}

// subtract \PI(0) derived from the moments evaluation
void
subtract_PI0( struct resampled **PIP ,
	      struct resampled **PADES ,
	      double **x ,
	      const struct input_params *INPARAMS ,
	      const size_t idx ,
	      const size_t i )
{
  printf( "--> Subtracting <--\n\n" ) ;
  size_t j ;
  for( j = 0 ; j < INPARAMS->NDATA[idx] ; j++ ) {
    subtract( &PIP[idx][j] , PADES[i][0] ) ;
  }
  // set the constant PADE term to 0
  equate_constant( &PADES[i][ 0 ] , 0.0 , 
		   PADES[i][ 0 ].NSAMPLES ,
		   PADES[i][ 0 ].restype ) ;
  return ;
}

// sensible matchpoint
int
sensible_matchpoint( const struct resampled *PADES ,
		     const struct resampled *data ,
		     const double *xavg ,
		     const size_t NDATA ,
		     const double ainverse , 
		     const int n , 
		     const int m )
{
  // option 1 pade lies below HVP at some point
  size_t match , k ;
  for( match = 1 ; match < NDATA ; match++ ) {
    struct x_descriptor X ; X.X = xavg[match] ; //sqrt( xavg[match] ) / ainverse ;
    size_t j ;
    double pade_hi = 0 , pade_lo = 1000 ;
    for( j = 0 ; j < data[match].NSAMPLES ; j++ ) {
      for( k = 0 ; k < n+m+1 ; k++ ) {
	pade_params[ k ] = PADES[k].resampled[j] ;
      }
      const double peval = pade_eval( pade_params , X , n + m + 1 ) ;
      if( peval > pade_hi ) pade_hi = peval ;
      if( peval < pade_lo ) pade_lo = peval ; 
    }
    printf( "HI comp :: %f > %f \n" , pade_hi , data[match].err_lo ) ; 
    printf( "LO comp :: %f < %f \n" , pade_lo , data[match].err_hi ) ; 
    if( pade_hi < data[match].err_lo || pade_lo > data[match].err_hi ) return match-1 ;
  }
  return 0 ;
}

// time moment eval
void
time_moment_eval( double **xavg ,
		  struct resampled **bootavg ,
		  struct mom_info **mominfo ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT ,
		  const bool renormalise )
{
  printf( "\n--> Time Moment evaluation <--\n" ) ;

  // set n and m
  int NPARAMS ;
  switch( INPARAMS->fittype ) {
  case PADE11 : pade_set_nm( 1 , 1 ) ; NPARAMS = 3 ; break ;
  case PADE12 : pade_set_nm( 1 , 2 ) ; NPARAMS = 4 ; break ;
  case PADE21 : pade_set_nm( 2 , 1 ) ; NPARAMS = 4 ; break ;
  case PADE31 : pade_set_nm( 3 , 1 ) ; NPARAMS = 5 ; break ;
  case PADE22 : pade_set_nm( 2 , 2 ) ; NPARAMS = 5 ; break ;
  case PADE13 : pade_set_nm( 1 , 3 ) ; NPARAMS = 5 ; break ;
  case PADE41 : pade_set_nm( 4 , 1 ) ; NPARAMS = 6 ; break ;
  case PADE32 : pade_set_nm( 3 , 2 ) ; NPARAMS = 6 ; break ; 
  case PADE23 : pade_set_nm( 2 , 3 ) ; NPARAMS = 6 ; break ;
  case PADE14 : pade_set_nm( 1 , 4 ) ; NPARAMS = 6 ; break ;
  case PADE33 : pade_set_nm( 3 , 3 ) ; NPARAMS = 7 ; break ;
  default :
    printf( "[TMOMENTS] pade not recognised, leaving \n" ) ;
    return ;
  }

  // renormalise tcorr
  printf( "\n--> Renormalising <--\n" ) ;

  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    size_t j ;
    for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
      #ifdef verbose
      printf( "[TMOMENTS] corr :: %zu %e %e \n" , j , bootavg[i][j].avg , bootavg[i][j].err ) ;
      #endif
      mult_constant( &bootavg[i][j] , INPARAMS->quarks[0].ZV ) ;
    }
  }

  // n,m pade
  int n , m ;
  pade_get_nm( &n , &m ) ;
  printf( "[TMOMENTS] ( %d %d ) pade ( %d ) \n" , n , m , NPARAMS ) ;

  // NMAX needs to be greater than n+m+1
  const int NMAX = 12 ; //10 > ( n + m + 1 ) ? 10 : n + m + 1 ;

  // compute polynomial coefficients
  struct resampled **PIS = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    PIS[i] = compute_PIS( bootavg[i] , 
			  NMAX , INPARAMS->dimensions[i][3] , 
			  HPQCD_MOMENTS ,
			  INPARAMS->mom_type ) ;
  }

  // compute the pade coefficients
  struct resampled **PADES = compute_PADE( PIS , NSLICES/4 , n , m ) ;

  // print the average ...
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    double pades[ n + m + 1 ] ;
    size_t k ;
    for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
      pades[ k ] = PADES[ i ][ k ].avg ;
    }
    // pretty print the pade
    struct x_descriptor X ;
    pade_description( "[TMOMENTS] pade-eval" , pades , X , n+m+1 ) ;
  }

  make_xmgrace_graph( "comparison.agr" , "|aq|" , "\\xP\\f{}" ) ;
  printf( "NSLICES :: %d \n" , NSLICES ) ;

  // allocate "pade params" array
  pade_params = malloc( ( n + m + 1 ) * sizeof( double ) ) ;

  // subtract \Pi(0)
  int matches[ NSLICES/4 ] ;
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    const size_t idx = 2*NSLICES/4 + i ;
    subtract_PI0( bootavg , PADES , xavg , INPARAMS , idx , i ) ;
    matches[i] = sensible_matchpoint( PADES[i] , bootavg[idx] , xavg[idx] , INPARAMS->NDATA[idx] , INPARAMS->quarks[i].ainverse , n , m ) ;
  }

  // initialise the xmgrace graph
  {
    for( i = 0 ; i < NSLICES/4 ; i++ ) {
      const size_t idx = 2*NSLICES/4 + i ;
      // plot the data first but keep the graph file open
      plot_data( bootavg[idx] , xavg[idx] , INPARAMS->NDATA[idx] ) ;

      int n , m ;
      pade_get_nm( &n , &m ) ;
      fitfunc fit ;
      fit.f = pade_eval ;
      const double ainv2 = INPARAMS->quarks[i].ainverse * INPARAMS->quarks[i].ainverse ;
      // plot the data we use
      plot_fitfunc2( PADES[i] , bootavg[idx] , xavg[idx] , mominfo[idx] , fit , n+m+1 ,
		     INPARAMS->quarks[i] , 0.0 , 
		     xavg[idx][ matches[i] ] ,
		     INPARAMS->dimensions[i][3] , INPARAMS -> NDATA[idx] , n+m+1 ) ;
      // plot some more of it
      plot_fitfunc2( PADES[i] , bootavg[idx] , xavg[idx] , mominfo[idx] , fit , n+m+1 ,
		     INPARAMS->quarks[i] , 
		     xavg[idx][ matches[i] ]  ,
		     xavg[idx][ find_idx( sqrt( INPARAMS -> fit_hi / ainv2 ) , xavg[idx] , 
					  INPARAMS->NDATA[idx] , 0 ) ] ,
		     INPARAMS->dimensions[i][3] , INPARAMS -> NDATA[2*i] , n+m+1 ) ;
    }
  }
  close_xmgrace_graph(  ) ;

  // numerically integrate the top part of the data
  struct resampled *HEAD = malloc( NSLICES/4 * sizeof( struct resampled ) ) ;
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    HEAD[i].resampled = malloc( PADES[i][0].NSAMPLES * sizeof( double ) ) ;
  }
  integrate_pade( HEAD , (const struct resampled**)PADES , 
		  (const double**)xavg ,
		  INPARAMS , 2*NSLICES/4 , n , m , NSLICES/4 ,
		  matches ) ;

  // numerically integrate the smooth tail
  struct resampled *TAIL = malloc( NSLICES/4 * sizeof( struct resampled ) ) ;
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    const size_t idx = 2*NSLICES/4 + i ;
    TAIL[i].resampled = malloc( bootavg[idx][0].NSAMPLES * sizeof( double ) ) ;
    integrate_tail( TAIL , (const struct resampled **)bootavg ,
		    (const double **)xavg , 
		    (const struct mom_info **)mominfo ,
		    INPARAMS , idx , i , matches[i] ) ;
  }

  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    printf( "[TMOMENTS] HEAD :: %e %e\n" , HEAD[ i ].avg , HEAD[i].err ) ;
    printf( "[TMOMENTS] TAIL :: %e %e\n" , TAIL[ i ].avg , TAIL[i].err ) ;
    add( &HEAD[i] , TAIL[i] ) ;
    mult_constant( &HEAD[ i ] , 4.0 / ( 137.03599878 * 137.03599878 ) ) ; 
    printf( "[TMOMENTS] amu :: %e %e\n" , HEAD[ i ].avg , HEAD[i].err ) ;
  }
  free_resampled( HEAD , NSLICES/4 ) ;
  free_resampled( TAIL , NSLICES/4 ) ;

  // free the integrations
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    free_resampled( PADES[i] , n+m+1 ) ;
    free_resampled( PIS[i] , NMAX ) ;
  }
  free( PIS ) ;
  free( PADES ) ;
  free( pade_params ) ;

  return ;
}

