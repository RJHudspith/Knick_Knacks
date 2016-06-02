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

// define PI0 somewhere
struct resampled *PI0 ;

// evaluate f(q^2)
static double
evaluate_fq_squared( const double Q2 )
{
  if( Q2 < 1E-32 ) return 0.0 ;
  const double musq = 0.1056583715 * 0.1056583715 ; 
  const double ZQ = -( 1.0 - sqrt( 1.0 + 4. * musq / Q2 ) ) / ( 2.0 * musq ) ;
  return musq * Q2 * pow( ZQ , 3.0 ) * ( 1.0 - Q2 * ZQ ) / ( 1.0 + musq * Q2 * ZQ * ZQ ) ;
}

// correlator length
static double ct[ 128 ] ;
static int Ltmp = 64 ;

// dft in some direction
double
arbq_indiv( const double Q2 )
{
  register double sum1 = 0.0 , sum2 = 0.0 , sum3 = 0.0 ;
  const double mom = 2.0 * asin( 0.5 * sqrt( Q2 ) ) ;
  int t ;
  for( t = 0 ; t < Ltmp ; t++ ) {
    sum1 += ct[t] * cos( mom * ( t < Ltmp/2 ? t : t-Ltmp ) ) ;
    sum2 += ct[t] ;
    sum3 += ct[t] * ( t < Ltmp/2 ? (double)t*t : (t-Ltmp)*(t-Ltmp) ) ;
  }
  return ( sum1 - sum2 ) / ( Q2 >0.0 ? Q2 : 1 ) + 0.5 * sum3 ;
}

// evaluate the integrand
double
evaluate_integrand_HVP( const double Q2 )
{
  return arbq_indiv( Q2 / ( 3.148 * 3.148 ) ) * evaluate_fq_squared( Q2 ) ;
}

// returns the DFT at some arbitrary q in one direction
static 
struct resampled arbq( double *smom ,
		       const struct resampled *corr , 
		       const double q ,
		       const int L )
{
  struct resampled res ;
  res.resampled = calloc( corr[0].NSAMPLES , sizeof( double ) ) ;
  res.NSAMPLES  = corr[0].NSAMPLES ;
  res.restype   = corr[0].restype ;
  // compute twiddle factors
  const double mom = q * 2.0 * M_PI / L ;
  *smom = 2.0 * sin( mom * 0.5 ) ;
  double coscache[ L ] ; 
  int t , j ;
  for( t = 0 ; t < L ; t++ ) {
    coscache[ t ] = cos( mom * ( t < L/2 ? t : -L + t ) ) ;
  }
  for( j = 0 ; j < res.NSAMPLES ; j++ ) {
    register double sum1 = 0.0 , sum2 = 0.0 , ct ;
    for( t = 0 ; t < L ; t++ ) {
      ct = corr[t].resampled[j] ;
      sum1 += ct * coscache[ t ] ;
      sum2 += ct ;
    }
    res.resampled[j] = ( sum1 - sum2 ) / ( q != 0.0 ? (*smom)*(*smom) : 1 ) ;
  }
  compute_err( &res ) ;
  return res ;
} 

// integrate arbitrary momentum up to maximum range
struct resampled 
integrate( const struct resampled *corr ,
	   const double max ,
	   const int L )
{
  struct resampled res ;
  res.resampled = calloc( corr[0].NSAMPLES , sizeof( double ) ) ;
  res.NSAMPLES  = corr[0].NSAMPLES ;
  res.restype   = corr[0].restype ;
  // compute twiddle factors
  size_t j ;
  for( j = 0 ; j < res.NSAMPLES ; j++ ) {
    Ltmp = L ;
    size_t t ;
    for( t = 0 ; t < Ltmp ; t++ ) {
      ct[t] = corr[t].resampled[j] ;
    }
    res.resampled[j] = adaptive_simpsons( evaluate_integrand_HVP ,
					  0.0 , 9 , 1E-9 ) ;
  }
  compute_err( &res ) ;
  mult_constant( &res , 4.0 / ( 137.03599878 * 137.03599878 ) ) ; 
  return res ;
}

// computes the t^2 moment
struct resampled
compute_PI0( struct resampled *corr ,
	     const size_t L )
{
  struct resampled temp , sum ;
  temp.resampled = malloc( corr[0].NSAMPLES * sizeof( double ) ) ;
  sum.resampled  = malloc( corr[0].NSAMPLES * sizeof( double ) ) ;
  equate_constant( &sum  , 0 , corr[0].NSAMPLES , corr[0].restype ) ;
  equate_constant( &temp , 0 , corr[0].NSAMPLES , corr[0].restype ) ;
  size_t j ;
  // add symmetrically
  for( j = 1 ; j < L/2 ; j++ ) {
    if( j != 0 ) {
      equate( &temp , corr[j] ) ;
      add( &temp , corr[L-j] ) ;
      mult_constant( &temp , 0.5*j*j ) ;
      add( &sum , temp ) ;
    }
  }
  // do the midpoint on its own to avoid double counting
  equate( &temp , corr[L/2] ) ;
  mult_constant( &temp , 0.5*j*j ) ;
  add( &sum , temp ) ;
  free( temp.resampled ) ;
  return sum ;
}

// time moment eval
void
arbitrary_mom_eval( double **xavg ,
		    struct resampled **bootavg ,
		    struct mom_info **mominfo ,
		    struct input_params *INPARAMS ,
		    const int NSLICES ,
		    const int LT ,
		    const bool renormalise )
{
  printf( "\n--> Arbitrary mom evaluation <--\n" ) ;

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

  // plot correlators
  make_xmgrace_graph( "corrs.agr" , "( aq )\\S2\\N" , "\\xP\\f{}" ) ;
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    plot_data( bootavg[i] , xavg[i] , INPARAMS->NDATA[i] ) ;
  }
  close_xmgrace_graph() ;

  // compute PI0
  PI0 = malloc( NSLICES/4 * sizeof( struct resampled ) ) ;
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    PI0[i] = compute_PI0( bootavg[i] , INPARAMS->dimensions[i][3] ) ;
    //divide( &PI0[i] ,  bootavg[i][INPARAMS->dimensions[i][3]/2] ) ;
    //printf( "%e %e \n" , pow( INPARAMS->dimensions[i][3] / INPARAMS -> quarks[0].ainverse , -4 ) , PI0[i].avg ) ;
    printf( "%e %e \n" , pow( INPARAMS->dimensions[i][3] * INPARAMS -> quarks[0].ainverse , 2 ) , PI0[i].avg ) ;
  }

  // test my silly ideas
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    int L ;
    for( L = 16 ; L <= INPARAMS->dimensions[i][3] ; L++ ) {
      struct resampled *PIS = compute_PIS( bootavg[i] , 12 , L , 
					   HPQCD_MOMENTS , 
					   INPARAMS->mom_type ) ;
      printf( "%e %e %e\n" , (double)L , PIS[0].avg , PIS[0].err ) ;
      /*
      struct resampled p0 = compute_PI0( bootavg[i] , L ) ;
      printf( "%e %e %e\n" , (double)L , p0.avg , p0.err ) ;
      free( p0.resampled ) ;
      */
    }
    printf( "\n" ) ;
  }

  printf( "\n--> Integrating amu <--\n" ) ;

  /*
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    struct resampled res = integrate( bootavg[i] , 
				      INPARAMS->dimensions[i][3]/2 -1 , 
				      INPARAMS->dimensions[i][3] ) ;
    printf( "[INTEGRATED] %e %e \n" , res.avg , res.err ) ;
    free( res.resampled ) ;
  }
  */

  struct resampled **data = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  double **x = malloc( NSLICES * sizeof( double* ) ) ;

  make_xmgrace_graph( "data.agr" , "( aq )\\S2\\N" , "\\xP\\f{}" ) ;
  
  for( i = 0 ; i < NSLICES/4 ; i++ ) {
    const double max = INPARAMS->dimensions[i][3]/2 -1 ;
    const int range = 300 ; // 100 measurements
    const double granularity = max / range ;
    x[i] = malloc( range * sizeof( double ) ) ;
    data[i] = malloc( range * sizeof( struct resampled ) ) ;
    int j ;
    #pragma omp parallel for private(j)
    for( j = 0 ; j < range ; j++ ) {
      double q = j * granularity , lmom ;
      data[i][j] = arbq( &lmom , bootavg[i] , q , 
				 INPARAMS->dimensions[i][3] ) ;
  
      add( &data[i][j] , PI0[i] ) ;
      x[i][j] = lmom * lmom ;
    }

    // plot him up
    plot_data( data[i] , x[i] , range ) ;
    plot_data( bootavg[2*NSLICES/4+i] , xavg[2*NSLICES/4+i] , INPARAMS->NDATA[2*NSLICES/4+i] ) ;
  }

  close_xmgrace_graph( ) ;
  
  return ;
}

