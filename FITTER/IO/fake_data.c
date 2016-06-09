#include "fitfunc.h"
#include "rng.h"
#include "fit_chooser.h"
#include "D0_fit.h"

//#define FAKE_VPF

////////////////////////////////////////////////////////////////////////////
// fake a bunch of exponentials with common mass parameter
struct resampled **
fake_data( double ***X , /// ewww! I feel pretty bad doing this
	   const double XRANGE ,
	   const fit_type fittype ,
	   const struct chiral *quarks ,
	   const int NSLICES ,
	   int *NDATA ,
	   const int NRAW ,
	   const int LT ,
	   const bool sim_params[12] )
{
  // make each NDATA a different length
  int slice ;
  for( slice = 0 ; slice < NSLICES ; slice++ ) {
    NDATA[slice] = rng_int( 50 ) + 50 ; 
    printf( "\nNSLICES :: %d NDATA :: %d NRAW :: %d \n" , NSLICES , NDATA[slice] , NRAW ) ;
  }
  printf( "\n" ) ;

  struct mom_info dummy_mom ;

  fitfunc fit ; // this is set, and allows us to call the eval
  int NPARAMS , NCOMMON = 0 ;
  initialise_f( fittype , &fit , &NPARAMS , &NCOMMON ) ;

  // set up the weights, weight 1 is for the common parameter, weight2
  // is for the simultaneous parameters ...
  double weight1 = 0.1 , weight2 = 0.5 ;

  // allocate the X data
  *X = (double**)malloc( NSLICES * sizeof( double* ) ) ;

  // allocate the RAW data
  struct resampled **RAW = malloc( NSLICES * sizeof( struct resampled* ) ) ;

  // init fit params
  double params[ NPARAMS ] , sim_weight[ NPARAMS ] ;

  // fix the simultaneous params first
  int NSIM = 0 , k ;
  for( k = 0 ; k < NPARAMS ; k++ ) {
    if( sim_params[k] == true ) {
      sim_weight[ NSIM ] = weight1 ;
      NSIM++ ;
    }
  }

  // loop number of datasets
  for( slice = 0 ; slice < NSLICES ; slice++ ) {

    // allocations
    (*X)[slice] = (double*)malloc( NDATA[slice] * sizeof( double ) ) ;
    RAW[slice] = malloc( NDATA[slice] * sizeof( struct resampled ) ) ;

    // set the others to pseudo random values
    NSIM = 0 ;
    for( k = 0 ; k < NPARAMS ; k++ ) {
      if( sim_params[k] == true ) {
	params[ k ] = sim_weight[ NSIM ] ;
	NSIM++ ;
      } else {
	params[ k ] = ( weight2 + weight2 ) * rng_double()  ;
      }
    }

    // provide a description of the function 
    struct x_descriptor temp = { 0.0 , quarks[slice] , dummy_mom , LT } ;
    #ifdef FAKE_VPF
    const double nparams[2] = { 0.095 , -0.18 } ;
    D0_description( "FAKED" , nparams , temp , NPARAMS ) ;
    #else
    fit.description( "FAKED" , params , temp , NPARAMS ) ;
    #endif
    //exit(1) ;
    int i ;
    for( i = 0 ; i < NDATA[slice] ; i++ ) {

      // set the x values
      #ifdef FAKE_VPF
      (*X)[ slice ][ i ] = ( rng_double() ) * XRANGE * quarks[slice].ainverse ;
      #else
      (*X)[ slice ][ i ] = ( rng_double() ) * XRANGE ;
      #endif

      // allocate the resampled data ...
      RAW[slice][i].resampled = malloc( NRAW * sizeof( double ) ) ;
      RAW[slice][i].NSAMPLES  = NRAW ;
      RAW[slice][i].restype   = RAWDATA ;

      // and make a fake data ...

      struct x_descriptor XX = { (*X)[slice][i] , quarks[slice] , dummy_mom , LT } ;
      #ifdef FAKE_VPF
      const double data = D0_eval( nparams , XX , NPARAMS ) ;
      #else
      const double data = fit.f( params , XX , NPARAMS ) ;
      #endif

      int j ;
      for( j = 0 ; j < NRAW ; j++ ) {
	RAW[slice][i].resampled[j] = ( data * ( rng_gaussian( 0.005 ) + 1 ) ) ; 
      }
    }
  }
  printf( "\n" ) ;

  return RAW ;
}

