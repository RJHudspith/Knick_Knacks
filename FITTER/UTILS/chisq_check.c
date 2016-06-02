/**
   @file chisq_check.c
   @brief compute the chisq
 */
#include "fitfunc.h"

#include "correlation.h"

// for checking the chisq
struct resampled
chisq_check( const struct resampled *fitparams ,
	     const struct resampled *ydata , 
	     const double *xdata , 
	     const int NPARAMS ,
	     const struct chiral *quarks ,
	     const struct mom_info *mominfo ,
	     const fitfunc fit ,
	     const int Ndata ,
	     const int LT )
{
  struct resampled chisq ;
  chisq.resampled = malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
  chisq.NSAMPLES  = fitparams[0].NSAMPLES ;
  chisq.restype   = fitparams[0].restype ;

  size_t j ;
  // loop boots
  #pragma omp parallel for private(j)
  for( j = 0 ; j < fitparams[0].NSAMPLES ; j++ ) {
    double data[ NPARAMS ] ;
    // switch to individual axis
    size_t k ;
    for( k = 0 ; k < NPARAMS ; k++ ) {
      data[ k ] = fitparams[k].resampled[j] ;
    }
    double eval[ Ndata ] ;
    for( k = 0 ; k < Ndata ; k++ ) {
      struct x_descriptor XX = { xdata[k] , quarks[k] , mominfo[k] , LT } ;
      eval[ k ] = ( ydata[k].resampled[j] - fit.f( data , XX , NPARAMS ) ) / ydata[k].err ;
    }
    // loop elements
    double chi = 0.0 ;
    size_t i ;
    for( i = 0 ; i < Ndata ; i++ ) {
      chi += eval[i] * eval[i] ;
    }
    chisq.resampled[j] = chi / ( Ndata - NPARAMS ) ;
  }
  compute_err( &chisq ) ;

  printf( "CHISQ :: %f +/- %f \n" , chisq.avg , chisq.err ) ;

  return chisq ;
}
