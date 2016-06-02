/**
   @file dft.c
   @brief perform the time-like DFT
 */
#include "fitfunc.h"

#include <complex.h>
#include "GLUdata_glueprop.h"

// timelike DFT of the data
struct resampled **
momspace_corr( struct mom_info ***mom , 
	       double ***x ,
	       const struct resampled **tcorr ,
	       const struct input_params *INPARAMS ,
	       const int NSLICES , 
	       const int LT )
{
  *mom = malloc( NSLICES * sizeof( struct mom_info* ) ) ;
  struct resampled **PIP = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    (*mom)[ i ] = malloc( LT * sizeof( struct mom_info ) ) ;
    PIP[ i ] = malloc( LT * sizeof( struct resampled ) ) ;
    // perform the DFT
    const double twiddle = 2.0 * M_PI / (double)( INPARAMS -> dimensions[ i ][ 3 ] ) ; 
    int pt ;
    #pragma omp parallel for private( pt )
    for( pt = -LT/2 ; pt < LT/2 ; pt++ ) {
      const int p = pt + LT / 2 ;
      PIP[ i ][ p ].resampled = malloc( tcorr[ i ][ 0 ].NSAMPLES * sizeof( double ) ) ;
      equate_constant( &PIP[ i ][ p ] , 0.0 , tcorr[ i ][ 0 ].NSAMPLES ,
		       tcorr[ i ][ 0 ].restype ) ;
      double complex epimu[ LT ] ;
      size_t t ;
      for( t = 0 ; t < LT ; t++ ) {
	epimu[ t ] = cos( twiddle * pt * t ) ; //- I * sin( twiddle * pt * t ) ;
      }
      // loop boots
      size_t k ;
      for( k = 0 ; k < tcorr[ i ][ 0 ].NSAMPLES ; k++ ) {
	size_t t ;
	register double complex sum = 0.0 ;
	for( t = 0 ; t < LT ; t++ ) {
	  sum += epimu[ t ] * tcorr[ i ][ t ].resampled[ k ] ;
	}
	PIP[ i ][ p ].resampled[ k ] = creal( sum ) ;
      }
      // do the average
      register double complex sum = 0.0 ;
      for( t = 0 ; t < LT ; t++ ) {
	sum += epimu[ t ] * tcorr[ i ][ t ].avg ;
      }
      PIP[ i ][ p ].avg = creal( sum ) ;
 
      // divide by q^2
      int nmu[ 4 ] = { 0 , 0 , 0 , pt } ;
      (*mom)[ i ][ p ] = fill_mominfo( 4 , nmu , INPARAMS->dimensions[i] , 
				       INPARAMS->mom_type ) ;

      mult_constant( &PIP[ i ][ p ] , 
		     1.0 / ( (*mom)[ i ][ p ].p2 < 1E-15 ? 1.0 : (*mom)[ i ][ p ].p2 ) ) ;
    }
  }

  // allocate x and set it to the sqrt of p2
  for( i = 0 ; i < NSLICES ; i++ ) {
    (*x)[ i ] = malloc( INPARAMS->NDATA[i] * sizeof( double ) ) ;
    size_t j ;
    for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
      (*x)[i][j] = sqrt( (*mom)[i][j].p2 ) ;
    }
  }
  return PIP ;
}
