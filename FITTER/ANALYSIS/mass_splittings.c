/**
   @file mass_splittings.c
   @brief computes ratios of masses
 */
#include "fitfunc.h"

void
mass_splittings_eval( double **xavg ,
		      struct resampled **bootavg ,
		      struct mom_info **mominfo ,
		      struct input_params *INPARAMS ,
		      const int NSLICES ,
		      const int LT )
{
  struct resampled a , result ;
  a.resampled = malloc( bootavg[0][0].NSAMPLES * sizeof( double ) ) ;
  equate( &a , bootavg[0][0] ) ;
  result.resampled = malloc( bootavg[0][0].NSAMPLES * sizeof( double ) ) ;

  printf( "%d %d \n" , NSLICES , INPARAMS->NDATA[0] ) ;

  size_t i ;
  for( i = 0 ; i < INPARAMS->NDATA[0] ; i++ ) {
    printf( "[SPLITTINGS] mass_%zu :: %e +/- %e \n" , 
	    i , bootavg[0][i].avg , bootavg[0][i].err ) ;
  }

  // possible combinations
  switch( INPARAMS -> NDATA[0] ) {
  case 2 : // computes ( a - b ) / a
    equate( &result , bootavg[0][0] ) ;
    subtract( &result , bootavg[0][1] ) ;
    printf( "[SPLITTING] ( a - b ) :: %e +/- %e \n" , 
	    result.avg , result.err ) ;
    divide( &result , a ) ;
    printf( "[SPLITTING] ( a - b ) / a :: %e +/- %e \n" , 
	    result.avg , result.err ) ;
    // computes a^2
    equate( &result , bootavg[0][0] ) ;
    raise( &result , 2 ) ;
    printf( "[POWER] ( a )^2 :: %e +/- %e \n" , 
	    result.avg , result.err ) ;
    // computes b^2
    equate( &result , bootavg[0][1] ) ;
    raise( &result , 2 ) ;
    printf( "[POWER] ( b )^2 :: %e +/- %e \n" , 
	    result.avg , result.err ) ;
    // computes a^2 / b^2
    equate( &result , bootavg[0][0] ) ;
    divide( &result , bootavg[0][1] ) ;
    raise( &result , 2 ) ;
    printf( "[RATIO] ( a / b )^2 :: %e +/- %e \n" , 
	    result.avg , result.err ) ;
    break ;
  default :
    break ;
  }
  free( a.resampled ) ;
  free( result.resampled ) ;
  return ;
}
