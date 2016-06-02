/**
   given a bunch of distributions we compute their correlation
 */
#include "fitfunc.h"
#include "correlation.h"

void
correlation( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT )
{
  printf( "CHECK :: %d %d \n", INPARAMS->NDATA[0] , NSLICES ) ;
  // correlation matrix
  double **correlation = malloc( INPARAMS->NDATA[0] * sizeof( double ) ) ;

  size_t i ;
  for( i = 0 ; i < INPARAMS->NDATA[0] ; i++ ) {
    correlation[i] = malloc( NSLICES * sizeof( double ) ) ;
  }
  
  // compute the correlation matrix
  correlations( correlation , bootavg[0] , INPARAMS->NDATA[0] ) ;

  // write out a mathematica-friendly file
  write_corrmatrix_mathematica( correlation , INPARAMS->NDATA[0] ) ;

  return ;
}
