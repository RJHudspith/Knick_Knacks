/**
   @file statpot.c
   @brief static potential analysis
 */
#include "fitfunc.h"
#include "fit_and_plot.h"

void
statpot_eval( double **xavg ,
	      struct resampled **bootavg ,
	      struct mom_info **mominfo ,
	      const struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT )
{
  int N = 0 ;
  switch( INPARAMS -> fittype ) {
  default : break ;
  }

  printf( "Static potential evaluation .... \n" ) ;
  // multiply by r?
  size_t i , j ;
  for( j = 0 ; j < NSLICES ; j++ ) {
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      xavg[j][i] = sqrt( xavg[j][i] ) ;

      mult_constant( &bootavg[j][i] , pow( xavg[j][i] , N ) ) ;
    }
  }
  // fit and plot are in here
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info**)mominfo , 
						    *INPARAMS , NSLICES , LT ) ;


  // compute r_0
  if( INPARAMS -> fittype == POLY2 ) {
    struct resampled res ; 
    res.resampled = malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
    equate( &res , fitparams[0] ) ;
    add_constant( &res , 1.65 ) ;
    divide( &res , fitparams[2] ) ;
    root( &res ) ;
    printf( "r_0 :: %f %f \n" , res.avg , res.err ) ;
  }
  // general Newton-Raphson or something?


#ifdef PRINT_DATA
  printf( "(r/a)^2   V(r)   Err  \n" ) ;
  for( j = 0 ; j < NSLICES ; j++ ) {
    printf( "T = %d \n" , 2 + j ) ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      printf( "%f %f %f \n" , xavg[j][i] , bootavg[j][i].avg , bootavg[j][i].err ) ;
    }
    printf( "\n" ) ;
  }
#endif
  return ;
}




