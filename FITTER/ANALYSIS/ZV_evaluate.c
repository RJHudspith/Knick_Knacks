
#include "fitfunc.h"
#include "fit_and_plot.h"

// expects the local current to be in NSLICE1
void
ZV_eval( double **xavg ,
	 struct resampled **bootavg ,
	 struct mom_info **mominfo ,
	 struct mom_info *moms ,
	 struct input_params *INPARAMS ,
	 const int NSLICES ,
	 const int LT )
{
  if( NSLICES != 2 ) {
    printf( "Expected NSLICES = 2\n" ) ;
  }

  // divide LL by CL
  size_t i , j ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    for( j = 0 ; j < INPARAMS -> NDATA[i] ; j++ ) {
      divide( &bootavg[0][j] , bootavg[1][j] ) ;
    }
  }

  // plot the correlator(s)
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , 1 , LT ) ;

  free( fitparams ) ;

  return ;
}
