/**
   Analysis for the VPF
 */
#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "dispersions.h"
#include "Utils.h"
#include "write_superdist.h"

// dispersions computation
void
speed_of_light_eval( double **xavg ,
		   struct resampled **bootavg ,
		   struct mom_info **mominfo ,
		   struct input_params *INPARAMS ,
		   const int NSLICES ,
		   const int LT )
{
  int j ;
  for( j = 0 ; j < NSLICES ; j++ ) {

    // compute the ground state mass, sorted so always bootavg[j][0]
    struct resampled mass ;
    mass.resampled = malloc( bootavg[j][0].NSAMPLES * sizeof( double ) ) ;

    // compute ( m.c(0)^2 ) ^2
    equate( &mass , bootavg[j][0] ) ;
    mult( &mass , mass ) ;

    // compute c^2 = ( E(p^2)^2 - ( M_0.c^2 )^2 ) / ( p^2 )
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      mult( &bootavg[j][i] , bootavg[j][i] ) ;
      subtract( &bootavg[j][i] , mass ) ;
      divide_constant( &bootavg[j][i] , xavg[j][i] > 0.0 ? xavg[j][i] : 1.0 ) ;
    } 

    free( mass.resampled ) ;
  }

  // plot only the flavour breaking difference
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NSLICES , LT ) ;

  printf( "[SOL] SPEED %1.12f %1.12f \n" , fitparams[0].avg , fitparams[0].err ) ;

  free( fitparams ) ;

  return ;
}

