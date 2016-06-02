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

#define LINEAR
//#define SUBTRACTED

// dispersions computation
void
dispersions_eval( double **xavg ,
		  struct resampled **bootavg ,
		  struct mom_info **mominfo ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT )
{
  struct resampled mass[ NSLICES ] ;

  int j ;
  for( j = 0 ; j < NSLICES ; j++ ) {
    mass[j].resampled = malloc( bootavg[j][0].NSAMPLES * sizeof( double ) ) ;
    equate( &mass[j] , bootavg[j][0] ) ;

    #ifdef SUBTRACTED
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      subtract( &bootavg[j][i] , mass[j] ) ;
      divide_constant( &bootavg[j][i] , xavg[j][i] > 0.0 ? xavg[j][i] : 1.0 ) ;
    }
    #elif defined LINEAR
    // subtract the mass
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      subtract( &bootavg[j][i] , mass[j] ) ;
    }
    #else
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      // square both sides
      mult( &bootavg[j][i] , bootavg[j][i] ) ;
      divide( &bootavg[j][i] , mass[j] ) ;
      subtract( &bootavg[j][i] , mass[j] ) ;
      divide_constant( &bootavg[j][i] , xavg[j][i] > 0.0 ? xavg[j][i] : 1.0 ) ;
    }
    #endif
  }

#if ( defined SUBTRACTED ) || ( defined LINEAR )
  for( j = 0 ; j < NSLICES ; j++ ) {
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      mult_constant( &bootavg[j][i] , 2 ) ;
    }
  }
#endif

  for( j = 0 ; j < NSLICES ; j++ ) {
    printf( "NDATA :: %d \n" , INPARAMS -> NDATA[j] ) ;
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      printf( "[%d,%d] %f %f \n" , j , i , xavg[j][i] , bootavg[j][i].avg ) ;
    }
  }

  // plot only the flavour breaking difference
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NSLICES , LT ) ;
  struct resampled kinetic_mass ;
  kinetic_mass.resampled = malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
#ifdef LINEAR
  equate( &kinetic_mass , fitparams[1] ) ;
#else
  equate( &kinetic_mass , fitparams[0] ) ;
#endif
  raise( &kinetic_mass , -1 ) ; 

  printf( "[DISPERSIONS] ZERO_MOM_MASS %1.12f %1.12f \n" , mass[0].avg , mass[0].err ) ;

  printf( "[DISPERSIONS] KINETIC_MASS %1.12f %1.12f \n" , kinetic_mass.avg , kinetic_mass.err ) ;

  free( kinetic_mass.resampled ) ;

  for( j = 0 ; j < NSLICES ; j++ ) {
    free( mass[j].resampled ) ;
  }

  return ;
}

