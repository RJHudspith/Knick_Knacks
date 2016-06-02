#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "dispersions.h"
#include "Utils.h"
#include "write_distribution.h"
#include "read_distribution.h"
#include "effmass.h"
#include "GLU_bswap.h"
#include "correlators.h"

#define COMPUTE_BINDING

#ifdef COMPUTE_BINDING

static void
flip_data( struct resampled *flipped ,
	   const struct resampled *unflipped ,
	   const int NDATA )
{
  size_t j ;
  equate( &flipped[0] , unflipped[0] ) ;
  for( j = 1 ; j < NDATA ; j++ ) {
    equate( &flipped[j] , unflipped[NDATA-j] ) ;
  }
}

static void
average_data( struct resampled *ave ,
	      const struct resampled *to_add ,
	      const int NDATA )
{
  size_t j ;
  for( j = 0 ; j < NDATA ; j++ ) {
    add( &ave[j] , to_add[j] ) ;
    mult_constant( &ave[j] , 0.5 ) ;
  }
}

#endif

// tetra computation
void
tetra_eval( double **xavg ,
	    struct resampled **bootavg ,
	    struct mom_info **mominfo ,
	    struct mom_info *moms ,
	    struct input_params *INPARAMS ,
	    const int NSLICES ,
	    const int LT )
{
  // Diquark-Diquark is in Slice #0 and Dimeson is in #1
  // Vector Meson is in #2
  // Pseudoscalar Meson is in #3
  size_t j ;
#ifdef COMPUTE_BINDING
  if( NSLICES == 4 ) {
    // take product of vector and pseudoscalar put into 2
    for( j = 0 ; j < INPARAMS -> NDATA[2] ; j++ ) {
      mult( &bootavg[2][j] , bootavg[3][j] ) ;
      mult_constant( &bootavg[2][j] , -1 ) ;
    }
  } else if( NSLICES == 8 ) {
    // take product of vector and pseudoscalar
    for( j = 0 ; j < INPARAMS -> NDATA[2] ; j++ ) {
      mult( &bootavg[2][j] , bootavg[3][j] ) ;
      mult_constant( &bootavg[2][j] , -1 ) ;
      
      mult( &bootavg[6][j] , bootavg[7][j] ) ;
      mult_constant( &bootavg[6][j] , -1 ) ;
    }

    // backwards data needs to b time flipped
    flip_data( bootavg[3] , bootavg[4] , INPARAMS -> NDATA[4] ) ;
    flip_data( bootavg[4] , bootavg[5] , INPARAMS -> NDATA[4] ) ;
    flip_data( bootavg[5] , bootavg[6] , INPARAMS -> NDATA[4] ) ;

    // average the data
    average_data( bootavg[0] , bootavg[3] , INPARAMS -> NDATA[0] ) ;
    average_data( bootavg[1] , bootavg[4] , INPARAMS -> NDATA[0] ) ;
    average_data( bootavg[2] , bootavg[5] , INPARAMS -> NDATA[0] ) ;
  }

  for( j = 0 ; j < INPARAMS -> NDATA[2] ; j++ ) {
    // divide correlator #0 by this result
    divide( &bootavg[0][j] , bootavg[2][j] ) ;
    divide( &bootavg[1][j] , bootavg[2][j] ) ;
  }
  
  // evaluate the correlator now
  correlator_eval( xavg , bootavg , mominfo , moms , 
		   INPARAMS , 2 , LT ) ;
#else
  // take product of vector and pseudoscalar put into 2
  for( j = 0 ; j < INPARAMS -> NDATA[2] ; j++ ) {
    mult( &bootavg[2][j] , bootavg[3][j] ) ;
    mult_constant( &bootavg[2][j] , -1 ) ;
  }
  // evaluate the correlator now
  correlator_eval( xavg , bootavg , mominfo , moms , 
		   INPARAMS , 3 , LT ) ;
#endif

  return ;
}
