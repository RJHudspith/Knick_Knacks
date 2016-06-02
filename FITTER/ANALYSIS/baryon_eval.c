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


// tetra computation
void
baryon_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct mom_info *moms ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT )
{
  // if we have NSLICES == 4, average L4 and L5
  if( NSLICES == 2 ) {
    size_t j ;
    for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
      add( &bootavg[0][j] , bootavg[1][j] ) ;
      mult_constant( &bootavg[0][j] , 0.5 ) ;
    }
    // evaluate the correlator now
    correlator_eval( xavg , bootavg , mominfo , moms , 
		     INPARAMS , 1 , LT ) ;
  } else if( NSLICES == 4 ) {
    size_t j ;
    for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
      add( &bootavg[0][j] , bootavg[1][j] ) ;
      // do the walls
      equate( &bootavg[1][j] , bootavg[2][j] ) ;
      add( &bootavg[1][j] , bootavg[3][j] ) ;
      // divide by 2
      mult_constant( &bootavg[0][j] , 0.5 ) ;
      mult_constant( &bootavg[1][j] , 0.5 ) ;
    }
    // evaluate the correlator now
    correlator_eval( xavg , bootavg , mominfo , moms , 
		     INPARAMS , 2 , LT ) ;
  } else {
    // evaluate the correlator now
    correlator_eval( xavg , bootavg , mominfo , moms , 
		     INPARAMS , NSLICES , LT ) ;
  }
  return ;
}
