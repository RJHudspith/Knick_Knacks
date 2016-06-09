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

//#define COMPUTE_BINDING

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
tetra_eval_meson( double **xavg ,
		  struct resampled **bootavg ,
		  struct mom_info **mominfo ,
		  struct mom_info *moms ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT )
{
  size_t j ;
  for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
    mult( &bootavg[0][j] , bootavg[1][j] ) ;
    //equate( &bootavg[1][j] , bootavg[2][j] ) ;
    //mult( &bootavg[1][j] , bootavg[3][j] ) ;
  }

  // evaluate the correlator now
  correlator_eval( xavg , bootavg , mominfo , moms , 
		   INPARAMS , 1 , LT ) ;

  return ;
}

// tetra computation
void
tetra_eval_corr( double **xavg ,
		 struct resampled **bootavg ,
		 struct mom_info **mominfo ,
		 struct mom_info *moms ,
		 struct input_params *INPARAMS ,
		 const int NSLICES ,
		 const int LT )
{
  printf( "TETRA CORR EVAL" ) ;
  // evaluate the correlator now
  correlator_eval( xavg , bootavg , mominfo , moms , 
		   INPARAMS , NSLICES , LT ) ;
  return ;
}
