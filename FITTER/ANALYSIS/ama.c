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

//#define ONLY_CORRECTED

// ama computation
void
ama_eval( double **xavg ,
	  struct resampled **bootavg ,
	  struct mom_info **mominfo ,
	  struct mom_info *moms ,
	  struct input_params *INPARAMS ,
	  const int NSLICES ,
	  const int LT )
{
  // Sloppy on the same place should be first, then exact
  // then the symmetry-shifted
  const int NCHANNELS = INPARAMS -> NCHANNELS ;

  const int tru = (NSLICES/NCHANNELS)/3 ;

#ifdef ONLY_CORRECTED
  struct resampled **corrected = malloc( tru * NCHANNELS * sizeof( struct resampled* ) ) ;
  int i , j , k ;
  for( k = 0 ; k < tru ; k++ ) {
    for( i = 0 ; i < NCHANNELS ; i++ ) {
      const int idx = i + 3*NCHANNELS*k ;
      const int cidx = i+k*NCHANNELS ;
      corrected[cidx] = malloc( INPARAMS -> NDATA[idx] * sizeof( struct resampled ) ) ;
      for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
	corrected[cidx][j].resampled = malloc( bootavg[idx+NCHANNELS ][j].NSAMPLES* sizeof( double ) ) ;
	equate( &corrected[cidx][j] , bootavg[idx][j] ) ;
	// computes O^{rest} in bootavg[0]
	subtract( &corrected[cidx][j] , bootavg[idx+NCHANNELS][j] ) ;
	mult_constant( &corrected[cidx][j] , -1 ) ;
	// adds the bias correction to the sloppys
	add( &corrected[cidx][j] , bootavg[idx+2*NCHANNELS][j] ) ;
      }
    }
  }
  correlator_eval( xavg , corrected , mominfo , moms , 
		   INPARAMS , tru*NCHANNELS , LT ) ;
  free_resampled_dp( corrected , tru*NCHANNELS , INPARAMS->NDATA ) ;

#else
  int i , j , k ;
  for( k = 0 ; k < tru ; k++ ) {
    for( i = 0 ; i < NCHANNELS ; i++ ) {
      const int idx = i + 3*NCHANNELS*k ;
      for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
	// computes O^{rest} in bootavg[0]
	subtract( &bootavg[idx][j] , bootavg[idx+NCHANNELS][j] ) ;
	mult_constant( &bootavg[idx][j] , -1 ) ;

	//#if 0
	divide( &bootavg[idx][j] , bootavg[idx+NCHANNELS][j] ) ;
	mult_constant( &bootavg[idx][j] , 100 ) ;
	//#endif
	// adds the bias correction to the sloppys
	//add( &bootavg[idx][j] , bootavg[idx+2*NCHANNELS][j] ) ; 
      }
    }
  }
  // evaluate the correlator now
  correlator_eval( xavg , bootavg , mominfo , moms , 
		   INPARAMS , NSLICES , LT ) ;
#endif
  return ;
}
