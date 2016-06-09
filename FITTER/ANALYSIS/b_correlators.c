/**
   @file b_correlators.c
   @brief average the forward and backward correlators
 */
#include "fitfunc.h"

#include "correlators.h"

void
average_pp( struct resampled *fwd ,
	    struct resampled *bwd ,
	    const size_t i ,
	    const size_t j1 ,
	    const size_t j2 )
{
  add( &fwd[j1] , bwd[j2] ) ;
  mult_constant( &fwd[j1] , 0.5 ) ; 
}

// tetra computation
void
bmeson_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct mom_info *moms ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT )
{
  // average forward and time-reversed backwards props
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < NSLICES/2 ; i++ ) {
    size_t j ;
    average_pp( bootavg[i] , bootavg[i+NSLICES/2] , 
		i , 0 , 0 ) ;

    for( j = 1 ; j < INPARAMS -> NDATA[i] ; j++ ) {
      average_pp( bootavg[i] , bootavg[i+NSLICES/2] , 
		  i , j , INPARAMS -> NDATA[0] - j ) ;
    }
  }

  // evaluate the correlator now
  correlator_eval( xavg , bootavg , mominfo , moms , 
		   INPARAMS , NSLICES/2 , LT ) ;

  return ;
}
