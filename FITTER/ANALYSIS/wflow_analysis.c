/**
   Small analysis code for the wilson flow
 */

#include "fitfunc.h"
#include "graph_data.h"

#if 0
void
analyse_w0( struct resampled **RAW , // need to resample this and root it ...
	    double *X ,              // have this as the am_l values
	    struct input_params INPARAMS ,// chiral data is in INPARAMS
	    const int NSLICES
	    )
{
  // binning
  struct resampled **BINNED = malloc( NSLICES * sizeof( struct resampled* ) ) ;

  int i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    BINNED[i] = malloc( INPARAMS.NDATA * sizeof( struct resampled ) ) ;
    int j ;
    for( j = 0 ; j < INPARAMS.NDATA ; j++ ) {
      BINNED[i][j] = bin_the_data( RAW[i][j] , INPARAMS.binning[j] ) ;
    }
  }

  // psq average or whatever ?
  struct resampled **BOOT = (struct resampled*)malloc( NSLICES * sizeof( struct resampled * ) ) ;

  // resample into BOOT
  for( i = 0 ; i < NSLICES ; i++ ) {
    BOOT[i] = resample_data( BINNED[i] , INPARAMS.NDATA , 
			     INPARAMS.resample , 
			     INPARAMS.NBOOTS ) ;
  }

  // take the square root
  double **sigma = malloc( NSLICES * sizeof( double* ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    sigma[i] = (double*)malloc( INPARAMS.NDATA * sizeof( double ) ) ;
    int j ;
    for( j = 0 ; j < INPARAMS.NDATA ; j++ ) {
      root( &BOOT[i][j] ) ;
      printf( "Result :: %f %f %f \n" , X[j] , BOOT[i][j].avg , BOOT[i][j].err ) ;
      sigma[i][j] = BOOT[i][j].err ;
    }
  }

  int *NPARAMS = malloc( NSLICES * sizeof( int ) ) ;

  // fit ?
  // and make sure the fit works
  INPARAMS.fit_lo = ( INPARAMS.fit_lo < 0 ) ? 0 : INPARAMS.fit_lo ;
  INPARAMS.fit_hi = ( INPARAMS.fit_hi > INPARAMS.NDATA - 1 ) ?	\
    INPARAMS.NDATA-1 : INPARAMS.fit_hi ;

  // have a nice little function for this
  fitfunc fit ;
  struct resampled **fitparams = perform_bootfit( &fit ,
						  NPARAMS ,
						  (const struct resampled**)BOOT , // the bootstrapped data 
						  X , // the (sorted) x-axis data
						  (const double**)sigma , // the s.d of the y-data
						  INPARAMS.quarks , // chiral data
						  INPARAMS.fittype , // the type of fit(s) we are using ...
						  NSLICES , // number of boots to fit
						  INPARAMS.NDATA , // the length of the X-array
						  INPARAMS.fit_hi ,
						  INPARAMS.fit_lo ) ;

  // and graph that badboy
  make_xmgrace_graph( INPARAMS.graph_name ) ;
  for( i = 0 ; i < NSLICES ; i++ ) { 
    plot_data( BOOT[i] , X , INPARAMS.NDATA ) ;
    plot_fitfunc( fitparams[i] , fit , NPARAMS[i] ,
		  INPARAMS.quarks[i] ,
		  0.0 , 
		  X[INPARAMS.fit_hi] ,
		  0 ) ;
  }
  close_xmgrace_graph(  ) ;
  printf( "Graph %s\n" , INPARAMS.graph_name) ;

  // free memory
  for( i = 0 ; i < NSLICES ; i++ ) {
    int j ;
    for( j = 0 ; j < INPARAMS.NDATA ; j++ ) {
      free( BOOT[i][j].resampled ) ;
      free( BINNED[i][j].resampled ) ;
    }
    free( BINNED[i] ) ;
    free( BOOT[i] ) ;
    // free the fit paramaters
    for( j = 0 ; j < NPARAMS[i] ; j++ ) {
      free( fitparams[i][j].resampled ) ;
    }
    free( fitparams[i] ) ;
  }
  free( BINNED ) ;
  free( BOOT ) ;
  free( fitparams ) ;
  free( NPARAMS ) ;
  free( sigma ) ;
 
  return ;
}

#endif
