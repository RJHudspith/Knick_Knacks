/**
   individual function for generic fitting
 */

#include "fitfunc.h"
#include "fitdata.h"
#include "Utils.h"
#include "graph_data.h"
#include "Utils.h"
#include "correlation.h"
#include "derivatives.h"

#define EX1

// correlation matrix
struct cmatrix *C = NULL ;

// fit and perhaps plot, returning the fit results
struct resampled *
fit_data_plot_data( const struct resampled **BOOT ,
		    const double **X ,
		    const struct mom_info **mominfo ,
		    const struct input_params INPARAMS ,
		    const int NSLICES ,
		    const int LT )
{
  int i ;

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS.graph_name ,
		      INPARAMS.graph_xaxis , 
		      INPARAMS.graph_yaxis ) ;

  // plot the data first but keep the graph file open
  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_data( BOOT[i] , X[i] , INPARAMS.NDATA[i] ) ;
  }
  graph_reset_color( ) ;

  struct resampled *fitparams = NULL ;

  // perform the fit if it is required
  if( INPARAMS.fittype != NOFIT ) {

    int sum = 0 ;
    for( i = 0 ; i < NSLICES ; i++ ) {
      sum += INPARAMS.NDATA[i] ;
    }
    const int NFLAT[1] = { sum } ;

    // number of common fit parameters
    int NCOMMON = 0 ;
    for( i = 0 ; i < 12 ; i++ ) {
      if( INPARAMS.sim_params[ i ] == true ) NCOMMON++ ;
    }

    // compute a flattend x-axis index ...
    double *flatx = flatten_double_array( X , NSLICES , INPARAMS.NDATA ) ;

    // flatten the mom info
    struct mom_info *flatmom = flatten_mom_array( mominfo , NSLICES , INPARAMS.NDATA ) ;

    // flatten BOOTS into a single array ...
    struct resampled *RES = flatten_resampled_array( BOOT , NSLICES , INPARAMS.NDATA ) ;

    // compute a flattened sigma
    double *sigma = resampled_to_double( (const struct resampled*)RES , NFLAT[0] , ERR ) ;

    #ifdef VERBOSE
    // print out flattened data and sigmas
    {
      int j ;
      printf( "{" ) ;
      for( j = 0 ; j < NFLAT[0] ; j++ ) {
	printf( "{ %e , %e , %e } ,\n" , flatx[j] , RES[j].avg , RES[j].err ) ;
      }
      printf( "}\n" ) ;

      printf( "{" ) ;
      for( j = 0 ; j < NFLAT[0] ; j++ ) {
	printf( " %1.12e ,\n" , 1.0 / ( sigma[j] * sigma[j] ) ) ;
      }
      printf( "}\n" ) ;
    }
    #endif

    // perform the fit?
    int *NPARAMS = malloc( NSLICES * sizeof( int ) ) ;
    fitfunc fit ;

    // tell us the fit ranges
    printf( "\n" ) ;
    int j ;
    for( j = 0 ; j < NSLICES ; j++ ) {
      printf( "FIT RANGE %d :: %d %d \n" , j , 
	      find_idx( INPARAMS.fit_lo , X[j] , 
			INPARAMS.NDATA[j] , 0 ) ,
	      find_idx( INPARAMS.fit_hi , X[j] , 
			INPARAMS.NDATA[j] , 0 ) ) ;
    }
    printf( "\n" ) ;

    fitparams = perform_bootfit( &fit ,
				 NPARAMS ,
				 (const struct resampled*)RES , // the bootstrapped data 
				 (const double*)flatx , // the (sorted) x-axis data
				 (const double*)sigma , // the s.d of the y-data
				 (const struct mom_info*)flatmom , // flattened momenta
				 INPARAMS.quarks , // chiral data
				 INPARAMS.fittype , // the type of fit(s) we are using ...
				 INPARAMS.NDATA , // the length of the X-array
				 INPARAMS.fit_hi ,
				 INPARAMS.fit_lo ,
				 NSLICES ,
				 LT ,
				 INPARAMS.sim_params ,
				 NCOMMON ) ;

    // plot and stuff
    // unpack the simultaneous fit into something more palatable
    struct resampled *fit1 = malloc( NPARAMS[0] * sizeof( struct resampled ) ) ;

    for( i = 0 ; i < NPARAMS[0] ; i++ ) {
      fit1[i].resampled = (double*)malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
    }

    // loop slices plotting fitfunctions
    for( i = 0 ; i < NSLICES ; i++ ) {

      int k ;
      if( i > 0 ) {
	int check = i*( NPARAMS[0] - NCOMMON ) + NCOMMON ;
	for( k = 0 ; k < NPARAMS[0] ; k++ ) {
	  if( INPARAMS.sim_params[ k ] == true ) {
	    equate( &fit1[ k ] , fitparams[ k ] ) ;
	  } else {
	    equate( &fit1[ k ] , fitparams[ check ] ) ;
	    check++ ;
	  }
	}
      } else {
        #pragma omp parallel for private(k)
	for( k = 0 ; k < NPARAMS[0] ; k++ ) {
	  equate( &fit1[ k ] , fitparams[ k ] ) ;
	}
      }
      switch( INPARAMS.fittype ) {
      case PPAA_WW : case PPAA : case VV_WW : case VV_VTVT_WW :
	plot_fitfunc2( fit1 , BOOT[i] , X[i] , mominfo[i] , fit , NPARAMS[0] ,
		       INPARAMS.quarks[i] , INPARAMS.fit_lo , INPARAMS.fit_hi ,
		       LT , INPARAMS.NDATA[i] , i ) ;
	break ;
      default :
	plot_fitfunc2( fit1 , BOOT[i] , X[i] , mominfo[i] , fit , NPARAMS[0] ,
		       INPARAMS.quarks[i] , INPARAMS.fit_lo , INPARAMS.fit_hi ,
		       LT , INPARAMS.NDATA[i] , NPARAMS[i] ) ;
	break ;
      }
    }
    // free the number of fit params array
    free( NPARAMS ) ;

    // free the flattened arrays
    free( flatx ) ;
    free( sigma ) ;
    free_resampled( RES , NFLAT[0] ) ;
    free( flatmom ) ;
  }

  // close up the graph
  close_xmgrace_graph(  ) ;

  printf( "\n---> Graph plotted to %s <---\n" , INPARAMS.graph_name) ;

  return fitparams ;
}
