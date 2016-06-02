/**
   @file mass_splittings.c
   @brief computes ratios of masses
 */
#include "fitfunc.h"
#include "fit_and_plot.h"
#include "resampled_ops.h"

#include "svd.h"
#include "write_distribution.h"
#include "graph_data.h"
#include "fit_chooser.h"

static const double MPISQ = 0.01821868254756 ;

void
mass_extrap2( double **xavg ,
	      struct resampled **bootavg ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT )
{
  int NPARAMS ;
  fitfunc fit ;
  initialise_f( INPARAMS->fittype , &fit , &NPARAMS , 0 ) ;

  printf( "NPARAMS :: %d\n" , NPARAMS ) ;
  /*
  if( INPARAMS->fittype != POLY1 ) {
    printf( "[SS-EXTRAP] expected to be linear!\n" ) ;
    return ;
  }
  */

  size_t i ;
  for( i = 0 ; i < INPARAMS->NDATA[0] ; i++ ) {
    xavg[0][i] = ( INPARAMS -> quarks[i].ml - MPISQ ) ;
    printf( "[MASSextrap] %e :: %e +/- %e \n" , 
	    xavg[0][i] , bootavg[0][i].avg , bootavg[0][i].err ) ;
  }

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS -> graph_name ,
		      INPARAMS -> graph_xaxis , 
		      INPARAMS -> graph_yaxis ) ;

  // plot the data first but keep the graph file open
  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_data( bootavg[i] , xavg[i] , INPARAMS -> NDATA[i] ) ;
  }

  // allocate fit params
  struct resampled *fparams = malloc( NPARAMS * sizeof( struct resampled ) ) ;

  // linearly solve using the svd
  {
    size_t j ;

    // loop NSLICES
    for( i = 0 ; i < 1 ; i++ ) {

      const int range = INPARAMS -> NDATA[ i ] ;

      // set fparams to the same as bootparams
      for( j = 0 ; j < NPARAMS ; j++ ) {
	fparams[j].resampled = malloc( bootavg[i][0].NSAMPLES * sizeof( double ) ) ;
	fparams[j].restype = bootavg[i][0].restype ;
	fparams[j].NSAMPLES = bootavg[i][0].NSAMPLES ;
      }

      // set the average first
      {
	// set the data
	double ydata[ range ] , xs[ range ] , sigma[ range ] , coeffs[ NPARAMS ] , chisq ;
	int k ;
	for( k = 0 ; k < range ; k++ ) {
	  ydata[ k ] = bootavg[ i ][ k ].avg ;
	  sigma[ k ] = bootavg[ i ][ k ].err ;
	  xs[ k ] = xavg[ i ][ k ] ;
	}
	//
	compute_coefficients( coeffs , &chisq , ydata , sigma , 
			      xs , range , NPARAMS ) ;
	for( k = 0 ; k < NPARAMS ; k++ ) {
	  fparams[ k ].avg = coeffs[ k ] ;
	}
      }

      // loop boots
#pragma omp parallel for private(j)
      for( j = 0 ; j < bootavg[i][0].NSAMPLES ; j++ ) {

	// set the data
	double ydata[ range ] , xs[ range ] , sigma[ range ] , coeffs[ NPARAMS ] , chisq ;

	int k ;
	for( k = 0 ; k < range ; k++ ) {
	  ydata[ k ] = bootavg[ i ][ k ].resampled[ j ] ;
	  sigma[ k ] = bootavg[ i ][ k ].err ;
	  xs[ k ] = xavg[ i ][ k ] ;
	}
	//
	compute_coefficients( coeffs , &chisq , ydata , sigma , 
			      xs , range , NPARAMS ) ;

	for( k = 0 ; k < NPARAMS ; k++ ) {
	  fparams[ k ].resampled[ j ] = coeffs[ k ] ;
	}
      }
      for( j = 0 ; j < NPARAMS ; j++ ) {
	compute_err( &fparams[j] ) ;
	printf( "RESULT_%zu :: %f +/- %f\n" , j , fparams[j].avg , fparams[j].err ) ;
      }
    }
  }

  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_fitfunc2( fparams , bootavg[i] , xavg[i] , mominfo[i] , fit , 
		   NPARAMS ,
		   INPARAMS -> quarks[i] , 
		   INPARAMS -> fit_lo , 
		   INPARAMS -> fit_hi ,
		   LT , INPARAMS -> NDATA[i] , NPARAMS ) ;
  }

  double xx[ 1 ] = { 0.0 } ;
  plot_data( fparams , xx , 1 ) ;
  graph_reset_color( ) ;

  free( fparams ) ;

  close_xmgrace_graph(  ) ;

  return ;
}
