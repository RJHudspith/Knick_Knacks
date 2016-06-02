/**
   @file bk_ss_extrap.c
   @brief linear extrapolation in the strange, write out the result
 */
#include "fitfunc.h"

#include "fit_and_plot.h"
#include "GLU_bswap.h"
#include "svd.h"
#include "write_distribution.h"
#include "graph_data.h"
#include "fit_chooser.h"

void
ss_extrap_eval( double **xavg ,
		struct resampled **bootavg ,
		struct mom_info **mominfo ,
		struct mom_info *moms ,
		struct input_params *INPARAMS ,
		const int NSLICES ,
		const int LT )
{
  if( INPARAMS->fittype != POLY1 ) {
    printf( "[SS-EXTRAP] expected to be linear!\n" ) ;
    return ;
  }

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS -> graph_name ,
		      INPARAMS -> graph_xaxis , 
		      INPARAMS -> graph_yaxis ) ;

  // plot the data first but keep the graph file open
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_data( bootavg[i] , xavg[i] , INPARAMS -> NDATA[i] ) ;
  }

  // allocate fit params
  struct resampled *fparams = malloc( 2 * sizeof( struct resampled ) ) ;

  // linearly solve using the svd
  {
    size_t j ;

    // loop NSLICES
    for( i = 0 ; i < 1 ; i++ ) {

      const int range = INPARAMS -> NDATA[ i ] ;

      // set fparams to the same as bootparams
      for( j = 0 ; j < 2 ; j++ ) {
	fparams[j].resampled = malloc( bootavg[i][0].NSAMPLES * sizeof( double ) ) ;
	fparams[j].restype = bootavg[i][0].restype ;
	fparams[j].NSAMPLES = bootavg[i][0].NSAMPLES ;
      }

      // set the average first
      {
	// set the data
	double ydata[ range ] , xs[ range ] , sigma[ range] , coeffs[ 2 ] , chisq ;
	int k ;
	for( k = 0 ; k < range ; k++ ) {
	  ydata[ k ] = bootavg[ i ][ k ].avg ;
	  sigma[ k ] = bootavg[ i ][ k ].err ;
	  xs[ k ] = xavg[ i ][ k ] ;
	}
	//
	compute_coefficients( coeffs , &chisq , ydata , sigma , 
			      xs , range , 2 ) ;
	fparams[ 0 ].avg = coeffs[ 0 ] ;
	fparams[ 1 ].avg = coeffs[ 1 ] ;
      }

      // loop boots
      for( j = 0 ; j < bootavg[i][0].NSAMPLES ; j++ ) {

	// set the data
	double ydata[ range ] , xs[ range ] , sigma[ range ] , coeffs[ 2 ] , chisq ;

	int k ;
	for( k = 0 ; k < range ; k++ ) {
	  ydata[ k ] = bootavg[ i ][ k ].resampled[ j ] ;
	  sigma[ k ] = bootavg[ i ][ k ].err ;
	  xs[ k ] = xavg[ i ][ k ] ;
	}
	//

	compute_coefficients( coeffs , &chisq , ydata , sigma , 
			      xs , range , 2 ) ;

	fparams[ 0 ].resampled[ j ] = coeffs[ 0 ] ;
	fparams[ 1 ].resampled[ j ] = coeffs[ 1 ] ;
      }
      for( j = 0 ; j < 2 ; j++ ) {
	compute_err( &fparams[j] ) ;
	printf( "RESULT_%zu :: %f +/- %f\n" , j , fparams[j].avg , fparams[j].err ) ;
      }
    }
  }

  /*
  for( i = 0 ; i < NSLICES ; i++ ) {
    fitfunc fit ;
    initialise_f( POLY1 , &fit , INPARAMS->NDATA , 0 ) ;
    plot_fitfunc2( fparams , bootavg[i] , xavg[i] , mominfo[i] , fit , 
		   INPARAMS -> NDATA[0] ,
		   INPARAMS -> quarks[i] , 
		   INPARAMS -> fit_lo , 
		   INPARAMS -> fit_hi ,
		   LT , INPARAMS -> NDATA[i] , i ) ;
  }
  */
  
  double xx[ 1 ] = { 0.0 } ;
  plot_data( fparams , xx , 1 ) ;
  graph_reset_color( ) ;

  printf( "\n--> Writing out the extrap to %s <--\n" , INPARAMS->output_file ) ;

  // write out the mass to a file
  FILE *outfile = fopen( INPARAMS->output_file , "wb" ) ;

  // write the distribution and its size
  write_singledist( fparams[0] , outfile ) ;

  fclose( outfile ) ;

  free( fparams ) ;

  close_xmgrace_graph(  ) ;

  return ;
}
