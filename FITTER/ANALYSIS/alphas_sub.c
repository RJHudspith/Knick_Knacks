/**
   Subtract PT and look at data
 */
#include "fitfunc.h"

#include "coefficients.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "D0_diff.h"
#include "cylinder.h"
#include "equivalents.h"
#include "momentum_stuff.h"
#include "Utils.h"
#include "graph_data.h"

// Alpha_s computation
void
alphas_sub( double **xavg ,
	    struct resampled **bootavg ,
	    struct mom_info **mominfo ,
	    struct input_params *INPARAMS ,
	    const int NSLICES ,
	    const int LT ,
	    const bool renormalise )
{
  printf( "\n--> Alpha_s Evaluation <--\n" ) ;

  // set mu
  const double renscale = INPARAMS -> quarks[0].mu ;
  set_mu_D0_diff( renscale ) ;

  printf( "\n--> Recylinder <--\n" ) ;

  const double widths[3] = { 0.25 , 0.25 , 0.25 } ;
  recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , 
	      LT , 4 , widths ) ;

  // BOOT and X are freed in momavg
  struct mom_info **mavg ;
  struct resampled **BAVG = momentum_average( &mavg , INPARAMS , 
					      (const struct mom_info**)mominfo , 
					      (const struct resampled**)bootavg , 
					      NSLICES , PSQ_AVERAGE ) ;
  // set up x-avg
  double **x = malloc( NSLICES * sizeof( double* ) ) ;
  size_t i , j ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    x[ i ] = ( double* )malloc( INPARAMS -> NDATA[i] * sizeof( double ) ) ;
    int k ;
    for( k = 0 ; k < INPARAMS -> NDATA[ i ] ; k++ ) {
      x[ i ][ k ] = mavg[ i ][ k ].p2 ;
      printf( "%f %f %f \n" , mavg[i][k].p2 , BAVG[i][k].avg , BAVG[i][k].err ) ;
    }
  }

  printf( "\n--> Converting to physical momenta <--\n" ) ;

  convert_to_physmom( mavg , x , INPARAMS , NSLICES ) ;

  printf( "\n--> Performing PT subtraction <--\n" ) ;

  // renormalise 
  for( j = 0 ; j < NSLICES ; j++ ) {

    const int Q2_idx = find_idx( INPARAMS -> quarks[j].mu , x[j] , 
				 INPARAMS -> NDATA[j]-1 , 0 ) ;

    // fix this! hack for now
    set_q2q3_D0_diff( mavg[j][Q2_idx].p2 , j ) ;

    //
    struct resampled tmp ;
    tmp.resampled = malloc( BAVG[j][Q2_idx].NSAMPLES * sizeof( double ) ) ;
    equate( &tmp , BAVG[j][Q2_idx] ) ;

    // set the subtraction point
    const double XXref = x[j][Q2_idx]*x[j][Q2_idx] ;

    const double tref = log( ( XXref ) / ( renscale * renscale ) ) ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      subtract( &BAVG[j][i] , tmp ) ;

      const double t = log( ( x[j][i]*x[j][i] ) / ( renscale * renscale ) ) ;
      
      const double a = 0.28 / M_PI + 0.1 * pow( INPARAMS -> quarks[j].ainverse , -2.0 ) 
	+ 10 * pow( INPARAMS -> quarks[j].ainverse , -4.0 ) ;
      subtract_constant( &BAVG[j][i] ,
			 ( D0_OPE( t , a , 0 , 5 )
			   - D0_OPE( tref , a , 0 , 5 ) ) ) ;
      
    }

    // free him up
    free( tmp.resampled ) ;
  }

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS -> graph_name ,
		      INPARAMS -> graph_xaxis , 
		      INPARAMS -> graph_yaxis ) ;

  // plot the data first but keep the graph file open
  for( i = 0 ; i < NSLICES ; i++ ) {
    const int Q2_idx = find_idx( INPARAMS -> quarks[i].mu , x[i] , 
   				 INPARAMS -> NDATA[i]-1 , 0 ) ;
    const double XXref = pow( x[i][Q2_idx] , 2 ) ;
    for( j = 0 ; j < INPARAMS -> NDATA[i] ; j++ ) {
      //res_log( &BAVG[i][j] ) ;
      x[i][j] = ( x[i][j] * x[i][j] - XXref ) ;
    }
    plot_data( BAVG[i] , x[i] , INPARAMS -> NDATA[i] ) ;
  }
  graph_reset_color( ) ;

  // close up the graph
  close_xmgrace_graph(  ) ;

  return ;
}
