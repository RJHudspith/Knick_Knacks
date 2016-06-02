/**
   Analysis for the VPF fit ...
   code uses the full subtraction 

   -4\pi^2 ( \Pi(q^2) - \Pi(q_1^2) ) / log( q_1^2 / q_2^2 ) - 1
 */
#include "fitfunc.h"

#include "fit_and_plot.h"
#include "run_boots.h"
#include "D0_diff_v1.h"
#include "Utils.h"
#include "correlation.h"
#include "cylinder.h"
#include "equivalents.h"
#include "momentum_stuff.h"
#include "coefficients.h"

#define PRINTDATA

// Alpha_s computation
void
alphas_v1_eval( double **xavg ,
		struct resampled **bootavg ,
		struct mom_info **mominfo ,
		struct input_params *INPARAMS ,
		const int NSLICES ,
		const int LT ,
		const bool renormalise )
{
  printf( "\n--> Renormalisation scale <--\n" ) ;

  set_mu_D0_diff_v1( INPARAMS -> quarks[0].mu ) ;

  printf( "\n--> Recylinder <--\n" ) ;

  int j , i ;
  const double widths[3] = { 0.25 , 0.19 , 0.27 } ;
  recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , LT , 4 , widths ) ;

  // BOOT and X are freed in momavg
  struct mom_info **mavg ;
  struct resampled **BAVG = momentum_average( &mavg , INPARAMS , 
					      (const struct mom_info**)mominfo , 
					      (const struct resampled**)bootavg , 
					      NSLICES , PSQ_AVERAGE ) ;
  // set up x-avg
  double **x = malloc( NSLICES * sizeof( double* ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    x[ i ] = ( double* )malloc( INPARAMS -> NDATA[i] * sizeof( double ) ) ;
    int k ;
    for( k = 0 ; k < INPARAMS -> NDATA[ i ] ; k++ ) {
      x[ i ][ k ] = mavg[ i ][ k ].p2 ;
    }
  }

  printf( "\n--> Converting to physical momenta <--\n" ) ;

  convert_to_physmom( mavg , x , INPARAMS , NSLICES ) ;

  printf( "\n--> Creating Delta <--\n" ) ;

  // renormalise 
  for( j = 0 ; j < NSLICES ; j++ ) {

    // is this guaranteed to be lower or just closest
    const int Q2_idx = find_idx( INPARAMS -> quarks[j].mu , x[j] , 
				 INPARAMS -> NDATA[j]-1 , 0 ) ;

    //
    struct resampled tmp ;
    tmp.resampled = malloc( BAVG[j][Q2_idx].NSAMPLES * sizeof( double ) ) ;

    equate( &tmp , BAVG[j][Q2_idx] ) ;
    set_q2q3_D0_diff_v1( mavg[j][Q2_idx].p2 , j ) ;

    //const double xref = x[j][Q2_idx] ;
    const double xxref = x[j][Q2_idx]*x[j][Q2_idx] ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      // perform a subtraction
      subtract( &BAVG[j][i] , tmp ) ;
      // set t1 - t2 factor
      const double t1mt2 = log( x[j][i] * x[j][i] / xxref ) ;
      const double factor = 4.0 * M_PI * M_PI / ( fabs( t1mt2 ) > 1E-12 ? t1mt2 : 1.0 ) ;
      mult_constant( &BAVG[j][i] , factor ) ;
      //subtract_constant( &BAVG[j][i] , 1.0 ) ;
      //x[j][i] -= xref ;
    }

    // free him up
    free( tmp.resampled ) ;
  }
  
  // fit and plot are in here
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)BAVG , 
						    (const double**)x , 
						    (const struct mom_info **)mavg ,
						    *INPARAMS , NSLICES , LT ) ;

  // use the fit result to get alpha_s at Mz
  if( INPARAMS->fittype != NOFIT ) {
      
    size_t mu ;
    for( mu = 0  ; mu < NSLICES ; mu++ ) {
      
      // should be sorted?
      if( fitparams[0].avg < 0.39 && 
	  fitparams[0].err_hi < 0.42 && 
	  fitparams[0].err_lo > 0.0 ) {
	printf( "--> Running results <--\n" ) ;
	struct resampled alpha = init_dist( fitparams + mu , fitparams[0].NSAMPLES , fitparams[0].restype ) ;
	struct resampled amz = init_dist( NULL , fitparams[0].NSAMPLES , fitparams[0].restype ) ;
	
	equate( &alpha , fitparams[ mu ] ) ;
	
	printf( "a(%1.2fGeV,3f) %f (%f) \n" , INPARAMS->quarks[mu].mu , alpha.avg , alpha.err ) ;
	
	amz = boot_run_MS_quick( INPARAMS->quarks[mu].mu , alpha , 5 , 3 , true ) ;
	
	printf( "a(M_Z,5f)_%zu %f %f %f %f \n" , mu ,
		INPARAMS -> fit_lo , INPARAMS -> fit_hi , 
		amz.avg , amz.err ) ;
	
	free( alpha.resampled ) ;
	free( amz.resampled ) ;
      } else {
	printf( "a(M_Z,5f)_%zu %f %f 1.0 1.0 \n" , mu ,
		INPARAMS -> fit_lo , INPARAMS -> fit_hi ) ;
      }
    }
  }

  free( fitparams ) ;

  return ;
}



