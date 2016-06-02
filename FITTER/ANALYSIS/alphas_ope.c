/**
   Analysis for the VPF fit ...
 */
#include "fitfunc.h"

#include "fit_and_plot.h"
#include "run_boots.h"
#include "D0_diff.h"
#include "D0_diff_multi.h"
#include "Utils.h"
#include "correlation.h"
#include "cylinder.h"
#include "equivalents.h"
#include "momentum_stuff.h"

#define PRINTDATA

// Alpha_s computation
void
alphas_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ,
	     const bool renormalise )
{
  printf( "\n--> Renormalisation scale <--\n" ) ;

  set_mu_D0_diff( INPARAMS -> quarks[0].mu ) ;

  printf( "\n--> Recylinder <--\n" ) ;

  int j , i ;
  const double widths[3] = { 0.25 , 0.25 , 0.25 } ;
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

  printf( "\n--> Computing correlation matrix <--\n" ) ;
#ifdef CORR_MATRIX
  double **correlation = malloc( INPARAMS -> NDATA[0] * sizeof( double* ) ) ;
  for( i = 0 ; i < INPARAMS -> NDATA[0] ; i++ ) {
    correlation[i] = malloc( INPARAMS -> NDATA[0] * sizeof( double ) ) ;
  }  
  correlations( correlation , BAVG[0] , INPARAMS -> NDATA[0] ) ;

  write_corrmatrix_mathematica( correlation , INPARAMS -> NDATA[0] ) ;

  write_corrmatrix_mathematica( correlation , INPARAMS -> NDATA[0] ) ;
  exit(1) ; 
#endif

  printf( "\n--> Converting to physical momenta <--\n" ) ;

  convert_to_physmom( mavg , x , INPARAMS , NSLICES ) ;

  printf( "\n--> Creating Delta <--\n" ) ;

  // renormalise 
  for( j = 0 ; j < NSLICES ; j++ ) {

    // is this guaranteed to be lower or just closest
    const int Q2_idx = find_idx( INPARAMS -> fit_lo , x[j] , 
				 INPARAMS -> NDATA[j]-1 , 0 ) ;

    //
    struct resampled tmp ;
    tmp.resampled = malloc( BAVG[j][Q2_idx].NSAMPLES * sizeof( double ) ) ;

    equate( &tmp , BAVG[j][Q2_idx] ) ;
    set_q2q3_D0_diff( mavg[j][Q2_idx].p2 , j ) ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      subtract( &BAVG[j][i] , tmp ) ;
    }

    printf( "{" ) ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      printf( "%1.10f, \n" , BAVG[j][i].err ) ;
    }
    printf( "}\n" ) ;

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
    printf( "--> Running results <--\n" ) ;

    size_t mu ;
    for( mu = 0  ; mu < NSLICES ; mu++ ) {
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
    }
  }

  free( fitparams ) ;

  return ;
}


















    // test the lerp again
    /*
    equate( &tmp , BAVG[j][Q2_idx+1] ) ;
    subtract( &tmp , BAVG[j][Q2_idx+1] ) ;
    mult_constant( &tmp , 
		   ( INPARAMS -> quarks[j].mu - x[j][Q2_idx] ) /
		   ( x[j][Q2_idx+1] - x[j][Q2_idx] ) ) ;
    add( &tmp , BAVG[j][Q2_idx] ) ;
    set_q2q3_D0_diff( INPARAMS -> quarks[j].mu , j ) ;
    */







    // kims original idea
#if 0
    size_t step ;
    for( step = 5 ; step < 6 ; step++ ) {
      for( i = 2 ; i < 25 ; i++ ) {
	struct resampled temp = init_dist( &BAVG[j][i] , BAVG[j][i].NSAMPLES , BAVG[j][i].restype ) ;
	subtract( &temp , BAVG[j][i+step] ) ;
	divide_constant( &temp , log( (x[j][i]*x[j][i])/(x[j][i+step]*x[j][i+step])) ) ;
	mult_constant( &temp , 4.0*M_PI*M_PI ) ;
	subtract_constant( &temp , 1.0 ) ;

	const double mu = 2 ; //x[j][i] + ( x[j][i] + x[j][i+step] ) * 0.5 ;

	struct resampled alpha = alpha_solve( temp , log(x[j][i]*x[j][i]/(mu*mu) )  ,
					      log(x[j][i+step]*x[j][i+step]/(mu*mu) ) ) ;

	//struct resampled amz = init_dist( &alpha , alpha.NSAMPLES , alpha.restype ) ;
    
	//amz = boot_run_MS_quick( mu , alpha , 5 , 3 , true ) ;
	printf( "%f %f %f\n" , x[j][i] , alpha.avg , alpha.err ) ;
      }
      printf( "\n" ) ;
    }


	//printf( "{%1.10f,%1.10f}, \n" , x[j][i] , BAVG[j][i].avg ) ;
     /*
      divide_constant( &BAVG[j][i] , x[j][i] * x[j][i] * 
		       ( i != Q2_idx ? log( x[j][i]*x[j][i] / ( xxref ) ) : 1.0 ) 
		       ) ;
      mult_constant( &BAVG[j][i] , 4.0*M_PI*M_PI ) ;
      subtract_constant( &BAVG[j][i] , 1.0 ) ;
      x[j][i] = x[j][i] - xref ;
      */

  /*
  // subtract a from b
#pragma omp parallel for private(i)
  for( i = 0 ; i < INPARAMS -> NDATA[0] ; i++ ) {
    subtract( &BAVG[0][i] , BAVG[1][i] ) ;
  }
  */
#endif


