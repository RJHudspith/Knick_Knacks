/**
   Analysis for the VPF
 */
#include "fitfunc.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "gevp.h"
#include "Utils.h"

#define GEVP_EFFMASS

// dispersions computation
void
gevp_eval( double **xavg ,
	   struct resampled **bootavg ,
	   struct mom_info **mominfo ,
	   struct input_params *INPARAMS ,
	   const int NSLICES ,
	   const int N , // correlation matrix length
	   const int LT )
{
  if( NSLICES%(N*N) != 0 ) {
    printf( "[GEVP] For this code N must mod into NSLICES \n" ) ;
    printf( "[GEVP] Nslices :: %d , N :: %d \n" , NSLICES , N ) ;
    return ;
  }

  // allocate eigenvalues
  const int NSIMS = NSLICES / ( N * N ) ;
  const int NEVALUES = ( N * NSIMS ) ;

  struct resampled **evalues = malloc( NEVALUES * sizeof( struct resampled *) ) ;
  int k ;
  for( k = 0 ; k < NEVALUES ; k++ ) {
    evalues[ k ] = malloc( INPARAMS -> NDATA[0] * sizeof( struct resampled ) ) ;
    int t ;
    for( t = 0 ; t < INPARAMS -> NDATA[0] ; t++ ) {
      evalues[k][t].resampled = malloc( bootavg[0][0].NSAMPLES * sizeof( double ) ) ;
      // set to zero
      equate_constant( &evalues[k][t] , 0.0 , bootavg[0][0].NSAMPLES , 
		       bootavg[0][0].restype ) ;
    }
  }
 
  int t ;
  #pragma omp parallel for private(t)
  for( t = 0 ; t < INPARAMS->NDATA[0] ; t++ ) {

    int i ;
    for( i = 0 ; i < NSIMS ; i++ ) {

      // temporary matrix space
      double C1[ N * N ] ;
      double C0[ N * N ] ;
      
      int j , n ;
      for( j = 0 ; j < bootavg[0][t].NSAMPLES ; j++ ) {
	
	double re_evalues[ N ] ;
	
	// set C1 and C0
	for( n = 0 ; n < N*N ; n++ ) {
	  C0[ n ] = bootavg[ n + i*N*N ][t].resampled[j] ;
	  C1[ n ] = bootavg[ n + i*N*N ][0].resampled[j] ;
	}

	// compute eigenvalues
	solve_gevp( re_evalues , C0 , C1 , N ) ;
	
	// poke into solution
	for( n = 0 ; n < N ; n++ ) {
	  evalues[ n + i*N ][ t ].resampled[ j ] = re_evalues[ n ] ;
	}
	
	// and that is it
      }
    
      // compute the error on our distribution
      for( n = 0 ; n < N ; n++ ) {
	compute_err( &evalues[ n + i*N ][ t ] ) ;
      }
      //
    } // NSIMS 
  }

#ifdef GEVP_EFFMASS
  // fit effective mass
  struct resampled **effmass = effective_mass( (const struct resampled**)evalues ,
					       INPARAMS -> NDATA ,
					       NEVALUES , ATANH_EFFMASS ) ;

  // fit and plot
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)effmass , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NEVALUES , LT ) ;

  free( fitparams ) ;
  free_resampled_dp( effmass , NEVALUES , INPARAMS -> NDATA ) ;
#else

  // fit and plot
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)evalues , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NEVALUES , LT ) ;

  free( fitparams ) ;

#endif
  free_resampled_dp( evalues , NEVALUES , INPARAMS -> NDATA ) ;

  return ;
}

