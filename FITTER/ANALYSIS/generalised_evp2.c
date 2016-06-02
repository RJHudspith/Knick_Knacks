/**
   Analysis for the VPF
 */
#include "fitfunc.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "gevp.h"
#include "Utils.h"

#define GEVP_EFFMASS

// compute the GEVP using ss, ls , sl , ll and maybe wl and ww i.e.
//
//  | ss  sl  0   0  |
//  | ls  ll  0   0  |
//  | 0   0   wl  0  |
//  | 0   0   0   ww |
//
void
gevp2_eval( double **xavg ,
	    struct resampled **bootavg ,
	    struct mom_info **mominfo ,
	    struct input_params *INPARAMS ,
	    const int NSLICES ,
	    const int N ,  // correlation matrix length
	    const int LT )
{
  const int N2 = NSLICES - N*N ;

  if( N2 < 0 ) {
    printf( "[GEVP2] input doesn't make sense %d < %d \n" , NSLICES , N*N ) ;
    return ;
  }

  printf( "[GEVP2] %d blocked and %d diagonal \n" , N , N2 ) ;

  // allocate eigenvalues
  const int NEVALUES = ( N + N2 ) ;

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
  //#pragma omp parallel for private(t)
  for( t = 0 ; t < INPARAMS->NDATA[0] ; t++ ) {

    // temporary matrix space
    double C1[ NEVALUES*NEVALUES ] ;
    double C0[ NEVALUES*NEVALUES ] ;
      
    size_t j , n ;
    for( j = 0 ; j < bootavg[0][t].NSAMPLES ; j++ ) {

      // initialise to 0
      for( n = 0 ; n < NEVALUES*NEVALUES ; n++ ) {
	C1[ n ] = C0[ n ] = 0.0 ;
      }
	
      // square matrix bit
      for( n = 0 ; n < N ; n++ ) {
	size_t n2 ;
	for( n2 = 0 ; n2 < N ; n2++ ) {
	  C0[ n2 + n * NEVALUES ] = bootavg[ n2 + n * N ][t].resampled[j] ;
	  C1[ n2 + n * NEVALUES ] = bootavg[ n2 + n * N ][0].resampled[j] ;
	}
      }

      // excess bit is diagonal
      size_t idx = N*N ;
      for( n = N ; n < NEVALUES ; n++ ) {
	C0[ n * ( NEVALUES + 1 ) ] = bootavg[ idx ][t].resampled[j] ;
	C1[ n * ( NEVALUES + 1 ) ] = bootavg[ idx ][0].resampled[j] ;
	idx++ ;
      }

      // compute eigenvalues
      double re_evalues[ NEVALUES ] ;

      solve_gevp( re_evalues , C0 , C1 , NEVALUES ) ;
	
      // poke into solution
      for( n = 0 ; n < NEVALUES ; n++ ) {
	evalues[ n ][ t ].resampled[ j ] = re_evalues[ n ] ;
      }
	
      // and that is it
    }
    
    // compute the error on our distribution
    for( n = 0 ; n < N+N2 ; n++ ) {
      compute_err( &evalues[ n ][ t ] ) ;
    }
    //
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

