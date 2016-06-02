/**
   Analysis for the VPF
 */
#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"

// tauvus computation
void
flavour_combination_eval( double **xavg ,
			  struct resampled **bootavg ,
			  struct mom_info **mominfo ,
			  struct input_params *INPARAMS ,
			  const int NSLICES ,
			  const int LT ,
			  const bool renormalise )
{
  printf( "\n--> Flavour combination evaluation <--\n" ) ;

  // tell us if we have given the wrong arguments
  if( NSLICES != 4 ) {
    printf( "Expected V (ss) VV(ls) A(ll) A(ls)\n" ) ;
    return ;
  }

  printf( "\n--> Converting to physical momenta <--\n" ) ;
  // OK, momentum first
  int j , i ;
  for( j = 0 ; j < NSLICES ; j++ ) {
    // ainverse multipliers
    const double aI  = INPARAMS -> quarks[j].ainverse ;
    const double aI2 = pow( INPARAMS -> quarks[j].ainverse , 2 ) ;
    const double aI4 = pow( INPARAMS -> quarks[j].ainverse , 4 ) ;
    const double aI6 = pow( INPARAMS -> quarks[j].ainverse , 6 ) ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {

      // set to physical lattice spacing xavg[j][i] is now |p| in GeV!!
      xavg[j][i] = aI * sqrt( xavg[j][i] ) ;
      // and the moms
      mominfo[j][i].p2 *= aI2 ;
      mominfo[j][i].p4 *= aI4 ;
      mominfo[j][i].p6 *= aI6 ;

      // check for some consistency
      if( fabs( mominfo[j][i].p2 - xavg[j][i]*xavg[j][i] ) > 1E-12 ) {
	printf( "P2 Broken !! %e \n" , mominfo[j][i].p2 - xavg[j][i]*xavg[j][i] ) ;
	exit(1) ;
      } 
    }
  }

  // multiply by Z_V
  printf( "\n--> Renormalising <--\n" ) ;

  char str[ 256 ] ;
  sprintf( str , "m%g_m%g.dat" , INPARAMS -> quarks[1].ml , INPARAMS -> quarks[1].ms ) ;

  FILE *file = fopen( str , "w" ) ;

  // so the idea is to form all contributions into bootavg[0]
  const double ZV = INPARAMS -> quarks[0].ZV ;

  fprintf( file , "ZV :: %f\n" , INPARAMS -> quarks[0].ZV ) ;
  fprintf( file , "a^{-1} :: %f GeV\n" , INPARAMS -> quarks[0].ainverse ) ;
  fprintf( file , "Combination :: ( ss - ls )_V + ( ll - ls )_A (0+1) \n" ) ;
  fprintf( file , "q [GeV]        Pi(q^2)             Err\n" ) ;

  printf( "Renormalising with %f \n" , ZV ) ;

  // form flavour breaking combination
  for( j = 0 ; j < INPARAMS->NDATA[0] ; j++ ) {

    // form V(0+1) -(ss-ls) put in index 0
    subtract( &bootavg[0][j] , bootavg[1][j] ) ;

    // form A (0+1)-(ll-ls)  and put in index 2
    subtract( &bootavg[2][j] , bootavg[3][j] ) ;

    // and compute the flavour breaking difference
    add( &bootavg[0][j] , bootavg[2][j] ) ;

    // and then renormalise
    mult_constant( &bootavg[0][j] , ZV ) ;

    fprintf( file , "%f %1.12E %1.12E \n" , xavg[0][j] , bootavg[0][j].avg , bootavg[0][j].err ) ;
  }

  // plot only the flavour breaking difference
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , 1 , LT ) ;
  
  fclose( file ) ;
  free( fitparams ) ;

  return ;
}
