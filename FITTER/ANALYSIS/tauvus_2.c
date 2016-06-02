/**
   Analysis for the VPF
 */
#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "Utils.h"

#define VplusA

// tauvus computation
void
tauvus2_eval( double **xavg ,
	       struct resampled **bootavg ,
	       struct mom_info **mominfo ,
	       struct input_params *INPARAMS ,
	       const int NSLICES ,
	       const int LT ,
	       const bool renormalise )
{
  printf( "\n--> J0 evaluation <--\n" ) ;

  // tell us if we have given the wrong arguments
  if( NSLICES != 2 ) {
    printf( "Expected V(ls) A(ls)\n" ) ;
    return ;
  }

  printf( "\n--> Converting to physical momenta <--\n" ) ;

  int j , i ;
  // OK, momentum first
  for( j = 0 ; j < NSLICES ; j++ ) {
    // ainverse multipliers
    const double aI  = 1.0 ; //INPARAMS -> quarks[j].ainverse ;
    const double aI2 = pow( aI , 2 ) ;
    const double aI4 = pow( aI , 4 ) ;
    const double aI6 = pow( aI , 6 ) ;

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

  for( i = 0 ; i < NSLICES ; i++ ) {
    const double ZV = INPARAMS -> quarks[i].ZV ;
    for( j = 0 ; j < INPARAMS -> NDATA[ i ] ; j++ ) {
      mult_constant( &bootavg[i][j] , ZV ) ;
    }
  }

  char str[ 256 ] ;
  sprintf( str , "m%g_m%g.dat" , INPARAMS -> quarks[0].ml , INPARAMS -> quarks[0].ms ) ;

  FILE *file = fopen( str , "w" ) ;

  sprintf( str , "m%g_m%g_frac.dat" , INPARAMS -> quarks[0].ml , INPARAMS -> quarks[0].ms ) ;
  FILE *fracfile = fopen( str , "w" ) ;

  // so the idea is to form all contributions into bootavg[0]
  const double ZV = INPARAMS -> quarks[0].ZV ;

  fprintf( file , "ZV :: %f\n" , INPARAMS -> quarks[0].ZV ) ;
  fprintf( file , "a^{-1} :: %f GeV\n" , INPARAMS -> quarks[0].ainverse ) ;
  fprintf( file , "Combination :: ls(%g %g) (V+A;0+1)\n" , 
	   INPARAMS -> quarks[0].ml , INPARAMS -> quarks[0].ms ) ;
  fprintf( file , "q [GeV]        Pi(q^2)             Err\n" ) ;

  printf( "Renormalising with %f \n" , ZV ) ;


#ifdef VplusA
  // form flavour breaking combination
  for( j = 0 ; j < INPARAMS->NDATA[0] ; j++ ) {

    // form V+A (0+1) -ls and put in index 2
    add( &bootavg[0][j] , bootavg[1][j] ) ;

    // compute fractional error
    fprintf( fracfile , "FRAC :: %1.12E %1.12E \n" , xavg[0][j] , 100. * fabs( bootavg[0][j].err / bootavg[0][j].avg ) ) ;
    fprintf( file , "%f %1.12E %1.12E \n" , xavg[0][j] , bootavg[0][j].avg , bootavg[0][j].err ) ;
  }
  // plot only the flavour breaking difference
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , 1 , LT ) ;
#else
  /*
  for( i = 0 ; i < NSLICES ; i++ ) {
    for( j = 0 ; j < INPARAMS->NDATA[0] ; j++ ) {
      fprintf( fracfile , "FRAC :: %1.12E %1.12E \n" , xavg[0][j] , 100. * fabs( bootavg[i][j].err / bootavg[i][j].avg ) ) ;
      // multiply by Q^2 and ZV
      mult_constant( &bootavg[i][j] , xavg[0][j] * xavg[0][j] ) ;
      fprintf( file , "%f %1.12E %1.12E \n" , xavg[0][j] , bootavg[i][j].avg , bootavg[i][j].err ) ;
    }
    fprintf( fracfile , "FRAC :: \n" ) ;
    fprintf( file , "\n" ) ;
  }
  */

  // plot only the flavour breaking difference
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NSLICES , LT ) ;
#endif

  // this should be a loop over all possible lopos & hipos
  double a_k[ INPARAMS->NDATA[0] ] ;
  const size_t NPOLES = 2 ;

  // compute a_k
  for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
    register long double prod = 1.0 ;
    size_t i ;
    for( i = j ; i < (j+NPOLES) ; i++ ) {
      prod *= ( i != j ) ? ( xavg[0][i] * xavg[0][i] - xavg[0][j] * xavg[0][j] ) : 1.0 ;
      //printf( "%zu %e \n" , i , (double)prod ) ;
    }
    a_k[ j ] = 1.0 / (double)( prod ) ;
    //printf( "%f %e \n" , xavg[0][j] , a_k[j] ) ;
    mult_constant( &bootavg[0][j] , a_k[j] ) ;
    printf( "%f %f %f \n" , xavg[0][j] , bootavg[0][j].avg , bootavg[0][j].err ) ;
  } 

  
  fclose( fracfile ) ;
  fclose( file ) ;
  free( fitparams ) ;

  return ;
}

  /*
  // get the low and high pos
  const int lopos = find_idx( INPARAMS->fit_lo , xavg[0] , INPARAMS->NDATA[0] , 0 ) ;
  const int hipos = find_idx( INPARAMS->fit_hi , xavg[0] , INPARAMS->NDATA[0] , 0 ) ;

  printf( "lo,hi :: %d %d\n" , lopos, hipos) ;

  // compute a_k
  register double sum2 = 0.0 ;
  for( j = lopos ; j < hipos ; j++ ) {
    register long double prod = 1.0 ;
    size_t i ;
    for( i = lopos ; i < hipos ; i++ ) {
      prod *= ( i != j ) ? ( xavg[0][i] * xavg[0][i] - xavg[0][j] * xavg[0][j] ) : 1.0 ;
    }
    a_k[ j ] = 1.0 / ( prod ) ;
    printf( "%1.12f \n" , a_k[j] ) ;
    sum2 += ( a_k[ j ] ) * pow( xavg[0][i] , 0 ) ;
  } 
  printf( "\n--> should be zero :: %e <--\n\n" , sum2 ) ;

  equate_constant( &sum , 0.0 , bootavg[0][0].NSAMPLES , bootavg[0][0].restype ) ;
  for( j = lopos ; j < hipos ; j++ ) {
    equate( &tmp , bootavg[0][j] ) ;
    printf( "%f %e %e \n" , xavg[0][j]*xavg[0][j] , a_k[j] , tmp.avg ) ;
    mult_constant( &tmp , a_k[j] ) ;
    add( &sum , tmp ) ;
  }

  printf( "SUM :: %e %e \n" , sum.avg , sum.err ) ; 
  */
