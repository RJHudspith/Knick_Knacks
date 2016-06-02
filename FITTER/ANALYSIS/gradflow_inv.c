/**
   @file gradflow_inv.c
   @brief infinite volume extrap of gradient flow
 */
#include "fitfunc.h"

#include "graph_data.h"
#include "svd.h"
#include "Utils.h"

#if 0
int
gradflow_inv( double **X ,
	      struct resampled **BOOT ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT )
{
  // for each momentum extrapolate linearly in 1/V
 
  // precompute volume factors
  double _V[ NSLICES ] ;
  size_t i ;
  for( i = 0; i < NSLICES ; i++ ) {
    _V[ i ] = 1.0 / ( INPARAMS -> dimensions[ i ][ 0 ] *
		      INPARAMS -> dimensions[ i ][ 1 ] * 
		      INPARAMS -> dimensions[ i ][ 2 ] * 
		      INPARAMS -> dimensions[ i ][ 3 ] ) ;
    printf( "INVOLS :: %f \n" , _V[i] ) ;
  }

  struct resampled *infvol = malloc( INPARAMS -> NDATA[ 0 ] * sizeof( struct resampled ) ) ;


  // set the range based on fit_hi and fit_lo
  const int lopos = find_idx( INPARAMS->fit_lo , _V , NSLICES , 0 ) ;
  const int hipos = find_idx( INPARAMS->fit_hi , _V , NSLICES , 0 ) ;
  const int range = hipos - lopos + 1 ;

  printf( "LOPOS :: %d %f \n" , lopos , _V[ lopos ] ) ;
  printf( "HIPOS :: %d %f \n" , hipos , _V[ hipos ] ) ;
  printf( "range :: %d \n" , range ) ;

  if( range < 2 ) {
    printf( "Fit range does not allow for linear extrapolation\n" ) ;
    return FAILURE ;
  }

  // loop samples
  size_t j ;
#pragma omp parallel for private(j)
  for( j = 0 ; j < INPARAMS -> NDATA[ 0 ] ; j++ ) {
    infvol[ j ].resampled = malloc( BOOT[0][j].NSAMPLES * sizeof( double ) ) ;
    equate_constant( &infvol[j] , 0.0 , BOOT[0][j].NSAMPLES , BOOT[0][j].restype ) ;

    double ydata[ range ] , xdata[ range ] , sigma[ range ] ;
    double coeffs[ 2 ] ;

    // loop samples
    size_t l , k , idx ;
    for( l = 0 ; l < BOOT[0][j].NSAMPLES ; l++ ) {
      // set y data
      idx = 0 ;
      for( k = lopos ; k <= hipos ; k++ ) {
	ydata[ idx ] = BOOT[k][j].resampled[l] ;
	sigma[ idx ] = BOOT[k][j].err ;
	xdata[ idx ] = _V[ k ] ;
	idx ++ ;
      }
      compute_coefficients( coeffs , ydata , sigma , xdata , range , 2 ) ;
      infvol[j].resampled[l] = coeffs[0] ;
    }
    
    // set y data
    idx = 0 ;
    for( k = lopos ; k <= hipos ; k++ ) {
      ydata[ idx ] = BOOT[k][j].avg ;
      sigma[ idx ] = BOOT[k][j].err ;
      xdata[ idx ] = _V[ k ] ;
      idx++ ;
    }
    compute_coefficients( coeffs , ydata , sigma , xdata , range , 2 ) ;
    infvol[j].avg = coeffs[0] ;

    // make sure we get the errors correct
    compute_err( &infvol[j] ) ;
  }
  printf( "Graphs?\n" ) ;

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS -> graph_name ,
		      INPARAMS -> graph_xaxis , 
		      INPARAMS -> graph_yaxis ) ;

  // plot the data first but keep the graph file open
  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_data( BOOT[i] , X[i] , INPARAMS -> NDATA[i] ) ;
  }
  plot_data( infvol , X[0] , INPARAMS -> NDATA[0] ) ;

  close_xmgrace_graph(  ) ;

  return SUCCESS ;
}
#endif

int
gradflow_inv( double **X ,
	      struct resampled **BOOT ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT )
{
  int DOF = 1 ;
  // for each momentum extrapolate linearly in 1/V
  switch( INPARAMS -> fittype ) {
  case POLY0 :
    DOF = 1 ; printf( "CONSTANT FIT \n" ) ;
    break ;
  case POLY1 :
    DOF = 2 ; printf( "LINEAR FIT \n" ) ;
    break ;
  case POLY2 :
    DOF = 3 ; printf( "QUADRATIC FIT \n" ) ;
    break ;
  case POLY3 :
    DOF = 4 ; printf( "CUBIC FIT \n" ) ;
    break ;
  default :
    printf( "POLY not recognised\n" ) ;
    return FAILURE ;
  }
 
  // precompute volume factors
  double _V[ NSLICES ] ;
  size_t i ;
  for( i = 0; i < NSLICES ; i++ ) {
    _V[ i ] = 1.0 / ( INPARAMS -> quarks[i].ainverse ) ;
    printf( "INVOLS :: %f \n" , _V[i] ) ;
  }

  struct resampled *infvol = malloc( INPARAMS -> NDATA[ 0 ] * sizeof( struct resampled ) ) ;
  struct resampled *chi = malloc( INPARAMS -> NDATA[ 0 ] * sizeof( struct resampled ) ) ;

  // set the range based on fit_hi and fit_lo
  const int lopos = find_idx( INPARAMS->fit_lo , _V , NSLICES , 0 ) ;
  const int hipos = find_idx( INPARAMS->fit_hi , _V , NSLICES , 0 ) ;
  const int range = hipos - lopos + 1 ;

  printf( "LOPOS :: %d %f \n" , lopos , _V[ lopos ] ) ;
  printf( "HIPOS :: %d %f \n" , hipos , _V[ hipos ] ) ;
  printf( "range :: %d \n" , range ) ;

  if( range < 2 ) {
    printf( "Fit range does not allow for linear extrapolation\n" ) ;
    return FAILURE ;
  }

  // loop samples
  size_t j ;
  //#pragma omp parallel for private(j)
  for( j = 0 ; j < INPARAMS -> NDATA[ 0 ] ; j++ ) {
    infvol[ j ].resampled = malloc( BOOT[0][j].NSAMPLES * sizeof( double ) ) ;
    equate_constant( &infvol[j] , 0.0 , BOOT[0][j].NSAMPLES , BOOT[0][j].restype ) ;

    chi[ j ].resampled = malloc( BOOT[0][j].NSAMPLES * sizeof( double ) ) ;
    equate_constant( &chi[j] , 0.0 , BOOT[0][j].NSAMPLES , BOOT[0][j].restype ) ;

    double ydata[ range ] , xdata[ range ] , sigma[ range ] , chisq = 0.0 ;
    double coeffs[ DOF ] ;

    // set y data
    size_t k , idx ;
    idx = 0 ;
    for( k = lopos ; k <= hipos ; k++ ) {
      ydata[ idx ] = BOOT[k][j].avg ;
      sigma[ idx ] = BOOT[k][j].err ;
      xdata[ idx ] = _V[ k ] ;
      idx++ ;
    }
    compute_coefficients( coeffs , &chisq , ydata , xdata , 
			  sigma , range , DOF ) ;
    printf( "WHAT ? %e \n" , chisq ) ;
    chi[j].avg = chisq / DOF ;
    infvol[j].avg = coeffs[0] ;

    size_t l ;
    // loop samples
    for( l = 0 ; l < BOOT[0][j].NSAMPLES ; l++ ) {
      // set y data
      idx = 0 ;
      for( k = lopos ; k <= hipos ; k++ ) {
	ydata[ idx ] = BOOT[k][j].resampled[l] ;
	xdata[ idx ] = _V[ k ] ;
	idx ++ ;
      }
      compute_coefficients( coeffs , &chisq , ydata , sigma , 
			    xdata , range , DOF ) ;
      chi[j].resampled[l] = chisq / DOF ;
      infvol[j].resampled[l] = coeffs[0] ;
    }

    // make sure we get the errors correct
    compute_err( &chi[j] ) ;
    compute_err( &infvol[j] ) ;

    printf( "CHISQ :: %f %f +/- %f \n" , X[0][j] , chi[j].avg , chi[j].err ) ;
    printf( "EXTRAP :: %f %f +/- %f \n" , X[0][j] , infvol[j].avg , infvol[j].err ) ;
  }

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS -> graph_name ,
		      INPARAMS -> graph_xaxis , 
		      INPARAMS -> graph_yaxis ) ;

  // plot the data first but keep the graph file open
  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_data( BOOT[i] , X[i] , INPARAMS -> NDATA[i] ) ;
  }
  plot_data( infvol , X[0] , INPARAMS -> NDATA[0] ) ;

  /*
  for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
    printf( "%f \n" , X[0][j] ) ;
    for( i = 0 ; i < NSLICES ; i++ ) {
      //res_log( &BOOT[i][j] ) ;
      printf( "%e %e %e \n" , log( 1.0 / INPARAMS->quarks[i].ainverse ) , 
	      BOOT[i][j].avg , BOOT[i][j].err ) ;
    }
    printf( "\n" ) ;
  }
  */

  close_xmgrace_graph(  ) ;

  return SUCCESS ;
}
