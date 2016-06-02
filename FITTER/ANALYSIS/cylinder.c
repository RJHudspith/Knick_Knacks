#include "fitfunc.h"

#include "correlation.h"
#include "equivalents.h"
#include "graph_data.h"
#include "GLU_timer.h"

enum{ DO_NOT_ADD , ADD_TO_LIST } list_creation ;

////////// Cylinder cutting procedurals //////////
// gets the body diagonal vectors for our lattice
static inline void
get_diagonal( n , i , DIMS )
     double n[ 4 ] ;
     const int i , DIMS ;
{
  int mu , subvol = 1 ;
  for( mu = 0 ; mu < 4 ; mu++ ) {
    if( mu < DIMS ) {
      n[ mu ] = ( ( i - i % subvol ) / subvol ) % 2 ;
      if( n[ mu ] == 0 ) {
	n[ mu ] -- ;
      }
      subvol *= 2 ;
    } else {// set it to 0?
      n[ mu ] = 0 ;
    }
  }
  return ;
}

// generic cylinder calc
static int 
cylinder_DF( const double q[ 4 ] ,
	     const int DIMS ,
	     const double cyl_width )
{
  // test that this satisfies the correct cylinder cutting procedure
  const double norm = 1. / ( 2. * sqrt( DIMS ) ) ;
  const int diagonals = 2 << ( 4 - 1 ) ;
 
  // generix for loop over the diagonals
  int mu ;
  double x[ 4 ] ;
  for( mu = 0 ; mu < diagonals ; mu ++ ) {
    // inline for getting a diagonal for lexi order
    get_diagonal( x , mu , DIMS ) ;

    double scalar_prod = 0. ; 
    int nu ;
    for( nu = 0 ; nu < 4 ; nu ++ ) {
      scalar_prod += q[ nu ] * x[ nu ] ;
    }
    scalar_prod *= norm ;

    double mod = 0.0 ; 
    for( nu = 0 ; nu < 4 ; nu ++ ) {
      register const double temp = q[ nu ] - scalar_prod * x[ nu ] ;
      mod += temp * temp ;
    }

    if( sqrt( mod ) <= cyl_width ) {
      return ADD_TO_LIST ;
    }
  }
  return DO_NOT_ADD ;
}

/**
   @brief dirty hack prints out the cylinder cut data
 */
void
cylinder_VmA( struct resampled **BOOTS ,
	      double **X ,
	      struct mom_info **mominfo ,
	      const struct input_params INPARAMS ,
	      const int NSLICES )
{
  if( NSLICES != 2 ) {
    printf( "Only set up for 2 files \n" ) ;
    return ;
  }

  // precompute twiddles
  double twiddles[ 4 ] ;
  int mu ;
  for( mu = 0 ; mu < 4 ; mu++ ) {
    twiddles[ mu ] = 2.0 * M_PI / (double)INPARAMS.dimensions[0][mu] ; 
  }

  int *list = malloc( INPARAMS.NDATA[0] * sizeof( int ) ) ;
  int ncut = 0 ;

  // cylinder width
  const double width = 2.0 * twiddles[ 0 ] ;

  printf( "\n--> Renormalising <--\n" ) ;

  // first pass to get the length of the cut
  int i ;
  for( i = 0 ; i < INPARAMS.NDATA[0] ; i++ ) {

    X[0][i] *= INPARAMS.quarks[0].ainverse * INPARAMS.quarks[0].ainverse ;

    // compute V-A
    subtract( &BOOTS[0][i] , BOOTS[1][i] ) ;

    // renormalise
    mult_constant( &BOOTS[0][i] , INPARAMS.quarks[0].ZV ) ;

    double q[ 4 ] ;
    int mu ;
    for( mu = 0 ; mu < 4 ; mu++ ) {
      q[ mu ] = mominfo[0][i].n[mu] * twiddles[ mu ] ;
    }
    // cut is a cylinder and 
    if( cylinder_DF( q , 4 , width ) == ADD_TO_LIST ) {
      if( sqrt( X[0][i] ) < 4.5 ) {
	list[ ncut ] = i ;
	ncut++ ;
      }
    }
  }

  print_time( ) ;

  printf( "\n--> Cutting <--\n" ) ;

  // cut data
  struct resampled *cut = malloc( ncut * sizeof( struct resampled ) ) ;
  struct mom_info *momcut = malloc( ncut * sizeof( struct mom_info ) ) ;

  for( i = 0 ; i < ncut ; i++ ) {
    cut[i].resampled = (double*)malloc( BOOTS[0][i].NSAMPLES * sizeof( double ) ) ;
    equate( &cut[i] , BOOTS[0][ list[i] ] ) ;

    #ifdef PION_POLE_SUBTRACTION
    add_constant( &cut[i] , ( INPARAMS.quarks[0].f_pi * INPARAMS.quarks[0].f_pi ) 
		  / ( xcut[i] + INPARAMS.quarks[0].m_pi * INPARAMS.quarks[0].m_pi ) ) ;
    #endif

    memcpy( &momcut[i] , &mominfo[0][list[i]] , sizeof( struct mom_info ) ) ;
  }

  printf( "\n--> Momentum Average <--\n" ) ;

  int NAVE = ncut ;
  struct mom_info *momavg ;
  struct resampled *cutave = momavg_single( &momavg , &NAVE , 
					    momcut , cut ) ;

  double *xcut = malloc( NAVE * sizeof( double ) ) ;
  for( i = 0 ; i < NAVE ; i++ ) {
    xcut[ i ] = momavg[ i ].p2 * 
      INPARAMS.quarks[0].ainverse * 
      INPARAMS.quarks[0].ainverse ;
    printf( "%1.15f %1.15f %1.15f\n" , xcut[i] ,
	    cut[i].avg , cut[i].err ) ;
  }

  print_time( ) ;

#if 0
  // write it out
  printf( "\n--> [CORR] Writing out an %d X %d correlation matrix <--\n" , 
	  NAVE , NAVE ) ;

  // Kim likes correlation matrices
  double **correlation = malloc( NAVE * sizeof( double* ) ) ;
  for( i = 0 ; i < NAVE ; i++ ) {
    correlation[i] = (double*)malloc( NAVE * sizeof( double ) ) ;
  }
  correlations( correlation , cutave , NAVE ) ;

  write_corrmatrix( correlation , NAVE ) ;

  // and free it up
  for( i = 0 ; i < NAVE ; i++ ) {
    free( correlation[i] ) ;
  }
  free( correlation ) ;
#endif

  print_time( ) ;

  // make a plot
  make_xmgrace_graph( INPARAMS.graph_name ,
		      INPARAMS.graph_xaxis , 
		      INPARAMS.graph_yaxis ) ;

  plot_data( cutave , xcut , NAVE ) ;

  close_xmgrace_graph(  ) ;

  printf( "\nGraph plotted to %s\n" , INPARAMS.graph_name) ;

  print_time( ) ;

  // memfrees
  free( list ) ;

  free( xcut ) ;

  for( i = 0 ; i < ncut ; i++ ) {
    free( cut[i].resampled ) ;
  }
  free( cut ) ;

  for( i = 0 ; i < NAVE ; i++ ) {
    free( cutave[i].resampled ) ;
  }
  free( cutave ) ;
  free( momavg ) ;

  return ;
}
