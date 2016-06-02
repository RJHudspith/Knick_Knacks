/**
   Analysis for the VPF fit ...
 */
#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "D0_diff.h"
#include "Utils.h"
#include "coefficients.h"
#include "graph_data.h"
#include "equivalents.h"
#include "shuffalgo.h"
#include "derivatives.h"
#include "fit_chooser.h"

#define PRINTDATA

enum{ DO_NOT_ADD , ADD_TO_LIST } list_creation ;

////////// Cylinder cutting procedurals //////////
// gets the body diagonal vectors for our lattice
static inline void
get_diagonal( ND , n , i , DIMS )
     const int ND ;
     double n[ ND ] ;
     const int i , DIMS ;
{
  int mu , subvol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
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
cylinder_DF( const int ND ,
	     const double q[ ND ] ,
	     const int DIMS ,
	     const double cyl_width )
{
  // test that this satisfies the correct cylinder cutting procedure
  const double norm = 1. / ( 2. * sqrt( DIMS ) ) ;
  const int diagonals = 2 << ( ND - 1 ) ;
 
  // generix for loop over the diagonals
  int mu ;
  double x[ ND ] ;
  for( mu = 0 ; mu < diagonals ; mu ++ ) {
    // inline for getting a diagonal for lexi order
    get_diagonal( ND , x , mu , DIMS ) ;

    double scalar_prod = 0. ; 
    int nu ;
    for( nu = 0 ; nu < ND ; nu ++ ) {
      scalar_prod += q[ nu ] * x[ nu ] ;
    }
    scalar_prod *= norm ;

    double mod = 0.0 ; 
    for( nu = 0 ; nu < ND ; nu ++ ) {
      register const double temp = q[ nu ] - scalar_prod * x[ nu ] ;
      mod += temp * temp ;
    }

    if( sqrt( mod ) <= cyl_width ) {
      return ADD_TO_LIST ;
    }
  }
  return DO_NOT_ADD ;
}

// set the cylinder according to momenta in the list
static void
set_cylinder( struct mom_info *momnew ,
	      const struct mom_info *momold ,
	      struct resampled *bootsnew ,
	      struct resampled *bootsold ,
	      double *xcyl ,
	      double *x ,
	      const int *list ,
	      const int num )
{
  int i ;
  for( i = 0 ; i < num ; i++ ) {
    bootsnew[ i ].resampled = (double*)malloc( bootsold[ list[ i ] ].NSAMPLES * sizeof( double ) ) ;
    memcpy( &momnew[ i ] , &momold[ list[ i ] ] , sizeof( struct mom_info ) ) ;
    memcpy( &bootsnew[ i ] , &bootsold[ list[ i ] ] , sizeof( struct resampled ) ) ;
    memcpy( &xcyl[ i ] , &x[ list[ i ] ] , sizeof( double ) ) ;
  }
}

// re-cylinder cut
static void
recylinder( struct input_params *INPARAMS ,    // change INPARAMS -> NDATA
	    struct resampled **boots ,
	    struct mom_info **mominfo ,    // pass by reference
	    double **x ,
	    const int NSLICES ,
	    const int LT ,
	    const int ND ,
	    const double width )
{
  int j ;

  // single pass to figure out what we can add
  int **list = malloc( NSLICES * sizeof( int* ) ) ;
  for( j = 0 ; j < NSLICES ; j++ ) {
    list[ j ] = (int*)malloc( INPARAMS -> NDATA[j] * sizeof( int ) ) ;
  }

  int num[ NSLICES ] ;
  for( j = 0 ; j < NSLICES ; j++ ) {
    double twiddles[ ND ] ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      twiddles[ mu ] = 2.0 * M_PI / (double)INPARAMS -> dimensions[j][mu] ; 
    }
    num[ j ] = 0 ;
    int i ;
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
      double q[ ND ] ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	q[ mu ] = mominfo[j][i].n[mu] * twiddles[ mu ] ;
      }
      // cut is a cylinder and 
      if( cylinder_DF( ND , q , ND , width ) == ADD_TO_LIST ) {
	list[ j ][ num[j] ] = i ;
	num[j]++ ;
      }
    }
  }

  printf( "After this %d vs. %d \n" , INPARAMS -> NDATA[ 0 ] , num[0] ) ;

#pragma omp parallel for private(j)
  for( j = 0 ; j < NSLICES ; j++ ) {

    // allocations
    struct mom_info *momcyl = malloc( num[j] * sizeof( struct mom_info ) ) ;
    struct resampled *cylinder = malloc( num[j] * sizeof( struct resampled ) ) ;
    double *xcyl = malloc( num[j] * sizeof( double ) ) ;

    // set the cylinder
    set_cylinder( momcyl , mominfo[j] ,
		  cylinder , boots[j] ,
		  xcyl , x[j] ,
		  list[j] , num[j] ) ;

    printf( "cylinder set .. \n" ) ;

    // reset the input parameters
    INPARAMS -> NDATA[j] = num[j] ;

    // and reallocate
    boots[ j ] = realloc( boots[ j ] , INPARAMS -> NDATA[j] * sizeof( struct resampled ) ) ;
    mominfo[ j ] = realloc( mominfo[ j ] , INPARAMS -> NDATA[j] * sizeof( struct mom_info ) ) ;
    x[ j ] = realloc( x[ j ] , INPARAMS -> NDATA[j] * sizeof( double ) ) ;

    // and copy over
    memcpy( boots[j] , cylinder , INPARAMS->NDATA[j] * sizeof( struct resampled ) ) ;
    memcpy( mominfo[j] , momcyl , INPARAMS->NDATA[j] * sizeof( struct mom_info ) ) ;
    memcpy( x[j] , xcyl , INPARAMS->NDATA[j] * sizeof( double ) ) ;

    // sort
    heapSort( x[j] , boots[j] , mominfo[j] , INPARAMS->NDATA[j] ) ;

    // and free
    free( xcyl ) ;
    free( cylinder ) ;
    free( momcyl ) ;
  }

  free( list ) ;

  return ;
}

// Alpha_s computation
void
adler_eval( double **xavg ,
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

  int j , i ;
  //recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , LT , 4 , 0.29 ) ;
  if( INPARAMS->dimensions[0][0] == 24 ) {
    recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , LT , 4 , 
		0.2*INPARAMS->quarks[0].ainverse ) ;
  } else {
    recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , LT , 4 , 
		0.15*INPARAMS->quarks[0].ainverse ) ;
  }

  for( i = 0 ; i < NSLICES ; i++ ) {
    xavg[ i ] = ( double* )malloc( INPARAMS -> NDATA[i] * sizeof( double ) ) ;
    int k ;
    for( k = 0 ; k < INPARAMS -> NDATA[ i ] ; k++ ) {
      xavg[i][k] = mominfo[i][k].p2 ;
    }
  }

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
      printf( "%f %f %f \n" , mavg[i][k].p2 , BAVG[i][k].avg , BAVG[i][k].err ) ;
    }
  }

  printf( "\n--> Converting to physical momenta <--\n" ) ;
  // OK, momentum first
  for( j = 0 ; j < NSLICES ; j++ ) {
    // ainverse multipliers
    //const double aI  = INPARAMS -> quarks[j].ainverse ;
    const double aI2 = pow( INPARAMS -> quarks[j].ainverse , 2 ) ;
    const double aI4 = pow( INPARAMS -> quarks[j].ainverse , 4 ) ;
    const double aI6 = pow( INPARAMS -> quarks[j].ainverse , 6 ) ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {

      // set to physical lattice spacing xavg[j][i] is now |p| in GeV!!
      x[j][i] = aI2 * x[j][i] ; //aI * sqrt( x[j][i] ) ;
      // and the moms
      mavg[j][i].p2 *= aI2 ;
      mavg[j][i].p4 *= aI4 ;
      mavg[j][i].p6 *= aI6 ;

      // check for some consistency
      if( fabs( mavg[j][i].p2 - x[j][i] ) > 1E-12 ) {
	printf( "P2 Broken !! %e \n" , mavg[j][i].p2 - x[j][i]*x[j][i] ) ;
	exit(1) ;
      } 
    }
  }
  
  // fit and plot are in here
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)BAVG , 
						    (const double**)x , 
						    (const struct mom_info **)mavg ,
						    *INPARAMS , NSLICES , LT ) ;

  fitfunc fit ;
  int NPARAMS , NCOMMON ;
  initialise_f( INPARAMS->fittype , &fit , &NPARAMS , &NCOMMON ) ;

  // unpack the simultaneous fit into something more palatable
  struct resampled *fit1 = malloc( NPARAMS * sizeof( struct resampled ) ) ;

  for( i = 0 ; i < NPARAMS ; i++ ) {
    fit1[i].resampled = (double*)malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
  }

  // open a derivative graph
  char tmp[256] ;
  sprintf( tmp , "der_%s" , INPARAMS->graph_name ) ;
  make_xmgrace_graph( tmp ,
		      INPARAMS->graph_xaxis , 
		      INPARAMS->graph_yaxis ) ;

  // compute the derivative
  for( i = 0 ; i < NSLICES ; i++ ) {
    int k ;
    if( i > 0 ) {
      int check = i*( NPARAMS - NCOMMON ) + NCOMMON ;
      for( k = 0 ; k < NPARAMS ; k++ ) {
	if( INPARAMS -> sim_params[ k ] == true ) {
	  equate( &fit1[ k ] , fitparams[ k ] ) ;
	} else {
	  equate( &fit1[ k ] , fitparams[ check ] ) ;
	  check++ ;
	}
      }
    } else {
      #pragma omp parallel for private(k)
      for( k = 0 ; k < NPARAMS ; k++ ) {
	equate( &fit1[ k ] , fitparams[ k ] ) ;
      }
    }

    for( k = 0 ; k < NPARAMS ; k++ ) {
      printf( "%f \n" , fit1[k].avg ) ;
    }

    int range = 500 ;
    struct resampled *y = malloc( range * sizeof( struct resampled ) ) ;
    double *xx = malloc( range * sizeof( double ) ) ;
    for( j = 0 ; j < range ; j++ ) {
      xx[j]  = INPARAMS->fit_lo + j * ( INPARAMS->fit_hi - INPARAMS->fit_lo ) / range ;
      y[j]   = fit_der( fit1 , xx[j] , NPARAMS , INPARAMS -> quarks[i] , 
			mavg[i][j] , fit , LT , NPARAMS ) ;
      mult_constant( &y[j] , xx[j] * 12.0 * M_PI * M_PI * 5.0 / 9.0 ) ;
    }

    plot_data( y , xx , range ) ;

    free( y ) ;
  }
    
  // close up the graph
  close_xmgrace_graph(  ) ;

  free( fitparams ) ;

  return ;
}
