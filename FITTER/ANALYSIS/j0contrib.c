/**
   Analysis for the VPF
 */
#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "shuffalgo.h"
#include "equivalents.h"
#include "Utils.h"
#include "correlation.h"

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
	    const double *widths )
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
      if( cylinder_DF( ND , q , ND , widths[j] ) == ADD_TO_LIST ) {
	list[ j ][ num[j] ] = i ;
	num[j]++ ;
      }
    }
    printf( "After this %d -> %d \n" , INPARAMS -> NDATA[ j ] , num[j] ) ;
  }

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

    printf( "cylinder %d set .. \n" , j ) ;

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
j0_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ,
	     const bool renormalise )
{
  // projs
  const char *proj[ 3 ] = { "Long" , "Trans" , "TransPLong" } ;

  size_t j , i ;

  printf( "\n--> J0 evaluation <--\n" ) ;

  printf( "\n--> Recylinder <--\n" ) ;

  const double widths[3] = { 0.32 , 0.32 , 0.32 } ;
  recylinder( INPARAMS , bootavg , mominfo , xavg , NSLICES , LT , 4 , widths ) ;

  // BOOT and X are freed in momavg
  struct mom_info **mavg ;
  struct resampled **bavg = momentum_average( &mavg , INPARAMS , 
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
  // OK, momentum first
  for( j = 0 ; j < NSLICES ; j++ ) {
    // ainverse multipliers
    const double aI  = INPARAMS -> quarks[j].ainverse ;
    const double aI2 = pow( INPARAMS -> quarks[j].ainverse , 2 ) ;
    const double aI4 = pow( INPARAMS -> quarks[j].ainverse , 4 ) ;
    const double aI6 = pow( INPARAMS -> quarks[j].ainverse , 6 ) ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {

      // set to physical lattice spacing xavg[j][i] is now |p| in GeV!!
      x[j][i] = aI * sqrt( x[j][i] ) ;
      // and the moms
      mavg[j][i].p2 *= aI2 ;
      mavg[j][i].p4 *= aI4 ;
      mavg[j][i].p6 *= aI6 ;

      // check for some consistency
      if( fabs( mavg[j][i].p2 - x[j][i]*x[j][i] ) > 1E-12 ) {
	printf( "P2 Broken !! %e \n" , mavg[j][i].p2 - x[j][i]*x[j][i] ) ;
	exit(1) ;
      } 
    }
  }

  // multiply by Z_V
  printf( "\n--> Renormalising <--\n" ) ;

  // output the data
  {
    // renormalise 
    for( j = 0 ; j < NSLICES ; j++ ) {

      char str[ 256 ] ;
      sprintf( str , "m%g_m%g.%s.dat" , INPARAMS -> quarks[j].ml , INPARAMS -> quarks[j].ms , proj[ j ] ) ;
      FILE *file = fopen( str , "w" ) ;

      // renormalise
      const double ZV = INPARAMS -> quarks[j].ZV ;
      fprintf( file , "ZV :: %f\n" , ZV ) ;
      fprintf( file , "a^{-1} :: %f GeV\n" , INPARAMS -> quarks[j].ainverse ) ;
      /*
#pragma omp parallel for private(i)
      for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
	mult_constant( &bavg[j][i] , ZV ) ;
      }
      */
      fprintf( file , "q [GeV]        Pi(q^2)             Err\n" ) ;
      // and print it out
      for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {
	fprintf( file , "%f %1.12E %1.12E \n" , x[j][i] , bavg[j][i].avg , bavg[j][i].err ) ;
      }
      fclose( file ) ;
    }
  }

  // compute correlations
  {
    for( i = 0 ; i < NSLICES ; i++ ) {

      char str[ 256 ] ;
      sprintf( str , "m%g_m%g.%s.cov" , INPARAMS -> quarks[i].ml , INPARAMS -> quarks[i].ms , proj[ i ] ) ;

      FILE *file = fopen( str , "w" ) ;

      double **correlation = malloc( INPARAMS -> NDATA[i] * sizeof( double* ) ) ;
      for( j = 0 ; j < INPARAMS -> NDATA[i] ; j++ ) {
	correlation[j] = malloc( INPARAMS -> NDATA[i] * sizeof( double ) ) ;
      }
      correlations( correlation , bavg[i] , INPARAMS -> NDATA[i] ) ;

      write_corrmatrix_to_file( file , correlation , INPARAMS -> NDATA[i] ) ;

      free_double_dp( correlation , INPARAMS -> NDATA[i] ) ;

      fclose( file ) ;
    }
  }

  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bavg , 
						    (const double**)x , 
						    (const struct mom_info **)mavg ,
						    *INPARAMS , NSLICES , LT ) ;
  
  free( fitparams ) ;

  free_resampled_dp( bavg , NSLICES , INPARAMS->NDATA ) ;
  free_double_dp( x , NSLICES ) ;

  return ;
}
