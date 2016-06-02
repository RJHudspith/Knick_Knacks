/**
   @file cylinder.c
   @brief cylinder cuts and what have you
 */
#include "fitfunc.h"
#include "equivalents.h"
#include "shuffalgo.h"

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
void
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
