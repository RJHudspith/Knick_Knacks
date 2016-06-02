/*
  I would like to put my merge sort here as it is lovely
  but alas the heapsort will have to do
 */

#include "fitfunc.h"
#include "stats.h"

static void
swap( double *swapped , const int idx1 , const int idx2 )
{
  const double temp = swapped[ idx1 ] ;
  swapped[ idx1 ] = swapped[ idx2 ] ;
  swapped[ idx2 ] = temp ;
  return ;
}

// need this for the y-data
static void
swap_resampled( struct resampled *swapped , const int idx1 , const int idx2 )
{
  struct resampled temp ;
  temp.resampled = malloc( swapped[0].NSAMPLES * sizeof( double ) ) ;
  memcpy( &temp , &swapped[idx1] , sizeof( struct resampled ) ) ;
  memcpy( &swapped[idx1] , &swapped[idx2] , sizeof( struct resampled ) ) ;
  memcpy( &swapped[idx2] , &temp , sizeof( struct resampled ) ) ;
  //free( temp.resampled ) ;
  return ;
}

static void
swap_mominfo( struct mom_info *swapped , const int idx1 , const int idx2 )
{
  struct mom_info temp ;
  memcpy( &temp , &swapped[idx1] , sizeof( struct mom_info ) ) ;
  memcpy( &swapped[idx1] , &swapped[idx2] , sizeof( struct mom_info ) ) ;
  memcpy( &swapped[idx2] , &temp , sizeof( struct mom_info ) ) ;
  return ;
}

/// heapsort
static void 
siftDown( double *x , 
	  struct resampled *y , 
	  struct mom_info *mominfo ,
	  int root , int bottom )
{
  int done = 0 , maxChild ; 

  while ( ( root*2 <= bottom ) && ( !done ) ) {
    if( root*2 == bottom ) {
      maxChild = root * 2 ; 
    } else if ( x[root * 2] > x[root * 2 + 1] ) {
      maxChild = root * 2 ; 
    } else {
      maxChild = root * 2 + 1 ; 
    }
    
    if ( x[root] < x[maxChild] ) {
      // do a swap on the x and y again
      swap( x , root , maxChild ) ;
      swap_resampled( y , root , maxChild ) ;
      swap_mominfo( mominfo , root , maxChild ) ;
      root = maxChild ; 
    } else {
      done = 1 ;
    } 
  }
  return ;
}

// heapsorts the data wrt to x index. y is along for the ride
void 
heapSort( double *x , 
	  struct resampled *y ,
	  struct mom_info *mominfo ,
	  const int array_size )
{
  int i ; 
  for ( i = ( array_size / 2 ) ;  i >= 0 ;  i-- ) {
    siftDown( x , y , mominfo ,  i ,  array_size - 1 ) ; 
  }
  
  // swapsies
  for ( i = array_size-1 ;  i >= 1 ;  i-- ) {
    // swap the x and the y
    swap( x , 0 , i ) ;
    swap_resampled( y , 0 , i ) ;
    swap_mominfo( mominfo , 0 , i ) ;
    // and go down the rabbit hole
    siftDown( x , y , mominfo , 0 ,  i-1 ) ; 
  }
  return ;
} 
