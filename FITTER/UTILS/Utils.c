/**
   Some general utilities for doing things erroneous code.
 */

#include "fitfunc.h"

/**
   divide and conquer finds closest X to the value target
   assumes X is a sorted array
   maybe one day I will comment this code
*/
int 
find_idx( const double target , const double *X , 
	  const int high , const int low )
{
  if( target > 0 ) {
    return ( ( high - low ) < 2 ) ?					\
      fabs( X[high] - target ) < fabs( X[low] - target ) ? high : low	\
      :									\
      X[ ( high + low ) >> 1 ] < target ? find_idx( target , X , high , ( high + low ) >> 1 ) : \
      find_idx( target , X , ( high + low ) >> 1 , low ) ;
  } else {
    size_t i ;
    int idx = low ;
    double best = fabs( target ) ;
    for( i = low ; i < high ; i++ ) {
      if( fabs( X[i] - target ) < best ) {
	best = fabs( X[i] - target ) ;
	idx = i ;
      }
    }
    return idx ;
  }
}

/**
   set up the params from the gsl_vector
   stops us from calling gsl_vector_get each time we want one of the
   params, reduces associated overhead
 */
double *
set_params( const gsl_vector *x , 
	    const int NPARAMS )
{
  double *params = malloc( NPARAMS * sizeof( double ) ) ;
  int i ;
  for( i = 0 ; i < NPARAMS ; i++ ) {
    params[ i ] = gsl_vector_get( x , i ) ;
  }
  return params ;
}

// free some resampled data
void
free_resampled( struct resampled *data , const int NDATA ) 
{
  int j ;
  for( j = 0 ; j < NDATA ; j++ ) {
    free( data[j].resampled ) ;
  }
  free( data ) ;
}

// free a double pointer of resampled data
void
free_resampled_dp( struct resampled **data , const int idx1 , const int *idx2 )
{
  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < idx1 ; i++ ) {
    free_resampled( data[i] , idx2[i] ) ;
  }
  free( data ) ;
  return ;
}

// free a double double pointer
void
free_double_dp( double **data , const int idx1 )
{
  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < idx1 ; i++ ) {
    free( data[i] ) ;
  }
  free( data ) ;
  return ;
}

// free the mominfo struct
void
free_mominfo( struct mom_info **data , const int idx1 )
{
  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < idx1 ; i++ ) {
    free( data[i] ) ;
  }
  free( data ) ;
  return ;
}

// bootstrap output
void
printboots( const struct resampled R )
{
  printf( "RESAMPLED \n" ) ;
  int i ;
  for( i = 0 ; i < R.NSAMPLES ; i++ ) {
    printf( "%d %f \n" , i , R.resampled[i] ) ;
  }
  printf( "\n" ) ;
  return ;
}

// flatten a resampled array R[NSLICES][NDATA] to an array RES[0][NSLICES*NDATA]
struct resampled*
flatten_resampled_array( const struct resampled **R ,
			 const int idx1 ,
			 const int *idx2 )
{
  ///
  int j , sum = 0 ;
  for( j = 0 ; j < idx1 ; j++ ) {
    sum += idx2[j] ;
  }

  int idx = 0 ;
  /// create an end-to-end array for the simultaneous fit
  struct resampled *RES = malloc( sum * sizeof( struct resampled ) ) ;

  for( j = 0 ; j < idx1 ; j++ ) {
    int k ;
    for( k = 0 ; k < idx2[j] ; k++ ) {
      // make space for it
      RES[idx].resampled = malloc( R[j][k].NSAMPLES * sizeof( double ) ) ;
      equate( &RES[idx++] , R[j][k] ) ;
    }
  }
  return RES ;
}

// flatten a double array R[NSLICES][NDATA] to an array RES[0][NSLICES*NDATA]
double*
flatten_double_array( const double **R ,
		      const int idx1 ,
		      const int *idx2 )
{
  // compute the length of him
  int j , sum = 0 ;
  for( j = 0 ; j < idx1 ; j++ ) {
    sum += idx2[j] ;
  }

  /// create an end-to-end array for the simultaneous fit
  double *RES = malloc( sum * sizeof( double ) ) ;

  int idx = 0 ; 
  //#pragma omp parallel for private(j)
  for( j = 0 ; j < idx1 ; j++ ) {
    int k ;
    for( k = 0 ; k < idx2[j] ; k++ ) {
      RES[idx++] = R[j][k] ;
    }
  }
  return RES ;
}

// flatten a mom array R[NSLICES][NDATA] to an array RES[0][NSLICES*NDATA]
struct mom_info*
flatten_mom_array( const struct mom_info **R ,
		   const int idx1 ,
		   const int *idx2 )
{
  // compute the length of him
  int j , sum = 0 ;
  for( j = 0 ; j < idx1 ; j++ ) {
    sum += idx2[j] ;
  }

  /// create an end-to-end array for the simultaneous fit
  struct mom_info *RES = malloc( sum * sizeof( struct mom_info ) ) ;

  int idx = 0 ; 
  for( j = 0 ; j < idx1 ; j++ ) {
    int k ;
    for( k = 0 ; k < idx2[j] ; k++ ) {
      memcpy( &RES[idx++] , &R[j][k] , sizeof( struct mom_info ) ) ;
    }
  }
  return RES ;
}

// create a double array from the resampled data built of some specific aspect
double*
resampled_to_double( const struct resampled *R ,
		     const int size ,
		     const int which )
{
  double *RES = malloc( size * sizeof( double* ) ) ;
  int j ;
#pragma omp parallel for private(j)
  for( j = 0 ; j < size ; j++ ) {
    switch( which ) {
      // for the sigma I use the distribution width
      // unless something gets really funky then I use the
      // average?
      case ERR : 
	RES[j] = R[j].err ;
	RES[j] = ( RES[j] <= 1E-20 ) ? 1E-10 : RES[j] ;
	break ;
    case HI :  RES[j] = R[j].err_hi ; break ;
    case LO :  RES[j] = R[j].err_lo ; break ;
    default :  RES[j] = R[j].avg ; break ;
    }
  }
  return RES ;
}


