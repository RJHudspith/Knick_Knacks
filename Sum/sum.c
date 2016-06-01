/**
  @file sum.c
  @brief sum a large array using methods with varying levels of success
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

static const double d = 4.656612873077393e-10 ;

static struct timeval GLUtimer ;
static double t1 ;

// print to stdout the time elapsed in seconds, minutes or hours
double
print_time( void )
{
  gettimeofday( &GLUtimer , NULL ) ;
  double diff = GLUtimer.tv_sec + ( GLUtimer.tv_usec/ 1E6 ) - t1 ;
  printf( "\n[TIMER] elapsed :: " ) ;
  if( diff > 60 ) {
    diff /= 60. ;
    if( diff > 60 ) { // cocky
      diff /= 60. ;
      printf( "%f (hours) \n", diff ) ;
    } else {
      printf( "%f (minutes) \n", diff ) ;
    }
  } else {
    printf( "%f (seconds) \n", diff ) ;
  }
  return diff ;
}

// start the timer
void
start_timer( void )
{
  gettimeofday( &GLUtimer , NULL ) ;
  t1 = GLUtimer.tv_sec + ( GLUtimer.tv_usec / 1E6 ) ;
  return ;
}

static void
build_array( float *a ,
	     const size_t N )
{
  size_t i ;
  for( i = 0 ; i < N ; i++ ) {
    a[ i ] = rand() * (float)d ;
  }
  return ;
}

static float
naive_sum( float *a ,
	   const size_t N )
{
  register float sum = 0.0 ;
  size_t i ;
  for( i = 0 ; i < N ; i++ ) {
    sum += a[i] ;
  }
  return sum ;
}

static double
double_sum( float *a ,
	   const size_t N )
{
  register double sum = 0.0 ;
  size_t i ;
  for( i = 0 ; i < N ; i++ ) {
    sum += (double)a[i] ;
  }
  return sum ;
}

// divide and conquer sum
static float
conker( float *a ,
	const size_t lo ,
	const size_t hi )
{
  return ( hi - lo ) < 2 ? a[lo] : \
    conker( a , lo , (hi+lo)/2 ) + \
    conker( a , (hi+lo)/2 , hi ) ;
}

// Kahan summation
static float
kahan( float *a ,
       const size_t N )
{
  size_t i ;
  register float sum = 0.0 , c = 0.0 , y , t ;
  for( i = 0 ; i < N ; i++ ) {
    y = a[i] - c ; // subtract correction
    t = sum + y ;
    c = ( t - sum ) - y ;
    sum = t ;
  }
  return sum ;
}

// parallel divide and conquer
static float
par_conker( float *a ,
	    const size_t N )
{
  size_t i , Nth = 1 ;
  #pragma omp parallel
  {
    Nth = omp_get_num_threads() ;
  }
  const size_t split = ( N + Nth - 1 ) / Nth ;
  register float sum = 0 ;
  #pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < Nth ; i++ ) {
    sum = sum + (float)conker( a , i * split , (i+1)*split>N?N:(i+1)*split ) ;
  }
  return sum ;
}

// run the tests
int main( void )
{
  srand( 1234 ) ;
  
  const size_t N = 4E8 ;

  float *a = malloc( N * sizeof( float ) ) ;

  build_array( a , N ) ;

  start_timer() ;

  //printf( "Naive :: %1.9f \n" , naive_sum( a , N )/N ) ;
  //printf( "Kahan  :: %1.9f \n" , kahan( a , N )/N ) ;
  //printf( "Conker :: %1.9f \n" , conker( a , 0 , N )/N ) ;
  printf( "ParCon :: %1.9f \n" , par_conker( a , N )/N ) ;
  //printf( "Double :: %1.9f \n" , double_sum( a , N )/N ) ;

  print_time() ;

  free( a ) ;

  return 0 ;
}
