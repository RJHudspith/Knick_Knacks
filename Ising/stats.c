#include "Ising.h"

// compute the average and returns the unnormalised variance
void
average( double *ave , double *err , 
	 const double *meas , const int N )
{
  int i ;
  *err = *ave = 0.0 ;
  for( i = 0 ; i < N ; i++ ) {
    const double delta = meas[i] - *ave ;
    *ave += delta / ((double)i+1.0) ;
    *err += delta * ( meas[i] - *ave ) ; 
  }
  return ;
}

// the jackknife error definition incorporates the code for the average and variance
void
jackknife_error( struct resampled *replicas )
{
  // average and variance come from the usual defs
  double ave , err ;
  average( &ave , &err , replicas -> resampled , replicas -> NSAMPLES ) ;

  err = sqrt( err * ( 1.0 - 1.0 / (double)replicas -> NSAMPLES ) ) ;
  replicas -> avg = ave ;
  replicas -> err_hi = ave + err ;
  replicas -> err_lo = ave - err ;
  replicas -> err    = err ;
  return ;
}

// perform a jackknife on the Raw data
// assumes Jackknife has had NRAW - 1 allocation
void
jackknife( struct resampled *Jackknife ,
	   const struct resampled Raw )
{
  const int N = Raw.NSAMPLES ;
  Jackknife -> NSAMPLES = N ; // should be set anyway ..
  int i ;

  double sum = 0.0 ;
  for( i = 0 ; i < N ; i++ ) {
    sum += Raw.resampled[i] ; // compute the sum
  }
  // subtract each element from the sum
  const double NORM = 1.0 / ( (double)N - 1.0 ) ;
  for( i = 0 ; i < N ; i++ ) {
    Jackknife -> resampled[ i ] = ( sum - Raw.resampled[i] ) * NORM ;
  }
  jackknife_error( Jackknife ) ;
  return ;
}

static void
NSAMPLES_ERR( const char *OP , const int a , const int b ) 
{
  printf( "[%s]\n" , OP ) ;
  printf( "the two distributions have different number of SAMPLES !!" ) ;
  printf( "%d vs %d \n" , a , b ) ;
  printf( "Leaving in disgust\n" ) ;
  exit(1) ;
}

static void
RESAMPLE_ERR( const char *OP , const int a , const int b ) 
{
  printf( "[%s]\n" , OP ) ;
  printf( "the two distributions have different sample types !!" ) ;
  printf( "%d vs %d \n" , a , b ) ;
  printf( "Leaving in disgust\n" ) ;
  exit(1) ;
}

// is always called
static void
resampled_checks( const char *OP , 
		  const struct resampled a , 
		  const struct resampled b ) 
{
  if( a.NSAMPLES != b.NSAMPLES ) { 
    NSAMPLES_ERR( OP , a.NSAMPLES , b.NSAMPLES ) ; 
  } else if( a.restype != b.restype ) {
    RESAMPLE_ERR( OP , a.restype , b.restype ) ; 
  }
  // and any other checks we might wish to perform ....
  return ;
}

// atomic, a -= b for the distribution
void
subtract( struct resampled *a , 
	  const struct resampled b ) 
{
  resampled_checks( "subtract" , *a , b ) ;
  int i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] -= b.resampled[i] ;
  }
  jackknife_error( a ) ;
  return ;
}

//
void
equate( struct resampled *a , 
	const struct resampled b )
{
  int l ;
  for( l = 0 ; l < b.NSAMPLES ; l++ ) {
    a -> resampled[l] = b.resampled[l] ;
  }
  a -> NSAMPLES = b.NSAMPLES ;
  a -> restype  = b.restype ;
  a -> avg = b.avg ;
  a -> err = b.err ;
  a -> err_hi = b.err_hi ;
  a -> err_lo = b.err_lo ;
  return ;
}

//
void
init_stats( struct resampled *raw , struct resampled *rawsq ,
	    struct resampled *jack , struct resampled *jacksq ,
	    struct resampled *susc , const int MEASUREMENTS )
{
  const int NBOOTS = MEASUREMENTS - 1 ;

  // array of measurements
  raw->resampled = (double*)malloc( MEASUREMENTS * sizeof( double ) ) ;
  raw->NSAMPLES = MEASUREMENTS ;
  raw->restype = RAWDATA ;

  rawsq->resampled = (double*)malloc( MEASUREMENTS * sizeof( double ) ) ;
  rawsq->NSAMPLES = MEASUREMENTS ;
  rawsq->restype = RAWDATA ;

  // array of measurements
  jack->resampled = (double*)malloc( MEASUREMENTS * sizeof( double ) ) ;
  jack->NSAMPLES = NBOOTS ;
  jack->restype = JACKDATA ;

  jacksq->resampled = (double*)malloc( MEASUREMENTS * sizeof( double ) ) ;
  jacksq->NSAMPLES = NBOOTS ;
  jacksq->restype = JACKDATA ;

  susc->resampled = (double*)malloc( MEASUREMENTS * sizeof( double ) ) ;
  susc->NSAMPLES = NBOOTS ;
  susc->restype = JACKDATA ;

  return ;
}
