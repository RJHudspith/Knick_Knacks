#include "fitfunc.h"
#include "rng.h" // for the bootstraps

// either we boot average or use the ensemble average ...
#define BOOT_AVERAGE

static double *boot_order ;

void
init_boot_order( const int NBOOTS , 
		 const int NRAW )
{
  rng_reseed( ) ;
  boot_order = (double*)malloc( NBOOTS * NRAW * sizeof( double ) ) ;
  int i ;
  for( i = 0 ; i < NBOOTS * NRAW  ; i++ ) {
    boot_order[ i ] = rng_double( ) ;
  }
  return ;
}

void
free_boot_order( void )
{
  free( boot_order ) ;
  return ;
}

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

// and fill up the raw resampled data ...
void
raw_err( struct resampled *replicas )
{
  double ave , err ;
  average( &ave , &err , replicas -> resampled , replicas -> NSAMPLES ) ;
  // use the standard error 
  err = sqrt( err / ( (double)replicas -> NSAMPLES * ( replicas -> NSAMPLES - 1.0 ) ) ) ;
  replicas -> avg = ave ;
  replicas -> err_hi = ave + err ;
  replicas -> err_lo = ave - err ;
  replicas -> err    = err ;
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

// qsort comparison function for the bootstrap
int 
comp( const void *elem1 , 
      const void *elem2 ) 
{
  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f > s) return  1 ;
  if (f < s) return -1 ;
  return 0 ;
}

// compute the error for the bootstrap assumes the data has
// been bootstrapped ...
void
bootstrap_error( struct resampled *replicas )
{
  double sorted[ replicas -> NSAMPLES ] ;
  memcpy( sorted , replicas -> resampled , sizeof( double ) * ( replicas -> NSAMPLES ) ) ;

  // sort the bootstraps
  qsort( sorted , replicas -> NSAMPLES , sizeof( double ) , comp ) ;

  // average is better behaved when sorted ...
  double ave = 0.0 ;
  int i ;
  for( i = 0 ; i < ( replicas -> NSAMPLES ) ; i++ ) {
    ave += sorted[ i ] ;
  }
  ave /= (double)( replicas -> NSAMPLES ) ;

  // confidence bounds are at 1 sigma
  const double confidence = 68.2689492 ;
  const double omitted = 0.5 * ( 100. - confidence ) ;
  const int bottom = (int)( ( omitted * replicas -> NSAMPLES ) / 100. ) ;
  const int top = ( ( replicas -> NSAMPLES ) - 1 - bottom ) ;

#ifdef BOOT_AVERAGE
  replicas -> avg = ave ;
#endif
  replicas -> err_hi = sorted[ top ] ;
  replicas -> err_lo = sorted[ bottom ] ;
  // symmetrized error
  replicas -> err    = 0.5 * ( replicas -> err_hi - replicas -> err_lo ) ;
  return ;
}

// and the bootstrap analysis , assumes we have NBOOTS of space in Bootstrap
void
bootstrap( struct resampled *Bootstrap ,
	   const struct resampled Raw )
{
  // reseed the rng ...
  rng_reseed( ) ;
  // loop the bootstrap index
  int i ;
  for( i = 0 ; i < Bootstrap -> NSAMPLES ; i++ ) {
    Bootstrap -> resampled[ i ] = 0.0 ;    
    // loop the data
    int j ;
    for( j = 0 ; j < Raw.NSAMPLES ; j++ ) {
      Bootstrap -> resampled[ i ] += 
	Raw.resampled[ (int)( boot_order[ j + Raw.NSAMPLES * i ] * ( Raw.NSAMPLES - 1 ) ) ] ;
    }
    Bootstrap -> resampled[ i ] /= (double)Raw.NSAMPLES ;
  }
#ifndef BOOT_AVERAGE
  Bootstrap -> avg = Raw.avg ;
#endif
  bootstrap_error( Bootstrap ) ;
  return ;
}

// a wrapper for the resampling
void
compute_err( struct resampled *replicas )
{
  switch( replicas -> restype ) {
  case RAWDATA :
    raw_err( replicas ) ;
    break ;
  case JACKDATA :
    jackknife_error( replicas ) ;
    break ;
  case BOOTDATA :
    bootstrap_error( replicas ) ;
    break ;
  default :
    printf( "<compute_err> Unrecognised resampling %d!! \n" , 
	    replicas -> restype ) ;
    exit(1) ;
    break ;
  }
  return ;
}

// bootstrap or jackknife or whatever
struct resampled *
resample_data( const struct resampled *RAW ,
	       const int N ,
	       const resample_type restype ,
	       const int NBOOTS )
{
  struct resampled *BOOT = malloc( N * sizeof( struct resampled ) ) ;
  const int NRAW = RAW[0].NSAMPLES ;
  int i ;

  printf( "[STATS] NBOOTS :: %d \n" , NBOOTS ) ;

  // and resample could be done in parallel because we have the rng
  // order as a look up table
  #pragma omp parallel for private(i)
  for( i = 0 ; i < N ; i++ ) {
    // resampling
    int k ;

    if( RAW[i].restype != RAWDATA ) {

      BOOT[i].resampled = (double*)malloc( RAW[i].NSAMPLES * sizeof( double ) ) ;
      equate( &BOOT[i] , RAW[i] ) ;

    } else {

      BOOT[i].restype = restype ;

      switch( restype ) {
      case RAWDATA :
	BOOT[i].resampled = (double*)malloc( NRAW * sizeof( double ) ) ;
	BOOT[i].NSAMPLES  = NRAW ;
	for( k = 0 ; k < NRAW ; k++ ) {
	  BOOT[ i ].resampled[ k ] = RAW[ i ].resampled[ k ] ;
	}
	raw_err( &BOOT[i] ) ;
	break ;
      case BOOTDATA :
	BOOT[i].resampled = (double*)malloc( NBOOTS * sizeof( double ) ) ;
	BOOT[i].NSAMPLES  = NBOOTS ;
	bootstrap( &BOOT[i] , RAW[i] ) ;
	break ;
      case JACKDATA :
	// I manipulate NBOOTS to be NRAW-1 in the input file
	BOOT[i].resampled = (double*)malloc( NBOOTS * sizeof( double ) ) ;
	BOOT[i].NSAMPLES  = NBOOTS ;
	jackknife( &BOOT[i] , RAW[i] ) ;
	break ;
      default :
	printf( "<resample_data> Unrecognised resampling %d !! \n" , restype ) ;
	exit(1) ;
	break ;
      }
      // break for the check on if it is rawdata or not
    }
  }
  return BOOT ;
}

// wraps the above into an array
struct resampled **
resample_array( const struct resampled **RAW ,
		const int NSLICES ,
		const int *NDATA ,
		const int resample , 
		const int NBOOTS )
{
  struct resampled **BOOT = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  int i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    BOOT[i] = resample_data( RAW[i] , 
			     NDATA[i] , 
			     resample , 
			     NBOOTS ) ;
  }
  return BOOT ;
}

// bins the raw data
struct resampled
bin_the_data( const struct resampled RAW ,
	      const int binning )
{
  const int NBINNED = (int)RAW.NSAMPLES / (int)binning ;
  struct resampled BINNED ;

  BINNED.resampled = malloc( NBINNED * sizeof( double ) ) ;
  BINNED.NSAMPLES = NBINNED ;
  BINNED.restype  = RAWDATA ;

  // OK, so we simply average within each bin
  int j ;
  for( j = 0 ; j < NBINNED ; j++ ) {
    
    // loop the elemnts within the bin, combining for an average
    double bin_ave = 0.0 ;
    int k ;
    for( k = 0 ; k < binning ; k++ ) {
      bin_ave += RAW.resampled[ k + j*binning ] ;
    }
    BINNED.resampled[j] = bin_ave / (double)binning ;
  }
  // and that is it

  return BINNED ;
}


// extrapolates a fit to the point "x"
struct resampled
extrapolate( const struct resampled *fitparams ,
	     const double extrap_point ,
	     const int NPARAMS ,
	     const struct chiral quarks ,
	     const struct mom_info mominfo ,
	     const fitfunc fit ,
	     const int LT ,
	     const int SLICE )
{
  struct resampled y ;
  y.resampled = malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
  y.NSAMPLES  = fitparams[0].NSAMPLES ;
  y.restype   = fitparams[0].restype ;

  // for each x, evaluate the fit function on average
  int j ; 
  #pragma omp parallel for private(j)
  for( j = 0 ; j < y.NSAMPLES ; j++ ) {
    double data[ NPARAMS ] ;
    // switch to individual axis
    int k ;
    for( k = 0 ; k < NPARAMS ; k++ ) {
      data[ k ] = fitparams[ k ].resampled[ j ] ;
    }
    
    // CCCCCC-Callback for the fitfunction
    struct x_descriptor XX = { extrap_point , quarks , mominfo , LT } ;
    y.resampled[ j ] = fit.f( data , XX , SLICE ) ;
  }

  // do the average?
  double data[ NPARAMS ] ;
  // switch to individual axis
  int k ;
  for( k = 0 ; k < NPARAMS ; k++ ) {
    data[ k ] = fitparams[ k ].avg ;
  }
    
  // CCCCCC-Callback for the fitfunction
  struct x_descriptor XX = { extrap_point , quarks , mominfo , LT } ;
  y.avg = fit.f( data , XX , SLICE ) ;

  // this has a switch for the resample type
  compute_err( &y ) ;

  return y ;
}
