/*
  2D Ising model code
 */
#include "Ising.h"

#include "geometry.h"
#include "RNG.h"
#include "stats.h"

// static arrays
struct general Latt ;

// fronteir memory are the size of the
// lattice because memory is cheap
int *newfront , *oldfront ;

// wolff cluster algorithm
//
// randomly select a site,
// for each neighbour with the same spin, build a fronteir if
// a random number is less than the clustering probability
static int
wolff( struct site *__restrict lat ,
       int *__restrict mag ,
       const double p )
{
  // some initial stuff
  int front_size = 1 , allfront = 1 , front_old = 1 ;
  int mu , i ;

  // have a linked list of the full fronteir
  oldfront[0] = ( (int)( KISS_dbl( ) * LVOLUME ) ) ;

  const int ref_spin = lat[ oldfront[0] ].spin ;
  const int flipped_value = -ref_spin ;
  register int thispos ;

  // flip the first one
  lat[ oldfront[0] ].spin = flipped_value ;

  register int prop ;
  // loop until we have no similar spins
  while( front_size > 0 ) {
    front_size = 0 ;
    // loop the old fronteir ...
    for( i = 0 ; i < front_old ; i++ ) {
      // stop rereading this
      thispos = oldfront[i] ;
      // loop directions?
      for( mu = 0 ; mu < 2*ND ; mu++ ) {
	// proposed direction
	prop = lat[ thispos ].neighbor[ mu ] ;
	if( lat[ prop ].spin == ref_spin ) {
	  if( p < KISS_dbl( ) ) {
	    // spin flip in here
	    newfront[ front_size ] = prop ;
	    lat[ prop ].spin = flipped_value ;
	    front_size++ ;
	  }
	}
      }
      // end of loop directions
    }
    // set the old one to the new front
    for( i = 0 ; i < front_size ; i++ ) {
      oldfront[ i ] = newfront[ i ] ;
    }
    // accumulate the full size of the cluster
    allfront += front_size ;
    // equate the old fronteir size to the new one
    front_old = front_size ;
  }
  // update mag
  *mag -= ( ref_spin * 2 * allfront ) ;

  return allfront ;
}

// perform a hot start
int
hot_start( struct site *lat )
{
  int i , MAG = 0 ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    lat[i].spin = ( KISS_dbl( ) < 0.5 ) ? -1 : 1 ;
    MAG += lat[i].spin ;
  }
  return MAG ;
}

// perform a cold start
int
cold_start( struct site *lat )
{
  int i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    lat[i].spin = 1 ;
  }
  return LVOLUME ;
}

// initialise lattice information
static void
init_latt( const int dims[ ND ] , const int maxiters ,
	   const int stride , const int thermalisation )
{
  // dimensions
  int mu ;
  // compute the volume
  Latt.vol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[mu] = dims[mu] ;
    Latt.vol *= Latt.dims[mu] ;
  }
  Latt.MAX_ITERS = maxiters ; //10000 ;
  Latt.NMEAS = stride ; //100 ;
  Latt.THERM = thermalisation ; //2000 ;
  return;
}

// every good program needs to tell you how to use it
static int
usage( void ) 
{
  return printf( "USAGE :: ./ISING L_0 L_1 ... L_N"
		 "(compiled for ND=%d)\n"
		 , ND ) ;
}

// main
int 
main( const int argc , char *argv[] ) 
{
  // check
  if( argc != ND+1 ) {
    return usage( ) ;
  }

  // dimensions are read from argv
  int DIMS[ ND ] , mu ;
  printf( "[DIMS] :: " ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    DIMS[ mu ] = (int)atoi( argv[ 1 + mu ] ) ;
    printf( " %d " , DIMS[ mu ] ) ;
  }
  printf( "\n" ) ;

  // initialisation
  init_latt( DIMS , 6750 , 4 , 750 ) ;

  // number of measurements being made
  const int MEASUREMENTS = ( Latt.MAX_ITERS - Latt.THERM ) / Latt.NMEAS ;

  // init navigation ...
  struct site *lat = (struct site*)malloc( LVOLUME * sizeof( struct site ) ) ;
  init_navig( lat ) ;

  // allocate newfront and oldfront
  newfront = (int*)malloc( LVOLUME * sizeof( int ) ) ;
  oldfront = (int*)malloc( LVOLUME * sizeof( int ) ) ;

  GLU_set_KISS_table( 1234 ) ;

  // set up the stats
  struct resampled raw , rawsq , jack , jacksq , susc ;
  init_stats( &raw , &rawsq , &jack , &jacksq , &susc , MEASUREMENTS ) ;

  // inverse volume
  const double INVOL = 1.0 / (double)LVOLUME ;
 
  // loop betas
  double beta ;
  for( beta = 1.8 ; beta < 2.8 ; beta+= 0.025 ) {
    
    // initialise measurement index
    int index = 0 ; 

    // initialise to a hot start
    //int mag = hot_start( lat ) ;
    int mag = cold_start( lat ) ;
    
    // precompute clustering probability
    const double p = exp( -2.0 / beta ) ;

    // approximately lattice number of spinflips
    for( int k = 0 ; k < Latt.MAX_ITERS ; k++ ) {
      // perform a lattice-wide cluster flip
      int iters = 0 ;
      while( iters < LVOLUME ) {
	iters += wolff( lat , &mag , p ) ;
      }
      // measure every stride
      if( ( k > Latt.THERM ) && ( k%Latt.NMEAS == 0 ) ) {
	raw.resampled[ index ] = abs( mag ) * INVOL ;
	rawsq.resampled[ index ] = ( mag * mag ) * ( INVOL * INVOL ) ;
	index++ ;
      }
      // 
    }

    // resample
    jackknife( &jack , raw ) ;
    jackknife( &jacksq , rawsq ) ;
    equate( &susc , jacksq ) ;

    printf( "%f %e %e " , beta , jack.avg , jack.err ) ;

    // square jack
    for( int i = 0 ; i < jack.NSAMPLES ; i++ ) {
      jack.resampled[ i ] *= jack.resampled[ i ] ;
    }

    subtract( &susc , jack ) ;

    // square jack
    for( int i = 0 ; i < jack.NSAMPLES ; i++ ) {
      susc.resampled[ i ] /= beta ;
    }
    
    printf( "%e %e\n" , susc.avg , susc.err ) ;
  }

  // free the rng and the spins
  GLU_free_KISS_table( ) ;

  // free allocated memory here
  free( lat ) ;
  free( oldfront ) ;
  free( newfront ) ;

  // free the stats
  free( raw.resampled ) ;
  free( rawsq.resampled ) ;
  free( jack.resampled ) ;
  free( jacksq.resampled ) ;

  return 0 ;
}
