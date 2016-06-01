/**
   @file 2DIsing.c
   @brief metropolis update of the 2D Ising model using MPI
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "stats.h"

#define ND 2

//#define HEATBATH

struct rng
{
  uint32_t a , b , c , d , e ;
} ;

static uint32_t
JKISS32( struct rng *__restrict RNG )
{
  RNG -> b ^= ( RNG -> b << 5 ) ;
  RNG -> b ^= ( RNG -> b >> 7 ) ;
  RNG -> b ^= ( RNG -> b << 22 ) ;
  const uint32_t t = RNG -> c + RNG -> d + RNG -> e ;
  RNG -> c = RNG -> d ;
  RNG -> e = t < 0 ;
  RNG -> d = t&2147483647 ;
  RNG -> a += 1411392427 ;
  return RNG -> a + RNG -> b + RNG -> d ;
}

double
KISS_dbl( struct rng *__restrict RNG ) {
#ifdef FULL_DOUBLE
  // to remove the type-punning error gcc always gives about strict aliasing
  union {
    uint32_t theInts[2] ;
    uint64_t ulint ;
    double x ;
  } lu2dbl ;
  // pack the bottom and top half of an unsigned long with uints
  lu2dbl.theInts[0] = JKISS32( RNG ) ;
  lu2dbl.theInts[1] = JKISS32( RNG ) ;
  // shift down so that only the mantissa has non-zero entries
  lu2dbl.ulint = ( lu2dbl.ulint >> 12 ) | 0x3FF0000000000000ULL ;
  // have to subtract one because of the implicit 2^0+mantissa
  return lu2dbl.x - 1.0 ;
#else
  return JKISS32( RNG ) * 2.3283064365386963e-10 ;
#endif
}

// RNG initialise
static void
init_rng( struct rng *RNG , const uint32_t seed )
{
  RNG -> a = seed ;
  RNG -> b = ( 1812433253UL * ( RNG -> a ^ ( RNG -> a >> 30)) + 1 ) ;
  RNG -> c = ( 1812433253UL * ( RNG -> b ^ ( RNG -> b >> 30)) + 2 ) ;
  RNG -> d = ( 1812433253UL * ( RNG -> c ^ ( RNG -> c >> 30)) + 3 ) ;
  RNG -> e = 0 ;
}

// five probabilities
static double probs[ 5 ] ;

// precompute
static void
precompute_props( const double J ) {
#ifdef HEATBATH
  probs[4] = exp( -J * 8. ) / ( exp( -J * 8. ) + 1. ) ;
  probs[3] = exp( -J * 4. ) / ( exp( -J * 4. ) + 1. ) ;
  probs[2] = 1. / 2. ;
  probs[1] = exp( J * 4. ) / ( exp( J * 4. ) + 1. ) ;
  probs[0] = exp( J * 8. ) / ( exp( J * 8. ) + 1. ) ;
#else
  probs[4] = exp( -J * 8. ) ;
  probs[3] = exp( -J * 4. ) ;
  probs[2] = 1.0 ;
  probs[1] = exp( J * 4. ) ;
  probs[0] = exp( J * 8. ) ;
#endif
}

// compute the geometry
static void 
get_mom_2piBZ( int x[ ND ] , 
	       const int i , 
	       const int sub_LATT[ND] )
{
  int mu , subvol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    x[ mu ] = ( ( i - i % subvol ) / subvol ) % sub_LATT[ mu ] ;
    subvol *= sub_LATT[ mu ] ;
  }
  return ;
}

// compute the local site
const int
loc_site( int x[ ND ] , const int sub_LATT[ND] )
{
  int res = x[ ND - 1 ] ;
  int mu ;
  for( mu = ND - 1 ; mu > 0 ; mu-- ) {
    res = sub_LATT[ mu - 1 ] * res + x[ mu - 1 ] ;
  }
  return res ;
}

// This is the new bona-fide generic shifting code
static int 
gen_shift( const int i , 
	   const int sub_LATT[ND] ,
	   const int dir )
{
  int x[ ND ] ; 
  get_mom_2piBZ( x , i , sub_LATT ) ;
  if( dir >= 0 ) {
    x[ dir ] = ( x[ dir ] + 1 ) % sub_LATT[ dir ] ;
  } else {
    register const int numu = -dir - 1 ;
    if( x[ numu ] == 0 ) {
      x[ numu ] = x[ numu ] - 1 + sub_LATT[ numu ] ;
    } else {
      x[ numu ] = ( x[ numu ] - 1 ) % sub_LATT[ numu ] ;
    }
  }
  return loc_site( x , sub_LATT ) ; 
}


// metropolis
static int
metropolis( int *__restrict lat ,
	    int *__restrict halo_up ,
	    int *__restrict halo_down ,
	    const int sub_LATT[2] ,
	    struct rng *RNG ,
	    const int site )
{
  // sum the neighbours of site in the "x" direction
  int sum = 0 ;

  sum += lat[ gen_shift( site , sub_LATT , 0 ) ] ;
  sum += lat[ gen_shift( site , sub_LATT , -1 ) ] ;

  // the y's are a bit more complicated
  const int temp1 = gen_shift( site , sub_LATT , 1 ) ;

  // we have wrapped around and so must communicate with the upper halo
  if( temp1 < site ) {
    sum += halo_up[ temp1 ] ;
  } else {
    sum += lat[ temp1 ] ;
  }

  const int temp2 = gen_shift( site , sub_LATT , -2 ) ;
  // we have wrapped round at the bottom
  if( temp2 > site ) {
    sum += halo_down[ site ] ;
  } else {
    sum += lat[ temp2 ] ;
  }
  sum *= lat[ site ] ;

  // spinflip dynamics
#ifdef HEATBATH
    if( KISS_dbl( RNG ) <= probs[sum/2+2] ) {
      goto SPIN_FLIP ;
    }
    return 0 ;
#else
  // metropolis
  if( sum > 0 ) {
    if( KISS_dbl( RNG ) <= probs[sum/2+2] ) {
      goto SPIN_FLIP ;
    }
    return 0 ;
  }
#endif

 SPIN_FLIP :
  lat[site] = -lat[site] ;
  return ( 2 * lat[site] ) ;
}

// metro
static int
metropolis_update( int *__restrict halo_up , 
		   int *__restrict halo_down , 
		   int *__restrict latt ,
		   struct rng *RNG , const int sub_LATT[2] , 
		   const int SUB_LVOLUME )
{
  // even sites
  int i , sum = 0 ;
  for( i = 0 ; i < SUB_LVOLUME ; i+=2 ) {
    sum += metropolis( latt , halo_up , halo_down , sub_LATT , RNG , i ) ;
  }
  MPI_Barrier( MPI_COMM_WORLD ) ;
  // odd sites
  for( i = 1 ; i < SUB_LVOLUME ; i+=2 ) {
    sum += metropolis( latt , halo_up , halo_down , sub_LATT , RNG , i ) ;
  }
  MPI_Barrier( MPI_COMM_WORLD ) ;
  return sum ;
}

//
static void
update_halos( int *halo_up , int *halo_down , int *latt , 
	      const int rank , const int numtasks , 
	      const int SUB_LVOLUME , const int *LATT )
{
  // communicate
  const int up = (rank+1)%numtasks ;
  const int down = (rank-1)<0?numtasks-1:rank-1 ;

  // ok so we have a halo above and a halo below
  MPI_Status status ;
  MPI_Sendrecv( latt , LATT[0] , MPI_INT , down , 0 ,
		halo_up , LATT[0] , MPI_INT , up , 0 ,
		MPI_COMM_WORLD , &status ) ;

  MPI_Sendrecv( ( latt + SUB_LVOLUME - LATT[0] ) , LATT[0] , MPI_INT , up , 0 ,
		halo_down , LATT[0] , MPI_INT , down , 0 ,
		MPI_COMM_WORLD , &status ) ;

  // wait until they are all here
  MPI_Barrier( MPI_COMM_WORLD ) ;
  return ;
}

// start hot
static int
hotstart( int *latt , struct rng *RNG , const int SUB_LVOLUME )
{
  int i , loc_sum = 0 ;
  for( i = 0 ; i < SUB_LVOLUME ; i++ ) {
    latt[i] = KISS_dbl( RNG ) > 0.5 ? 1.0 : -1.0 ;
    loc_sum += latt[i] ;
  }
  MPI_Barrier( MPI_COMM_WORLD ) ;
  return loc_sum ;
}

// start cold
static int
coldstart( int *latt , const int SUB_LVOLUME )
{
  int i , loc_sum = 0 ;
  for( i = 0 ; i < SUB_LVOLUME ; i++ ) {
    latt[i] = 1.0 ;
    loc_sum += latt[i] ;
  }
  MPI_Barrier( MPI_COMM_WORLD ) ;
  return loc_sum ;
}

// mainprog
int 
main( int argc , 
      char *argv[] ) 
{
  // if the number of arguments are wrong
  if( argc != 2 ) {
    return -1 ;
  }

  // some initialisations
  const int therm = 500 ;
  const double T = atof( argv[1] ) ;
  const double J = 1.0 / T ;

  // sweet, have the props
  precompute_props( J ) ;

  // some stuff
  const int NPROCS = 4 ;
  const int LATT[2] = { 1536 , 1536 } ;
  const int sub_LATT[2] = { LATT[0] , LATT[1]/NPROCS } ;
  const int SUB_LVOLUME = sub_LATT[0] * sub_LATT[1] ;

  // measurements
  const int NMEAS = 250 ;
  const int STEP = 2 ;

  MPI_Init( &argc , &argv ) ;

  int numtasks , rank , len ;
  char hostname[MPI_MAX_PROCESSOR_NAME] ;

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(hostname, &len);

  // locally inited rng
  struct rng RNG ;
  init_rng( &RNG , 2*rank ) ;

  // local lattice initialised
  int *latt = malloc( SUB_LVOLUME * sizeof( int ) ) ;

  // halos are defined as above and below the slice we are on
  int *halo_up   = malloc( LATT[0] * sizeof( int ) ) ;
  int *halo_down = malloc( LATT[0] * sizeof( int ) ) ;

  //int local_sum = hotstart( latt , &RNG , SUB_LVOLUME ) ;
  int local_sum = coldstart( latt , SUB_LVOLUME ) ;

  // update the halos
  update_halos( halo_up , halo_down , latt , 
		rank , numtasks , SUB_LVOLUME , LATT ) ;

  // metropolis warm up
  int k ; 
  for( k = 0 ; k < therm ; k++ ) {
    local_sum += metropolis_update( halo_up , halo_down , latt , 
				    &RNG , sub_LATT , SUB_LVOLUME ) ;
    update_halos( halo_up , halo_down , latt , rank , 
		  numtasks , SUB_LVOLUME , LATT ) ;
  }

  // set up the stats
  struct resampled raw , rawsq , jack , jacksq , susc ;
  init_stats( &raw , &rawsq , &jack , &jacksq , &susc , NMEAS ) ;
  
  // loop measurements
  for( k = 0 ; k < NMEAS ; k++ ) {
    // don't measure every update
    int j ;
    for( j = 0 ; j < STEP ; j++ ) {
      // more updates
      local_sum += metropolis_update( halo_up , halo_down , latt , 
				      &RNG , sub_LATT , SUB_LVOLUME ) ;
      update_halos( halo_up , halo_down , latt , rank , 
		    numtasks , SUB_LVOLUME , LATT ) ;
    }
    // Reduce all of the local sums into the global sum
    int global_sum;
    MPI_Reduce( &local_sum , &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    raw.resampled[ k ] = abs( global_sum ) / (double)( numtasks * SUB_LVOLUME ) ;
    rawsq.resampled[ k ] = raw.resampled[ k ] * raw.resampled[ k ] ;
    //
  }

  // Print the result
  if ( rank == 0) {
    // resample
    jackknife( &jack , raw ) ;
    jackknife( &jacksq , rawsq ) ;
    equate( &susc , jacksq ) ;
    printf( "%f %e %e " , T , jack.avg , jack.err ) ;
    // square jack
    int i ;
    for( i = 0 ; i < jack.NSAMPLES ; i++ ) {
      jack.resampled[ i ] *= jack.resampled[ i ] ;
    }
    subtract( &susc , jack ) ;
    // square jack
    for( i = 0 ; i < jack.NSAMPLES ; i++ ) {
      susc.resampled[ i ] /= T ;
    }
    printf( "%e %e\n" , susc.avg , susc.err ) ;
  }

  // free the lattice
  free( latt ) ;
  free( halo_up ) ;
  free( halo_down ) ;

  // free the stats
  free( raw.resampled ) ;
  free( rawsq.resampled ) ;
  free( jack.resampled ) ;
  free( jacksq.resampled ) ;

  MPI_Finalize( ) ;

  return 0 ;
}
