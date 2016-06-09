/**
   @file read_tmoments.c
   @brief read the time-moment correlator C[ t ] = \sum_i c_ii(t)
 */
#include "fitfunc.h"
#include <stdint.h>
#include <complex.h>
#include "GLU_bswap.h"
#include "GLUdata_glueprop.h"

#define SUBZERO

// do we have to byte swap?
static bool must_swap = false ;

// for the correlator
typedef enum { TRANSVERSE = 0 , LONGITUDINAL = 1 } projtype ;

// read a single correlator
static int
read_data( struct resampled *sample ,
	   FILE *file ,
	   const int ND ,
	   const int LT ,
	   const int meas ,
	   const projtype proj ,
	   const int dir )
{
  // read magic number
  uint32_t magic[ 1 ] ;
  if( fread( magic , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    printf( "[IO] Magic number read failure \n" ) ;
    return FAILURE ;
  }

  if( magic[ 0 ] != 717685 ) {
    bswap_32( 1 , magic ) ;
    if( magic[ 0 ] != 717685 ) {
      printf( "[IO] magic number read failure \n" ) ;
      return FAILURE ;
    }
    must_swap = true ;
  }
  
  // read munu
  uint32_t in[ 1 ] ;
  if( fread( in , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    printf( "munu read failure \n" ) ;
    return FAILURE ;
  }
  if( must_swap == true ) bswap_32( 1 , in ) ;

  if( in[0] != ND*ND ) {
    printf( "[IO] unexpected munu length %d \n" , in[0] ) ;
    return FAILURE ;
  }

  // read data
  double data[ LT ] ;

  // initialise sample to zero
  {
    size_t t ;
    for( t = 0 ; t < LT ; t++ ) {
      sample[t].resampled[ meas ] = 0.0 ;
    }
  }

  // loop munu
  size_t munu ;
  for( munu = 0 ; munu < ND*ND ; munu++ ) {
    uint32_t lt[ 1 ] ;
    if( fread( lt , sizeof( uint32_t ) , 1 , file ) != 1 ) {
      printf( "LT read failure \n" ) ;
      return FAILURE ;
    }
    if( must_swap == true ) bswap_32( 1 , lt ) ;
    if( lt[ 0 ] != LT ) {
      printf( "[IO] LT mismatch %d %d \n" , lt[0] , LT ) ;
      return FAILURE ;
    }
    // read the data
    if( fread( data , sizeof( double ) , LT , file ) != LT ) {
      printf( "DATA read failure \n" ) ;
      return FAILURE ;
    }
    const size_t mu = munu / ND ;
    const size_t nu = munu % ND ;
    if( mu == nu ) {
      // switch the projections
      // longitudinal component is in the time direction
      size_t t ;
      switch( proj ) {
      case LONGITUDINAL :
	if( mu == dir ) {
	  if( must_swap == true ) bswap_64( LT , data ) ;
	  for( t = 0 ; t < LT ; t++ ) {
	    sample[t].resampled[ meas ] -= ( data[ t ] ) ; 
	  }
	}
	break ;
      case TRANSVERSE :
	if( mu != dir ) {
	  if( must_swap == true ) bswap_64( LT , data ) ;
	  for( t = 0 ; t < LT ; t++ ) {
	    sample[t].resampled[ meas ] += ( data[ t ] ) / (double)( ND-1 ); 
	  }
	}
	break ;
      }
    }
  }

  return SUCCESS ;
}

// read a raw GLU propagator file
struct resampled*
read_rawCORR( struct input_params *INPARAMS ,
	      const int nfile ,
	      const int ndata ,
	      const int ND , 
	      const int LT ,
	      const projtype proj )
{
  const int NMEAS = ( INPARAMS -> traj_end[nfile] - INPARAMS -> traj_begin[nfile] ) / 
    INPARAMS -> traj_increment[nfile] ;

  char filestr[ 512 ] ;
  sprintf( filestr , INPARAMS -> traj_file[nfile] , INPARAMS -> traj_begin[nfile] ) ;

  // allocate the resampled struct
  struct resampled *sample = malloc( LT * sizeof( struct resampled ) ) ;
  int i ;
  for( i = 0 ; i < ndata ; i++ ) {
    sample[i].resampled = malloc( NMEAS * sizeof( double ) ) ;
    sample[i].restype = RAWDATA ;
    sample[i].NSAMPLES = NMEAS ;
  }

  // loop trajectory files
  bool failure = false ;
  #pragma omp parallel for private(i) 
  for( i = INPARAMS -> traj_begin[nfile] ; i < INPARAMS -> traj_end[nfile] ;
       i += INPARAMS -> traj_increment[nfile] ) {

    const int meas = ( i - INPARAMS -> traj_begin[nfile] ) / 
      INPARAMS -> traj_increment[nfile] ;

    char loc_filestr[ 512 ] ;

    // create the file name
    sprintf( loc_filestr , INPARAMS -> traj_file[nfile] , i ) ;

    // read initial momlist
    FILE *loc_file = fopen( loc_filestr , "rb" ) ;
    if( loc_file == NULL ) {
      printf( "Cannot open %s \n" , loc_filestr ) ;
      failure = true ;
    }

    // set the raw data , timelike is first
    //if( read_data( sample , loc_file , ND , LT , meas , proj , ND-1-nfile ) == FAILURE ) {
    if( read_data( sample , loc_file , ND , LT , meas , proj , ND-1 ) == FAILURE ) {
      printf( "Read data failure %d \n" , i ) ;
      failure = true ;
    }

    // close the local file
    fclose( loc_file ) ;
  }
  if( failure == true ) return NULL ;

  return sample ;
}

// computes transverse correlator
struct resampled**
read_tmoments( double ***X ,
	       struct mom_info ***mominfo ,
	       struct input_params *INPARAMS ,
	       int *NSLICES )
{
  const size_t NDATA = 4 * INPARAMS -> NFILES ; // time corr, trans, long
  const size_t ND = 4 ;
  struct resampled **RAW = malloc( NDATA * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( NDATA * sizeof( double ) ) ;
  *mominfo = malloc( NDATA * sizeof( struct mom_info* ) ) ;
  // first set are transverse
  int nfile = 0 ;
  for( nfile = 0 ; nfile < INPARAMS -> NFILES ; nfile++ ) {
    const int LT = INPARAMS -> dimensions[ nfile ][ ND-1 ] ;
    INPARAMS -> NDATA[nfile] = LT ;
    RAW[nfile] = read_rawCORR( INPARAMS , nfile , LT , ND , LT , TRANSVERSE ) ;
    if( RAW[ nfile ] == NULL ) { exit(1) ; return NULL ; }
    (*X)[nfile] = (double*)malloc( INPARAMS -> NDATA[nfile] * sizeof( double ) ) ;
    (*mominfo)[nfile] = malloc( INPARAMS -> NDATA[nfile] * sizeof( struct mom_info ) ) ;
    int j ;
    #pragma omp parallel for private(j)
    for( j = 0 ; j < INPARAMS -> NDATA[nfile] ; j++ ) {
      (*X)[nfile][j] = j ;
      //const int nmu[ 4 ] = { 0 , 0 , 0 , j } ;
      (*mominfo)[ nfile ][ j ].p2 = j ;
    }
  }

  // second are longitudinal
  for( nfile = 0 ; nfile < INPARAMS -> NFILES ; nfile++ ) {
    const int LT = INPARAMS -> dimensions[ nfile ][ ND-1 ] ;
    const size_t fileno = INPARAMS -> NFILES + nfile ;
    INPARAMS -> NDATA[ fileno ] = LT ;
    RAW[ fileno ] = read_rawCORR( INPARAMS , nfile , LT , ND , LT , LONGITUDINAL ) ;
    if( RAW[ fileno ] == NULL ) { exit(1) ; return NULL ; }
    (*X)[ fileno ] = (double*)malloc( INPARAMS -> NDATA[nfile] * sizeof( double ) ) ;
    (*mominfo)[ fileno ] = malloc( INPARAMS -> NDATA[nfile] * sizeof( struct mom_info ) ) ;
    int j ;
    #pragma omp parallel for private(j)
    for( j = 0 ; j < INPARAMS -> NDATA[nfile] ; j++ ) {
      (*X)[fileno][j] = j ;
      (*mominfo)[fileno][ j ].p2 = j ;
    }
  }

  // third and fourth are the momentum-space DFT'd data
  for( nfile = 0 ; nfile < INPARAMS -> NFILES ; nfile++ ) {
    const int LT = INPARAMS -> dimensions[ nfile ][ ND-1 ] ;
    const size_t fileno  = 2 * INPARAMS->NFILES + nfile ;
    const size_t fileno2 = 3 * INPARAMS->NFILES + nfile ;
    INPARAMS -> NDATA[ fileno ] = INPARAMS -> NDATA[ fileno2 ] = LT ;
    RAW[ fileno  ] = malloc( INPARAMS -> NDATA[ nfile ] * sizeof( struct resampled ) ) ;
    RAW[ fileno2 ] = malloc( INPARAMS -> NDATA[ nfile ] * sizeof( struct resampled ) ) ;
    if( RAW[ fileno ] == NULL ) { exit(1) ; return NULL ; }
    if( RAW[ fileno2 ] == NULL ) { exit(1) ; return NULL ; }
    (*X)[ fileno ] = (double*)malloc( INPARAMS -> NDATA[nfile] * sizeof( double ) ) ;
    (*X)[ fileno2 ] = (double*)malloc( INPARAMS -> NDATA[nfile] * sizeof( double ) ) ;
    (*mominfo)[ fileno ] = malloc( INPARAMS -> NDATA[nfile] * sizeof( struct mom_info ) ) ;
    (*mominfo)[ fileno2 ] = malloc( INPARAMS -> NDATA[nfile] * sizeof( struct mom_info ) ) ;
    // perform the DFT
    const double twiddle = 2.0 * M_PI / (double)( INPARAMS -> dimensions[ nfile ][ 3 ] ) ; 

    int pt ;
    #pragma omp parallel for private( pt )
    for( pt = -LT/2 ; pt < LT/2 ; pt++ ) {
      const int p = pt + LT / 2 ;
      RAW[ fileno ][ p ].resampled = malloc( RAW[ nfile ][ 0 ].NSAMPLES * sizeof( double ) ) ;
      RAW[ fileno2 ][ p ].resampled = malloc( RAW[ nfile ][ 0 ].NSAMPLES * sizeof( double ) ) ;
      equate_constant( &RAW[ fileno ][ p ] , 0.0 , RAW[ nfile ][ 0 ].NSAMPLES ,  RAW[ nfile ][ 0 ].restype ) ;
      equate_constant( &RAW[ fileno2 ][ p ] , 0.0 , RAW[ nfile ][ 0 ].NSAMPLES , RAW[ nfile ][ 0 ].restype ) ;
      double complex epimu[ LT ] ;
      size_t t , k ;
      for( t = 0 ; t < LT ; t++ ) {
	epimu[ t ] = cos( twiddle * pt * t ) + I * sin( twiddle * pt * t ) ;
      }

      // loop boots performing zero sum subtraction
      register double complex sum , sum2 , zerosum , zerosum2 ;
      for( k = 0 ; k < RAW[ nfile ][ 0 ].NSAMPLES ; k++ ) {
	sum = sum2 = zerosum = zerosum2 = 0.0 ;
	for( t = 0 ; t < LT ; t++ ) {
	  sum += epimu[ t ] * RAW[ nfile ][ t ].resampled[ k ] ;
	  sum2 += epimu[ t ] * RAW[ nfile + INPARAMS -> NFILES ][ t ].resampled[ k ] ;

	  #ifdef SUBZERO
	  zerosum += RAW[ nfile ][ t ].resampled[ k ] ;
	  zerosum2 += RAW[ nfile + INPARAMS -> NFILES ][ t ].resampled[ k ] ;
	  #endif
	}
	RAW[ fileno ][ p ].resampled[ k ] = creal( sum - zerosum ) ;
	RAW[ fileno2 ][ p ].resampled[ k ] = creal( sum2 - zerosum2 ) ;
      }

      // fill in the momenta
      int nmu[ 4 ] = { 0 , 0 , 0 , pt } ;
      /*
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	if( mu == ND-nfile-1 ) { 
	  nmu[mu] = pt ;
	} else {
	  nmu[mu] = 0 ;
	}
      }
      */
      (*mominfo)[ fileno ][ p ] = fill_mominfo( 4 , nmu , INPARAMS->dimensions[nfile] , 
						INPARAMS->mom_type ) ;
      (*mominfo)[ fileno2 ][ p ] = fill_mominfo( 4 , nmu , INPARAMS->dimensions[nfile] , 
						 INPARAMS->mom_type ) ;
      
      // divide by q^2
      mult_constant( &RAW[ fileno ][ p ] , 1.0 / ( (*mominfo)[ fileno ][ p ].p2 > 1E-15 ? (*mominfo)[ fileno ][ p ].p2 : 1 ) ) ;
      mult_constant( &RAW[ fileno2 ][ p ] , 1.0 / ( (*mominfo)[ fileno ][ p ].p2 > 1E-15 ? (*mominfo)[ fileno ][ p ].p2 : 1 ) ) ;

      // set x
      (*X)[fileno][ p ] = sqrt( (*mominfo)[fileno][p].p2 ) ;
      (*X)[fileno2][ p ] = sqrt( (*mominfo)[fileno2][p].p2 ) ;
    }

    printf( "TEST!!!! %e should be 0 !! \n" , (*X)[fileno][ LT/2 ] ) ;
  }

  // set the number of datasets we have
  *NSLICES = NDATA ;

  return RAW ;
}
