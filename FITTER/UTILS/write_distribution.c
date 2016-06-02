/**
   Writes a UK-hadron readable file
 */

#include "fitfunc.h"

#include "GLU_bswap.h"

void
write_singledist( struct resampled data ,
		  FILE *file )
{
  // write out the resample type and the number of samples
  int deets[ 2 ] = { data.restype , data.NSAMPLES } ;
  if( !BigEndian ) bswap_32( 2 , deets ) ;
  fwrite( deets , sizeof( int ) , 2 , file ) ;
  // and then write out each individual set of bootstraps
  if( !BigEndian ) bswap_64( data.NSAMPLES , data.resampled ) ;
  fwrite( data.resampled , sizeof( double ) , data.NSAMPLES , file ) ;
  if( !BigEndian ) bswap_64( data.NSAMPLES , data.resampled ) ;
  double ave[ 1 ] = { data.avg } ;
  if( !BigEndian ) bswap_64( 1 , ave ) ;
  fwrite( ave , sizeof( double ) , 1 , file ) ;
  return ;
}

void
write_distribution_arr( struct resampled *data ,
			const char *filename ,
			const int NRES )
{
  FILE *file = fopen( filename , "wb" ) ;
  int ndata[1] = { NRES } ;
  if( !BigEndian ) bswap_32( 1 , ndata ) ;
  fwrite( ndata , sizeof( int ) , 1 , file ) ;

  int i ;
  for( i = 0 ; i < NRES ; i++ ) {
    write_singledist( data[i] , file ) ;
  }

  // at the end we should have a simple checksum like the CRC or something
  // and also specify the random number used perhaps?

  return ;
}

void
write_distribution_arr2( const double *X , 
			 struct resampled *data ,
			 const char *filename ,
			 const int NRES )
{
  FILE *file = fopen( filename , "wb" ) ;
  int ndata[1] = { NRES } ;
  if( !BigEndian ) bswap_32( 1 , ndata ) ;
  fwrite( ndata , sizeof( int ) , 1 , file ) ;

  // create a bullshit X array
  struct resampled *XFAKE = malloc( NRES * sizeof( struct resampled ) ) ;

  int i ;
  for( i = 0 ; i < NRES ; i++ ) {
    XFAKE[i].resampled = malloc( data[i].NSAMPLES * sizeof( double ) ) ;
    int j ;
    for( j = 0 ; j < data[i].NSAMPLES ; j++ ) {
      XFAKE[i].resampled[j] = X[i] ;
    }
    XFAKE[i].avg = X[i] ;
    XFAKE[i].NSAMPLES = data[i].NSAMPLES ;
    //printf( "%d \n" , XFAKE[i].NSAMPLES ) ;
    XFAKE[i].restype = data[i].restype ;

    // write each x
    write_singledist( XFAKE[i] , file ) ;
  }

  fwrite( ndata , sizeof( int ) , 1 , file ) ;

  // and write a list of the data
  for( i = 0 ; i < NRES ; i++ ) {
    write_singledist( data[i] , file ) ;
  }

  // free dat shit
  for( i = 0 ; i < NRES ; i++ ) {
    free( XFAKE[i].resampled ) ;
  }
  free( XFAKE ) ;

  // at the end we should have a simple checksum like the CRC or something
  // and also specify the random number used perhaps?

  return ;
}
