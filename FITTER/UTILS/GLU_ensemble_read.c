/**
   Reads in the data file
 */
#include "fitfunc.h"
#include "GLU_bswap.h"

// this is the reader for the rsq list ...
// if it fucks up we give RAW[0][0].NSAMPLES a value of -1 as the signalling error
struct resampled **
GLUdata( const char *filename , // filename (minus the number and .bin)
	 const int start ,      // trajectory start
	 const int end ,        // trajectory end
	 const int increment ,  // trajectory increment
	 const int NDATA   ,    // number of x's
	 const int NSLICES ,    // number of continguous arrays 
	 const bool POLYCORR ,  // are we measuring the POLY corr? 
	 const int ND           // dimensions of the MOM-list
	 )
{
  // compute the number of raw data
  const int NRAW = ( end - start ) / increment ;

  printf( "NRAW :: %d %d %d %d\n" , start , end , increment , NRAW ) ;

  // allocate the RAW data
  int i ;

  struct resampled **RAW = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    RAW[ i ] = ( struct resampled* )malloc( NDATA * sizeof( struct resampled ) ) ;
    int j ;
    for( j = 0 ; j < NDATA ; j++ ) {
      RAW[ i ][ j ].resampled = malloc( NRAW * sizeof( double ) ) ;
    }
  }

  // measurement index
  int meas = 0 ;

  // increment the trajectories
  int traj = start ;
  for( traj = start ; traj < end ; traj += increment ) {

    // OK so we always end the code with a .bin
    char str[ 256 ] ;
    sprintf( str , "%s.%d.bin" , filename , traj ) ;

    // this is the file we are reading in ...
    FILE *file = fopen( str , "rb" ) ;
    if( file == NULL ) {
      printf( "File opening failure :: %s \n" , str ) ;
      RAW[0][0].NSAMPLES = -1 ;
      return RAW ;
    }
    
    // read in the size ...
    int size[1] ;
    if( fread( size , ( sizeof(int) ) , 1 , file ) != 1 ) exit(1) ;
    if( !BigEndian ) { bswap_32( 1 , size ) ; }
    if( size[0] != NDATA ) {
      printf( "[INITIAL] File length disagreement !!! %d %d \n" , 
	      size[0] , NDATA ) ;
      RAW[0][0].NSAMPLES = -1 ;
      return RAW ;
    }
    
    // skip the remaining x-data, we already have these 
    for( i = 0 ; i < NDATA ; i++ ) {
      if( POLYCORR != 1 ) {
	int temp[1+ND] ;
	if( fread( temp , (sizeof(int)) , 1+ND , file ) != ND+1 ) exit(1) ;
      } else {
	// polyakov loop correlator is all doubles ...
	double temp[1] ;
	if( fread( temp , (sizeof(double)) , 1 , file ) != 1 ) exit(1) ;
      }
    }

    // this bit tells us how many t-correlators there are
    if( POLYCORR == 1 ) {
      if( fread( size , ( sizeof(int) ) , 1 , file ) != 1 ) exit(1) ;
      if( !BigEndian ) { bswap_32( 1 , size ) ; }
    }

    int t ;
    for( t = 0 ; t < NSLICES ; t++ ) {
      // read in the size again
      if( fread( size , (sizeof(int)) , 1 , file ) != 1 ) exit(1) ;
      if( !BigEndian ) { bswap_32( 1 , size ) ; }
      if( size[0] != NDATA ) {
	printf( "[SLICE %d] File length disagreement !!! %d %d \n" , 
		t , size[0] , NDATA ) ;
	RAW[0][0].NSAMPLES = -1 ;
	return RAW ;
      }
      
      // and finally read in one set of raw data
      double *yvalues = malloc( NDATA * sizeof( double ) ) ;
      if( fread( yvalues , sizeof( double ) , NDATA , file ) != NDATA ) exit(1);
      if( !BigEndian ) { bswap_64( NDATA , yvalues ) ; }
      
      // and stick into the resampled structure
      // horrible data access I know
      for( i = 0 ; i < NDATA ; i++ ) {
	RAW[ t ][ i ].resampled[ meas ] = yvalues[ i ] ;
	//printf( "%e \n" , RAW[ t ][ i ].resampled[ meas ] ) ;
      }
    }

    // and finally, that is it
    meas++ ;

    fclose( file ) ;
  }  

  // set the number of resamples
  int t ;
  for( t = 0 ; t < NSLICES ; t++ ) {
    for( i = 0 ; i < NDATA ; i++ ) {
      RAW[ t ][ i ].NSAMPLES = NRAW ;
    }
  }

  return RAW ;
}
