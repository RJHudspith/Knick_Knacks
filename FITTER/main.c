// is most important so goes first
#include "fitfunc.h"

#include <stdint.h>

#include "equivalents.h"
#include "GLU_bswap.h"
#include "rng.h"
#include "GLU_timer.h"
#include "read_inputs.h"

// for the fitfunction
#include "fitdata.h" 

// found in UTILS
#include "graph_data.h"
#include "Utils.h"
#include "stats.h" 
#include "shuffalgo.h"

#include "analysis_wrapper.h"
#include "io_wrapper.h"

#include "write_superdist.h"

// always nice to have a usage function
static int
usage( void )
{
  return printf( "USAGE :: ./GLUFIT infile\n" ) ;
}

// the thing that controls it all
int
main( const int argc , char *argv[] )
{
  // exit early if we don't know what's up
  if( argc != 2 ) { return usage() ; }

  struct input_params INPARAMS ;

  start_timer( ) ;

  printf( "\n--> Input File Data <--\n" ) ;

  // read that input file of ours ...
  if( read_inputs( &INPARAMS , argv[1] ) == FAILURE ) {
    return FAILURE ;
  }

  // check the endianness
  BigEndian = is_big_endian( ) ;
  if( BigEndian ) printf( "\n--> We are BIG endian <--\n" ) ;
  else printf( "\n--> We are LITTLE endian <--\n" ) ;

  // set the seed from the input file ...
  init_rng( INPARAMS.Seed ) ;

  // number of simultaneous fits we are performing ...
  int NSLICES = 1 ;
  const int LT = INPARAMS.dimensions[0][3] ;

  // fake data
  int i ;

  // general data 
  struct resampled **RAW = NULL ;
  struct mom_info **mominfo = NULL , *moms = NULL ;
  double **X = NULL ;

  read_data( &X , &RAW , &mominfo , &moms , &INPARAMS , &NSLICES , LT ) ;

  // leave
  if( RAW == NULL ) { printf( "[IO] RAWDATA empty \n" ) ;  return FAILURE ; }
  if( mominfo == NULL ) { printf( "[IO] moms empty \n" ) ;  return FAILURE ; }
  if( X == NULL ) { printf( "[IO] x empty \n" ) ;  return FAILURE ; }

  print_time( ) ;

  printf( "\n--> Sorting Data <--\n" ) ;

  // compute the error and sort the data 
  #pragma omp parallel for private(i)
  for( i = 0 ; i < NSLICES ; i++ ) {
    int k ;
    for( k = 0 ; k < INPARAMS.NDATA[i] ; k++ ) {
      compute_err( &RAW[i][k] ) ;
    }
    // sort the data
    //heapSort( X[i] , RAW[i] , mominfo[i] , INPARAMS.NDATA[i] ) ;
  }

  printf( "\n--> Momaveraging <--\n" ) ;

  // BOOT and X are freed in momavg
  struct mom_info **mominfoavg ;
  struct resampled **RAWAVG = momentum_average( &mominfoavg , &INPARAMS , 
						(const struct mom_info**)mominfo , 
						(const struct resampled**)RAW , 
						NSLICES , INPARAMS.momavg ) ;

  // set up x-avg
  double **xavg = malloc( NSLICES * sizeof( double* ) ) ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    xavg[ i ] = ( double* )malloc( INPARAMS.NDATA[i] * sizeof( double ) ) ;
    int k ;
    for( k = 0 ; k < INPARAMS.NDATA[ i ] ; k++ ) {
      xavg[ i ][ k ] = mominfoavg[ i ][ k ].p2 ;
    }
  }

  #if 0
  super_distribution( RAWAVG ,
		      (const struct mom_info**)mominfoavg , 
		      INPARAMS , NSLICES ,
		      4 , true ) ;
  exit(1) ;
  #endif

  print_time( ) ;

  printf( "\n--> Initialising bootstrap order <--\n" ) ;

  printf( "NSLICE :: %d \n" , NSLICES ) ;

  // initialise the boot order, only need the largest number of RAW data
  // points
  int nsamp = 0 ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    if( RAWAVG[i][0].NSAMPLES > nsamp ) {
      nsamp = RAWAVG[i][0].NSAMPLES ;
    }
  }
  init_boot_order( INPARAMS.NBOOTS , nsamp ) ;

  printf( "\n--> Resampling <--\n" ) ;

  // average and resample , should we do binning in here? maybe
  struct resampled **BOOT = resample_array( (const struct resampled**)RAWAVG , 
					    NSLICES ,
					    INPARAMS.NDATA ,
					    INPARAMS.resample , 
					    INPARAMS.NBOOTS ) ;

  // free the RAW data and the X data
  free_resampled_dp( RAWAVG , NSLICES , INPARAMS.NDATA ) ;

  print_time( ) ;

  printf( "\n--> Fitting <--\n" ) ;

  // analysis wrapper
  perform_analysis( xavg , BOOT , mominfoavg , moms , 
		    &INPARAMS , NSLICES , LT ) ;

  print_time( ) ;

  printf( "\n--> Cleaning Up <--\n" ) ;

  // and free the x-axis
  free_double_dp( xavg , NSLICES ) ;

  // free BOOTS
  free_resampled_dp( BOOT , NSLICES , INPARAMS.NDATA ) ;

  // free the momenta
  free_mominfo( mominfoavg , NSLICES ) ;

  // free the correlator momenta
  free( moms ) ;

  // free the rng
  free_rng( ) ;

  // free the allocated boot list
  free_boot_order( ) ;

  //
  print_time() ;

  return 0 ;
}

