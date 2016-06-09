/**
   IO for the gluon propagator ...
 */

#include "fitfunc.h"
#include <stdint.h>
#include "GLU_bswap.h"
#include "GLU_timer.h"

static bool must_swap = true ;

// want a simple reader that sanity checks the timeslice list
static void 
read_tlist( int *NDATA ,
	    FILE *file ,
	    const int *dimensions )
{
  uint32_t size[1] = { 0 } ;
  if( fread( size , sizeof(uint32_t) , 1 , file ) != 1 ) return ;
  if( must_swap ) { bswap_32( 1 , size ) ; }
  *NDATA = size[0] ;
  // time slice list ...
  uint32_t tlist[ *NDATA ] ;
  if( fread( tlist , sizeof(uint32_t) , *NDATA , file ) != *NDATA ) {
    printf( "Tlist read failure \n" ) ;
    exit(1) ;
  }
  printf( "Tlist read\n" ) ;
  return ;
}

// read in nmeas's data from the file
static void
read_data_tcorr( struct resampled *sample , 
		 int *LT ,
		 FILE *file ,
		 const int meas ,
		 const foldselection fold )
{
  int size[1] ;
  if( fread( size , (sizeof( int )) , 1 , file ) != 1 ) exit(1) ;
  if( must_swap ) bswap_32( 1 , size ) ;
  if( size[0] != *LT ) {
    printf( "NDATA-SIZE mismatch %d %d \n" , size[0] , *LT ) ;
    exit(1) ;
  }
  // read in all the doubles
  double *C = malloc( *LT * sizeof( double ) ) ;
  if( fread( C , (sizeof( double )) , *LT , file ) != *LT ) exit(1) ;
  if( must_swap ) bswap_64( *LT , C ) ;
  int t ;
  switch( fold ) {
  case NOFOLD :
    for( t = 0 ; t < *LT ; t++ ) {
      sample[t].resampled[meas] = ( C[t] ) ;
    }
    break ;
  case PlPl :
    sample[0].resampled[meas] = ( C[0] ) ;
    for( t = 1 ; t < *LT/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * ( C[t] + C[*LT-t] ) ;
    }
    break ;
  case PlMi :
    sample[0].resampled[meas] = ( C[0] ) ;
    for( t = 1 ; t < *LT/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * ( C[t] - C[*LT-t] ) ;
    }
    break ;
  case MiPl :
    sample[0].resampled[meas] = -( C[0] ) ;
    for( t = 1 ; t < *LT/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * ( -C[t] + C[*LT-t] ) ;
    }
    break ;
  case MiMi : 
    sample[0].resampled[meas] = -( C[0] ) ;
    for( t = 1 ; t < *LT/2 ; t++ ) {
      sample[t].resampled[meas] = -0.5 * ( C[t] + C[*LT-t] ) ;
    }
    break ;
  }
  free( C ) ;
  return ;
}

// read a raw GLU propagator file
struct resampled*
read_rawGLU_tcorr( struct mom_info **mominfo ,
		   struct input_params *INPARAMS ,
		   const int nfile ,
		   const momtype mom_type ,
		   const foldselection fold )
{
  const int NMEAS = ( INPARAMS -> traj_end[nfile] - INPARAMS -> traj_begin[nfile] ) / 
    INPARAMS -> traj_increment[nfile] ;

  char filestr[ 512 ] ;
  sprintf( filestr , INPARAMS -> traj_file[nfile] , INPARAMS -> traj_begin[nfile] ) ;

  printf( "Reading %s \n" , filestr ) ;

  // read initial momlist
  FILE *file = fopen( filestr , "rb" ) ;
  if( file == NULL ) {
    printf( "Cannot open %s \n" , filestr ) ;
    return NULL ;
  }

  int ndata = 0 ;
  read_tlist( &ndata , file , INPARAMS->dimensions[nfile] ) ;

  print_time() ;

  // set this
  // set this
  if( fold == NOFOLD ) {
    INPARAMS -> NDATA[ nfile ] = ndata ;
  } else {
    INPARAMS -> NDATA[ nfile ] = ndata/2 ;
  }

  fclose( file ) ;

  // allocate the resampled struct
  struct resampled *sample = malloc( ndata * sizeof( struct resampled ) ) ;
  int i ;
  for( i = 0 ; i < ndata ; i++ ) {
    sample[i].resampled = malloc( NMEAS * sizeof( double ) ) ;
    sample[i].restype = RAWDATA ;
    sample[i].NSAMPLES = NMEAS ;
  }

  // loop trajectory files
  //#pragma omp parallel for private(i) 
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
      exit(1) ;
    }

    // read in the momentum list
    int tndata = ndata ;
    read_tlist( &tndata , loc_file , INPARAMS->dimensions[nfile] ) ;

    // set the raw data
    read_data_tcorr( sample , &tndata , loc_file , meas , fold ) ;

    fclose( loc_file ) ;
  }

  return sample ;
}

struct resampled**
read_GLUprop_config( struct mom_info ***mominfo ,
		     double ***X , // is the momentum^2
		     struct input_params *INPARAMS ,
		     int *NSLICES ,
		     const momtype mom_type )
{
  struct resampled **RAW = malloc( INPARAMS -> NFILES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( INPARAMS -> NFILES * sizeof( double ) ) ;
  *mominfo = (struct mom_info**)malloc( INPARAMS -> NFILES * sizeof( struct mom_info* ) ) ;
  int nfile = 0 ;
  for( nfile = 0 ; nfile < INPARAMS -> NFILES ; nfile++ ) {
    RAW[nfile] = read_rawGLU_tcorr( *mominfo+nfile , INPARAMS , nfile , mom_type , INPARAMS -> tfold == true ? PlPl : NOFOLD ) ;
    (*X)[nfile] = (double*)malloc( INPARAMS -> NDATA[nfile] * sizeof( double ) ) ;
    int j ;
    #pragma omp parallel for private(j)
    for( j = 0 ; j < INPARAMS -> NDATA[nfile] ; j++ ) {
      (*X)[nfile][j] = j ;
    }
  }
  *NSLICES = INPARAMS -> NFILES ;
  return RAW ;
}
