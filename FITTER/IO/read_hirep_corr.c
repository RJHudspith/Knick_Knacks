#include "fitfunc.h"
#include <stdint.h>

#include "GLU_timer.h"
#include "read_corr.h"

// ewwww! hirep correlator is plain text that I have to fish out
// of stdout, fucking disgusting!
static int
read_hicorrfile( struct resampled *sample ,
	       struct mom_info **mominfo ,
	       FILE *file ,
	       const int meas ,
	       const size_t LT ,
	       const foldselection fold )
{
  // loop file until we hit an EOF
  size_t t = 0 ;
  double C[ LT ] ;
  while( fscanf( file , " %le" , &C[t] ) != EOF ) {
    if( t > LT ) break ;
    t++ ;
  }
  if( t != (LT) ) {
    printf( "File of wrong length %zu vs. %zu \n" , t , LT-1 ) ;
    return FAILURE ;
  }
  switch( fold ) {
  case NOFOLD :
    for( t = 0 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = C[t] ;
    }
    break ;
  case PlPl :
    sample[0].resampled[meas] = fabs( ( C[0] ) ) ;
    for( t = 1 ; t < LT/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * ( C[t] + C[LT-t] ) ;
    }
    break ;
  case PlMi :
    sample[0].resampled[meas] = ( C[0] ) ;
    for( t = 1 ; t < LT/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * ( C[t] - C[LT-t] ) ;
    }
    break ;
  case MiPl :
    sample[0].resampled[meas] = -( C[0] ) ;
    for( t = 1 ; t < LT/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * ( -C[t] + C[LT-t] ) ;
    }
    break ;
  case MiMi : 
    sample[0].resampled[meas] = -( C[0] ) ;
    for( t = 1 ; t < LT/2 ; t++ ) {
      sample[t].resampled[meas] = -0.5 * ( C[t] + C[LT-t] ) ;
    }
    break ;
  }
  return SUCCESS ;
}

// read a raw correlator
static struct resampled *
read_hirawcorr( struct input_params *INPARAMS ,
		const char *filename ,
		const int nfile ,
		const foldselection fold ) 
{
  const size_t LT = INPARAMS -> dimensions[ nfile ][ 3 ] ;

  const int NMEAS = ( INPARAMS -> traj_end[nfile] - INPARAMS -> traj_begin[nfile] ) / 
    INPARAMS -> traj_increment[nfile] ;

  char filestr[ 256 ] ;
  sprintf( filestr , filename , INPARAMS -> traj_begin[nfile] ) ;

  // file open
  FILE *file = fopen( filestr , "rb" ) ;
  if( file == NULL ) {
    printf( "Cannot open %s \n" , filestr ) ;
    return NULL ;
  }
  print_time() ;

  // set this
  if( fold == NOFOLD ) {
    INPARAMS -> NDATA[ nfile ] = (int)LT ;
  } else {
    INPARAMS -> NDATA[ nfile ] = (int)LT/2 ;
  }

  fclose( file ) ;

  // allocate the resampled struct
  struct resampled *sample = malloc( INPARAMS -> NDATA[ nfile ]  * 
				     sizeof( struct resampled ) ) ;
  int i ;
  for( i = 0 ; i < INPARAMS -> NDATA[ nfile ] ; i++ ) {
    sample[i].resampled = malloc( NMEAS * sizeof( double ) ) ;
    sample[i].restype = RAWDATA ;
    sample[i].NSAMPLES = NMEAS ;
  }

  // loop trajectory files
  bool failure = false ;
  #pragma omp parallel for private(i) 
  for( i = INPARAMS -> traj_begin[nfile] ; 
       i < INPARAMS -> traj_end[nfile] ;
       i += INPARAMS -> traj_increment[nfile] ) {

    const int meas = ( i - INPARAMS -> traj_begin[nfile] ) / 
      INPARAMS -> traj_increment[nfile] ;

    char loc_filestr[ 256 ] ;

    // create the file name
    sprintf( loc_filestr , filename , i ) ;

    // read initial momlist
    FILE *loc_file = fopen( loc_filestr , "r" ) ;
    if( loc_file == NULL ) {
      printf( "Cannot open %s \n" , loc_filestr ) ;
      failure = true ;
    }

    // set the raw data
    struct mom_info *tmp = NULL ;
    read_hicorrfile( sample , &tmp , loc_file , meas , LT , fold ) ;
    free( tmp ) ;

    fclose( loc_file ) ;
  }

  return ( failure == false ) ? sample : NULL ;
}

// resampled prop
struct resampled**
read_hirep_corr( double ***X , // is the momentum^2
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold )
{
  const int NCHANNELS = 1 ;
  // use AP and PA ?, it is a bit noisy, set to 4 if you want
  *NSLICES = NCHANNELS * INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  expand_inparams( INPARAMS , NCHANNELS ) ;

  int nfiles , idx = 0 ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
    {
      RAW[idx] = read_hirawcorr( INPARAMS , str , nfiles , 
				 tfold == true ? PlPl : NOFOLD ) ;

      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;

      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
  }
  return RAW ;
}
