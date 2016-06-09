/**
   @file read_gevpcorr.c
   @brief read correlators for the gevp
 */
#include "fitfunc.h"

#include "read_corr.h"

// OK, simultaneous fit over ( P P ) ( A_t A_t ) ( P A_t )
struct resampled**
read_gevpcorr( double ***X , // is the tiem
	      struct input_params *INPARAMS ,
	      struct mom_info **mom ,
	      int *NSLICES ,
	      const bool tfold )
{
  const int NCHANNELS = 4 ;
  // use AP and PA ?, it is a bit noisy, set to 4 if you want
  *NSLICES = NCHANNELS * INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  expand_inparams( INPARAMS , NCHANNELS ) ;

  int nfiles , idx = 0 ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    sprintf( str , "%s" , INPARAMS -> traj_file[ idx ] ) ;
    struct mom_info *tmp = NULL ;
    {
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , idx , tfold , 5 , 5 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
    RAW[idx] = corrs( X , INPARAMS , mom , str , tfold , 5 , 9 , idx ) ; 
    idx++ ;
    RAW[idx] = corrs( X , INPARAMS , mom , str , tfold , 9 , 5 , idx ) ; 
    idx++ ;
    RAW[idx] = corrs( X , INPARAMS , mom , str , tfold , 9 , 9 , idx ) ; 
    idx++ ;
  }
  return RAW ;
}

// OK, use ( P P ) 
struct resampled**
read_gevpcorr2( double ***X , // is the tiem
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
    struct mom_info *tmp = NULL ;
    {
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , nfiles , tfold , 9 , 9 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
  }
  return RAW ;
}
