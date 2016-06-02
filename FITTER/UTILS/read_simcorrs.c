/**
   @file read_simcorrs.c
   @brief read simultaneous correlators
 */
#include "fitfunc.h"

#include "read_corr.h"

// resampled prop
struct resampled**
read_corr( double ***X , // is the momentum^2
	   struct input_params *INPARAMS ,
	   struct mom_info **mom ,
	   int *NSLICES ,
	   const int src ,
	   const int snk ,
	   const bool tfold )
{
  const int NCHANNELS = INPARAMS->NCHANNELS = 1 ;
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
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , nfiles , 
			       tfold == true ? MiMi : NOFOLD , 
			       src , snk ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
  }
  return RAW ;
}

// OK, simultaneous fit over ( P P ) ( A_t A_t ) ( P A_t )
struct resampled**
read_simcorr( double ***X , // is the tiem
	      struct input_params *INPARAMS ,
	      struct mom_info **mom ,
	      int *NSLICES ,
	      const bool tfold )
{
  const int NCHANNELS = INPARAMS->NCHANNELS = 4 ;
  // use AP and PA ?, it is a bit noisy, set to 4 if you want
  *NSLICES = NCHANNELS * INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  printf( "WHAT the living fuck! %d \n" , INPARAMS->NFILES ) ;
  // copy inparams at the start
  expand_inparams( INPARAMS , NCHANNELS ) ;

  int nfiles , idx = 0 ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
    struct mom_info *tmp = NULL ;
    {
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , idx , 
			       tfold == true ? PlPl : NOFOLD , 
			       5 , 5 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
    if( NCHANNELS > 1 ) {
      RAW[idx] = corrs( X , INPARAMS , mom , str , 
			tfold == true ? PlPl : NOFOLD , 
			9 , 9 , idx ) ; 
      idx++ ;
    }
    if( NCHANNELS > 2 ) {
      RAW[idx] = corrs( X , INPARAMS , mom , str , 
			tfold == true ? MiPl : NOFOLD , 
			5 , 9 , idx ) ; 
      idx++ ;
    }
    if( NCHANNELS > 3 ) {
      RAW[idx] = corrs( X , INPARAMS , mom , str , 
			tfold == true ? MiPl : NOFOLD , 
			9 , 5 , idx ) ; 
      idx++ ;
    }
  }
  return RAW ;
}

// OK, simultaneous fit over ( A_x A_x ) ( A_y A_y ) ( A_z A_z )
struct resampled**
read_simaxcorr( double ***X , // is the tiem
		struct input_params *INPARAMS ,
		struct mom_info **mom ,
		int *NSLICES ,
		const bool tfold )
{
  INPARAMS->NCHANNELS = 1 ;
  *NSLICES = INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  int nfiles ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
    struct mom_info *tmp = NULL ;
    {
      RAW[ nfiles ] = read_rawcorr( INPARAMS , &tmp , str , nfiles , 
				    tfold == true ? MiMi : NOFOLD , 6 , 6 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[ nfiles ] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[ nfiles ][j] = (double)j ; }
    }
    struct resampled *tmp1 = corrs( X , INPARAMS , mom , str , tfold == true ? MiMi : NOFOLD , 7 , 7 , nfiles ) ; 
    struct resampled *tmp2 = corrs( X , INPARAMS , mom , str , tfold == true ? MiMi : NOFOLD , 8 , 8 , nfiles ) ; 
    int j ;
    for( j = 0 ; j < INPARAMS -> NDATA[ nfiles ] ; j++ ) {
      add( &RAW[ nfiles ][ j ] , tmp1[ j ] ) ;
      add( &RAW[ nfiles ][ j ] , tmp2[ j ] ) ;
      mult_constant( &RAW[ nfiles ][ j ] , 1.0/3.0 ) ;
      free( tmp1[j].resampled ) ;
      free( tmp2[j].resampled ) ;
    }
    free( tmp1 ) ;
    free( tmp2 ) ;
  }
  return RAW ;
}

struct resampled**
read_simveccorr( double ***X , // is the tiem
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold )
{
  const int NCHANNELS = INPARAMS->NCHANNELS = 2 ;
  *NSLICES = NCHANNELS * INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  // expand the inparams
  expand_inparams( INPARAMS , NCHANNELS ) ;

  int nfiles , idx = 0 ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
    struct mom_info *tmp = NULL ;
    {
      RAW[ idx ] = read_rawcorr( INPARAMS , &tmp , str , nfiles , tfold == true ? MiMi : NOFOLD , 0 , 0 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[ idx ] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[ idx ][j] = (double)j ; }
    }
    struct resampled *tmp1 = corrs( X , INPARAMS , mom , str , tfold == true ? MiMi : NOFOLD , 1 , 1 , idx ) ; 
    struct resampled *tmp2 = corrs( X , INPARAMS , mom , str , tfold == true ? MiMi : NOFOLD , 2 , 2 , idx ) ; 
    int j ;
    for( j = 0 ; j < INPARAMS -> NDATA[ nfiles ] ; j++ ) {
      add( &RAW[ idx ][ j ] , tmp1[ j ] ) ;
      add( &RAW[ idx ][ j ] , tmp2[ j ] ) ;
      mult_constant( &RAW[ idx ][ j ] , 1.0/3.0 ) ;
      free( tmp1[j].resampled ) ;
      free( tmp2[j].resampled ) ;
    }
    free( tmp1 ) ;
    free( tmp2 ) ;
    idx++ ;

    if( NCHANNELS == 2 ) {
      // tensory ones
      {
	RAW[ idx ] = read_rawcorr( INPARAMS , &tmp , str , idx , tfold == true ? MiMi : NOFOLD , 12 , 12 ) ;
	(*mom)[0] = tmp[0] ;
	(*X)[ idx ] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
	int j ;
	for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[ idx ][j] = (double)j ; }
      }
      tmp1 = corrs( X , INPARAMS , mom , str , tfold == true ? MiMi : NOFOLD , 14 , 14 , idx ) ; 
      tmp2 = corrs( X , INPARAMS , mom , str , tfold == true ? MiMi : NOFOLD , 15 , 15 , idx ) ; 
      for( j = 0 ; j < INPARAMS -> NDATA[ nfiles ] ; j++ ) {
	add( &RAW[ idx ][ j ] , tmp1[ j ] ) ;
	add( &RAW[ idx ][ j ] , tmp2[ j ] ) ;
	mult_constant( &RAW[ idx ][ j ] , 1.0/3.0 ) ;
	free( tmp1[j].resampled ) ;
	free( tmp2[j].resampled ) ;
      }
      free( tmp1 ) ;
      free( tmp2 ) ;
      idx++ ;
    }
  }
  return RAW ;
}

// read simultaneous baryon corr
struct resampled**
read_simbarcorr( double ***X , // is the tiem
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold )
{
  const int NCHANNELS = INPARAMS->NCHANNELS = 1 ;
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
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , nfiles , tfold == true ? PlPl : NOFOLD , 5 , 5 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
  }
  return RAW ;
}

// read simultaneous baryon corr
struct resampled**
read_simvecbarcorr( double ***X , // is the tiem
		    struct input_params *INPARAMS ,
		    struct mom_info **mom ,
		    int *NSLICES ,
		    const bool tfold )
{
  *NSLICES = INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  int nfiles ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
    struct mom_info *tmp = NULL ;
    {
      RAW[ nfiles ] = read_rawcorr( INPARAMS , &tmp , str , nfiles , tfold == true ? PlPl : NOFOLD , 0 , 0 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[ nfiles ] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[ nfiles ][j] = (double)j ; }
    }
    struct resampled *tmp1 = corrs( X , INPARAMS , mom , str , tfold == true ? PlPl : NOFOLD  , 1 , 1 , nfiles ) ; 
    struct resampled *tmp2 = corrs( X , INPARAMS , mom , str , tfold == true ? PlPl : NOFOLD  , 2 , 2 , nfiles ) ; 
    size_t j ;
    for( j = 0 ; j < INPARAMS -> NDATA[ nfiles ] ; j++ ) {
      add( &RAW[ nfiles ][ j ] , tmp1[ j ] ) ;
      add( &RAW[ nfiles ][ j ] , tmp2[ j ] ) ;
      mult_constant( &RAW[ nfiles ][ j ] , 1.0/3.0 ) ;
    }
    // free the temps
    for( j = 0 ; j < INPARAMS -> NDATA[ nfiles ] ; j++ ) {
      free( tmp1[j].resampled ) ;
      free( tmp2[j].resampled ) ;
    }
    free( tmp1 ) ;
    free( tmp2 ) ;
  }
  return RAW ;
}
