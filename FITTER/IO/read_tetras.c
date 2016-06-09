/**
   @file read_simcorrs.c
   @brief read simultaneous correlators
 */
#include "fitfunc.h"

#include "read_corr.h"

static void
average_similar( struct resampled **RAW ,
		 double ***X , 
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 const int *opmap ,
		 char str[ 256 ] ,
		 const int idx ,
		 const bool tfold )
{
  struct resampled *tmp1 = corrs( X , INPARAMS , mom , 
				  str , tfold == true ? PlPl : NOFOLD , 
				  opmap[0] , 1 , idx ) ; 
  struct resampled *tmp2 = corrs( X , INPARAMS , mom , 
				  str , tfold == true ? PlPl : NOFOLD , 
				  opmap[1] , 2 , idx ) ; 
  size_t j ;
  for( j = 0 ; j < INPARAMS -> NDATA[ idx ] ; j++ ) {
    add( &RAW[ idx ][ j ] , tmp1[ j ] ) ;
    add( &RAW[ idx ][ j ] , tmp2[ j ] ) ;
    mult_constant( &RAW[ idx ][ j ] , 1.0/3.0 ) ;
    free( tmp1[j].resampled ) ;
    free( tmp2[j].resampled ) ;
  }
  free( tmp1 ) ;
  free( tmp2 ) ;
}

// OK, simultaneous fit over ( P P ) ( A_t A_t ) ( P A_t )
struct resampled**
read_tetras( double ***X , // is the tiem
	     struct input_params *INPARAMS ,
	     struct mom_info **mom ,
	     int *NSLICES ,
	     const bool tfold )
{
  if( INPARAMS -> NFILES == 3 ) {
    *NSLICES = 4 ;
  } else if( INPARAMS -> NFILES == 6 ) {
    *NSLICES = 8 ;
  } else if( INPARAMS -> NFILES == 1 ) {
    *NSLICES = 2 ; 
  } else {
    printf( "Must be a multiple of 3\n" ) ;
    exit(1) ;
  }
  //const size_t NCHANNELS = 4 ; // DiqDiq, Dimes, V , P
  const size_t NOPS = 2 ;
  const size_t opmap[ 2 ] = { 0 , 5 } ; 
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  // copy these right
  if( INPARAMS -> NFILES == 6 ) {
    copy_inparams( INPARAMS , 6 , 8 ) ;
    copy_inparams( INPARAMS , 5 , 7 ) ;
    copy_inparams( INPARAMS , 4 , 6 ) ;
    copy_inparams( INPARAMS , 4 , 5 ) ;
    copy_inparams( INPARAMS , 3 , 4 ) ;
    copy_inparams( INPARAMS , 2 , 3 ) ;
    copy_inparams( INPARAMS , 0 , 1 ) ;
    copy_inparams( INPARAMS , 0 , 0 ) ;
  } else if( INPARAMS -> NFILES == 3 ) {
    copy_inparams( INPARAMS , 3 , 4 ) ;
    copy_inparams( INPARAMS , 2 , 3 ) ;
    copy_inparams( INPARAMS , 0 , 1 ) ;
    copy_inparams( INPARAMS , 0 , 0 ) ;
  } else if( INPARAMS -> NFILES == 1 ) {
    copy_inparams( INPARAMS , 0 , 1 ) ;
    copy_inparams( INPARAMS , 0 , 0 ) ;
  }

  int nfiles = 0 , idx = 0 ;
  char str[ 256 ] ;
  {
    sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
    size_t op , nop ;
    for( nop = 0 ; nop < NOPS ; nop++ ) {
      op = opmap[ nop ] ;
      struct mom_info *tmp = NULL ;
      {
	RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , nfiles , 
				 tfold == true ? PlPl : NOFOLD , 
				 op , 0 ) ;
	(*mom)[0] = tmp[0] ;
	(*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
	
	int ops[2] = { op , op } ;
	average_similar( RAW , X , INPARAMS , mom , ops , str , idx , tfold ) ;
	idx++ ;
      }
    }
    // loop of operators
  }
  if( INPARAMS -> NFILES > 1 ) {
    // vector
    {
      nfiles = 1 ;
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
      int ops[2] = { 1 , 2 } ;
      average_similar( RAW , X , INPARAMS , mom , ops , str , idx , tfold ) ;
      idx++ ;
    }
    // pseudoscalar
    {
      sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
      struct mom_info *tmp = NULL ;
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , idx , 
			       tfold == true ? PlPl : NOFOLD , 
			       5 , 5 ) ;
      (*mom)[0] = tmp[0] ;
      (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
      int j ;
      for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
      idx++ ;
    }
    
    if( INPARAMS -> NFILES == 6 ) {
      nfiles = 3 ;
      {
	sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
	size_t op , nop ;
	for( nop = 0 ; nop < NOPS ; nop++ ) {
	  op = opmap[ nop ] ;
	  struct mom_info *tmp = NULL ;
	  {
	    RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , nfiles , 
				     tfold == true ? PlPl : NOFOLD , 
				     op , 0 ) ;
	    (*mom)[0] = tmp[0] ;
	    (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
	    int j ;
	    for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }

	    int ops[2] = { op , op } ;
	    average_similar( RAW , X , INPARAMS , mom , ops , str , idx , tfold ) ;

	    idx++ ;
	  }
	}
	// loop of operators
      }
      // vector
      nfiles = 4 ;
      {
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
	int ops[2] = { 1 , 2 } ;
	average_similar( RAW , X , INPARAMS , mom , ops , str , idx , tfold ) ;	

	idx++ ;
      }
      nfiles = 5 ;
      // pseudoscalar
      {
	sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
	struct mom_info *tmp = NULL ;
	RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , idx , 
				 tfold == true ? PlPl : NOFOLD , 
				 5 , 5 ) ;
	(*mom)[0] = tmp[0] ;
	(*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
	int j ;
	for( j = 0 ; j < INPARAMS -> NDATA[nfiles] ; j++ ) { (*X)[idx][j] = (double)j ; }
	idx++ ;
      }
    }
  }
  printf( "[Tetras] IO finished\n" ) ;
  return RAW ;
}

// OK mesonic part of the tetra on its own
struct resampled**
read_tetra_meson( double ***X , // is the tiem
		  struct input_params *INPARAMS ,
		  struct mom_info **mom ,
		  int *NSLICES ,
		  const bool tfold )
{
  *NSLICES = 2 * INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  copy_inparams( INPARAMS , 1 , 3 ) ;
  copy_inparams( INPARAMS , 1 , 2 ) ;
  copy_inparams( INPARAMS , 0 , 1 ) ;
  copy_inparams( INPARAMS , 0 , 0 ) ;

  int idx = 0 , nfiles ;
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
    int ops[2] = { 1 , 2 } ;
    average_similar( RAW , X , INPARAMS , mom , ops , str , idx , tfold ) ;
    idx++ ;
  
    // pseudoscalar
    {
      sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
      struct mom_info *tmp = NULL ;
      RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , idx , 
			       tfold == true ? PlPl : NOFOLD , 
			       5 , 5 ) ;
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
read_tetra_corr( double ***X , // is the tiem
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold )
{
  *NSLICES = 2 * INPARAMS -> NFILES ; 

  //const size_t NCHANNELS = 4 ; // DiqDiq, Dimes, V , P
  const size_t NOPS = 2 ;
  const size_t opmap[ 2 ] = { 0 , 5 } ; 
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double ) ) ;
  *mom = malloc( *NSLICES * sizeof( struct mom_info ) ) ;

  // copy these right
  copy_inparams( INPARAMS , 1 , 3 ) ;
  copy_inparams( INPARAMS , 1 , 2 ) ;
  copy_inparams( INPARAMS , 0 , 1 ) ;
  copy_inparams( INPARAMS , 0 , 0 ) ;

  int nfiles = 0 , idx = 0 ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    char str[ 256 ] ;
    {
      sprintf( str , "%s" , INPARAMS -> traj_file[ nfiles ] ) ;
      size_t op , nop ;
      for( nop = 0 ; nop < NOPS ; nop++ ) {
	op = opmap[ nop ] ;
	struct mom_info *tmp = NULL ;
	{
	  RAW[idx] = read_rawcorr( INPARAMS , &tmp , str , nfiles , 
				   tfold == true ? PlPl : NOFOLD , 
				   op , 0 ) ;
	  (*mom)[0] = tmp[0] ;
	  (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[nfiles] * sizeof( double ) ) ;
	  
	  int ops[2] = { op , op } ;
	  average_similar( RAW , X , INPARAMS , mom , ops , str , idx , tfold ) ;
	  idx++ ;
	}
      }
      // loop of operators
    }
  }
  printf( "[Tetras] IO finished\n" ) ;
  return RAW ;
}
