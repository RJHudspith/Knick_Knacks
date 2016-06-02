/**
   IO for the gluon propagator ...
 */

#include "fitfunc.h"
#include <stdint.h>
#include "GLU_bswap.h"
#include "GLU_timer.h"

static bool must_swap = true ;

// 
struct mom_info
fill_mominfo( const int ND ,
	      const int mom[ ND ] ,
	      const int *dimensions ,
	      const momtype type )
{
  struct mom_info mominfo ;
  mominfo.p8 = mominfo.p6 = mominfo.p4 = mominfo.p2 = 0.0 ;

  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mominfo.n[ mu ] = mom[ mu ] ;
    const double p = 2. * M_PI * (double)mom[ mu ] / (double)dimensions[ mu ] ;
    switch( type ) {
    case TWOSIN_MOM :
      mominfo.p2 += pow( 2.0 * sin( 0.5 * p ) , 2 ) ;
      mominfo.p4 += pow( 2.0 * sin( 0.5 * p ) , 4 ) ;
      mominfo.p6 += pow( 2.0 * sin( 0.5 * p ) , 6 ) ;
      mominfo.p8 += pow( 2.0 * sin( 0.5 * p ) , 8 ) ;
      break ;
    case PSQ_MOM :
      mominfo.p2 += pow( p , 2 ) ;
      mominfo.p4 += pow( p , 4 ) ;
      mominfo.p6 += pow( p , 6 ) ;
      mominfo.p8 += pow( p , 8 ) ;
      break ;
    case SIN_MOM :
      mominfo.p2 += pow( sin( p ) , 2 ) ;
      mominfo.p4 += pow( sin( p ) , 4 ) ;
      mominfo.p6 += pow( sin( p ) , 6 ) ;
      mominfo.p8 += pow( sin( p ) , 8 ) ;
      break ;
    case RSQ :
      mominfo.p2 += pow( mom[mu] , 2 ) ;
      mominfo.p4 += pow( mom[mu] , 4 ) ;
      mominfo.p6 += pow( mom[mu] , 6 ) ;
      mominfo.p8 += pow( mom[mu] , 8 ) ;
      /*
      mominfo.p2 += pow( mom[mu]<0?dimensions[ mu ]+mom[mu]:mom[mu] , 2 ) ;
      mominfo.p4 += 0.0 ;
      mominfo.p6 += 0.0 ; 
      mominfo.p8 += 0.0 ;
      */
      break ;
    }
  }
  return mominfo ;
}

// I want a reader that takes the n_\mu's and
// computes the p^2-terms and writes them out
// as double precision numbers ...
struct mom_info *
read_momlist( int *NDATA ,
	      int *ND ,
	      FILE *file ,
	      const int *dimensions ,
	      const momtype mom_type ,
	      const bool dummy )
{
  /*
  uint32_t magic[1] ;
  must_swap = false ;
  if( fread( magic , (sizeof( uint32_t )) , 1 , file ) != 1 ) {
    printf( "Unexpected EOF ..\n" ) ;
    return NULL ;
  }
  if( magic[0] != 717685 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 717685 ) {
      printf( "Magic number failure %u , leaving \n" , magic[0] ) ;
      return NULL ;
    }
    must_swap = true ;
  }
  */
  // if we are a dummy run, just skip the momlist
  if( dummy == true ) {
    fseek ( file , *NDATA * ( *ND + 1 ) * sizeof( uint32_t ) + 4 , SEEK_CUR );
    return NULL ;
  } else {

    // read in the size ...
    uint32_t size[1] ;
    if( fread( size , (sizeof( uint32_t )) , 1 , file ) != 1 ) return NULL ;
    if( must_swap ) { bswap_32( 1 , size ) ; }
  
    *NDATA = size[0] ;

    // allocate the psq
    struct mom_info *mominfo = malloc( *NDATA * sizeof( struct mom_info ) ) ;

    int NDtmp = 0 ;
    int i ;
    for( i = 0 ; i < *NDATA ; i++ ) {
      // ok so the momentum list is interesting
      int size[ 1 ] ;
      if( fread( size , (sizeof( int )) , 1 , file ) != 1 ) return NULL ;
      if( must_swap ) { bswap_32( 1 , size ) ; }
      NDtmp = size[0] ;

      // and read in the first ND ints
      int mom[ NDtmp ] ;
      if( fread( mom , (sizeof( int )) , NDtmp , file ) != NDtmp ) return NULL ;
      if( must_swap ) { bswap_32( NDtmp , mom ) ; }
    
      // and compute the momenta
      mominfo[ i ] = fill_mominfo( NDtmp , mom , dimensions , mom_type ) ;
    }

    *ND = NDtmp ;
    // and that's all folks

    return mominfo ;
  }
}

// read in nmeas's data from the file
static void
read_data( struct resampled *sample , 
	   int *NDATA ,
	   FILE *file ,
	   const int meas )
{
  int size[1] ;
  if( fread( size , (sizeof( int )) , 1 , file ) != 1 ) exit(1) ;
  if( must_swap ) bswap_32( 1 , size ) ;

  if( size[0] != *NDATA ) {
    printf( "NDATA-SIZE mismatch %d %d \n" , size[0] , *NDATA ) ;
    exit(1) ;
  }

  // read in all the doubles
  double *tmp = malloc( *NDATA * sizeof( double ) ) ;
  if( fread( tmp , (sizeof( double )) , *NDATA , file ) != *NDATA ) exit(1) ;
  if( must_swap ) bswap_64( *NDATA , tmp ) ;
  int i ;
  for( i = 0 ; i < *NDATA ; i++ ) {
    sample[i].resampled[meas] = tmp[i] ;
  }

  free( tmp ) ;
}

// read a raw GLU propagator file
struct resampled*
read_rawGLU( struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int nfile ,
	     const momtype mom_type )
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

  int ND , ndata = 0 ;
  if( ( *mominfo = read_momlist( &ndata , &ND , file , 
				 INPARAMS -> dimensions[ nfile ] , 
				 mom_type ,
				 false ) ) == NULL ) {
    printf( "MOMlist failure \n" ) ;
    return NULL ;
  }
  printf( "Momlist read \n" ) ;
  print_time() ;

  // set this
  INPARAMS -> NDATA[ nfile ] = ndata ;

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
    int tND = ND , tndata = ndata ;
    struct mom_info *tmominfo = read_momlist( &tndata , &tND , loc_file , 
					      INPARAMS -> dimensions[ nfile ] , 
					      true , mom_type ) ;

    // set the raw data
    read_data( sample , &tndata , loc_file , meas ) ;

    free( tmominfo ) ;

    fclose( loc_file ) ;
  }

  return sample ;
}

struct resampled**
read_GLUprop( struct mom_info ***mominfo ,
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
    RAW[nfile] = read_rawGLU( *mominfo+nfile , INPARAMS , nfile , mom_type ) ;
    (*X)[nfile] = (double*)malloc( INPARAMS -> NDATA[nfile] * sizeof( double ) ) ;
    int j ;
    #pragma omp parallel for private(j)
    for( j = 0 ; j < INPARAMS -> NDATA[nfile] ; j++ ) {
      (*X)[nfile][j] = (*mominfo)[nfile][j].p2 ;
    }
  }
  *NSLICES = INPARAMS -> NFILES ;
  return RAW ;
}

////////////////////// SPURER DISTRIBUTION ////////////////////////////////////
static struct resampled*
read_distGLU( struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int nfile ,
	      const momtype mom_type )
{
  // read the raw data do it all in this ....
  FILE *file = fopen( INPARAMS -> traj_file[nfile] , "rb" ) ;
  if( file == NULL ) {
    printf( "Super distribution %s not found!\n" , INPARAMS -> traj_file[nfile] ) ;
    exit(1) ;
  }

  // read the momlist
  int ND ;
  *mominfo = read_momlist( &( INPARAMS -> NDATA[nfile] ) , &ND , file , 
			   INPARAMS -> dimensions[nfile]  , mom_type , false ) ;

  struct resampled *sample = malloc( INPARAMS -> NDATA[nfile] * sizeof( struct resampled ) ) ;
  
  // loop the number of momenta
  int k ;
  for( k = 0 ; k < INPARAMS -> NDATA[nfile] ; k++ ) {
    uint32_t nsamp[1] = {} ;
    if( fread( nsamp , sizeof( uint32_t ) , 1 , file ) != 1 ) {
      printf( "NSAMP misread ... Leaving \n" ) ;
      return NULL ;
    }
    if( must_swap ) bswap_32( 1 , nsamp ) ;

    // allocate nsamples
    sample[k].resampled = malloc( nsamp[0] * sizeof( double ) ) ;
    sample[k].NSAMPLES = nsamp[0] ;

    // read in the data
    if( fread( sample[k].resampled , sizeof( double ) , nsamp[0] , file ) != nsamp[0] ) {
      printf( "samples %d misread ... Leaving \n" , k ) ;
      return NULL ;
    }
    if( must_swap ) bswap_64( nsamp[0] , sample[k].resampled ) ;
  }

  /*
  // is right at the bottom
  uint32_t resample[1] ;
  if( fread( resample , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    printf( "read failure \n" ) ;
    return NULL ;
  }
  if( must_swap ) bswap_32( 1 , resample ) ;
  */
  for( k = 0 ; k < INPARAMS -> NDATA[nfile] ; k++ ) {
    sample[k].restype = RAWDATA ; //(resample_type)resample[0] ;
    compute_err( &sample[k] ) ;
  }

  return sample ;
}

// read a super distribution
struct resampled**
read_superdist( struct mom_info ***mominfo ,
		double ***X , // is the momentum^2
		struct input_params *INPARAMS ,
		int *NSLICES ,
		const momtype mom_type )
{
  struct resampled **RAW = malloc( INPARAMS -> NFILES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( INPARAMS -> NFILES * sizeof( double ) ) ;
  *mominfo = (struct mom_info**)malloc( INPARAMS -> NFILES * sizeof( struct mom_info* ) ) ;

  int i ;
  for( i = 0 ; i < INPARAMS -> NFILES ; i++ ) {
    RAW[i] = read_distGLU( *mominfo+i , INPARAMS , i , mom_type ) ;
    (*X)[i] = (double*)malloc( INPARAMS -> NDATA[i] * sizeof( double ) ) ;
    int j ;
    #pragma omp parallel for private(j)
    for( j = 0 ; j < INPARAMS -> NDATA[i] ; j++ ) {
      (*X)[i][j] = (*mominfo)[i][j].p2 ;
    }
  }

  *NSLICES = INPARAMS -> NFILES ;
  return RAW ;
}
