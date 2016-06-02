/**
   @brief read a correlator from our contraction code
 */

#include "fitfunc.h"

#include <stdint.h>
#include <complex.h>

#include "GLUdata_glueprop.h"
#include "GLU_bswap.h"
#include "GLU_timer.h"

static bool must_swap = false ;

static int L0 ;

// read the magics, time and ngammas
static int
read_magic_gammas( struct mom_info **mom , 
		   FILE *file , 
		   uint32_t NGSRC[1] ,
		   uint32_t NGSNK[1] ,
		   uint32_t LT[1] )
{
  uint32_t magic[1] = {} ;
  if( fread( magic , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    return FAILURE ;
  }
  // check the magic number, tells us the edianness
  if( magic[0] != 67798233 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 67798233 ) {
      printf( "Magic number read failure\n" ) ;
      return FAILURE ;
    }
    must_swap = true ;
  }

  // read the momentum list
  uint32_t NMOM[ 1 ] ;
  if( fread( NMOM , sizeof(uint32_t) , 1 , file ) == FAILURE ) return FAILURE ;
  if( must_swap ) bswap_32( 1 , NMOM ) ;
  
  *mom = malloc( NMOM[0] * sizeof( struct mom_info ) ) ;

  size_t p , mu ;
  for( p = 0 ; p < NMOM[0] ; p++ ) {
    int n[ 4 ] ;
    if( fread( n , sizeof( uint32_t ) , 4 , file ) == FAILURE ) return FAILURE ;
    if( n[ 0 ] != 4-1 ) {
      printf( "[MOMLIST] %d should be %d \n" , n[ 0 ] , 3 ) ;
      return FAILURE ;
    }
    for( mu = 0 ; mu < 3 ; mu++ ) {
      mom[ p ] -> n[ mu ] = (int)n[ 1 + mu ] ;
    }
  }

  if( fread( NMOM , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , NMOM ) ;

  if( fread( NGSRC , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , NGSRC ) ;

  if( fread( NGSNK , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , NGSNK ) ;

  if( fread( LT , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , LT ) ;

  // set this
  L0 = LT[0] ;
  return SUCCESS ;
}

void
copy_inparams( struct input_params *INPARAMS ,
	       const int idx1 ,
	       const int idx2 )
{
  if( idx1 > idx2 ) { 
    printf( "Must copy idx1 into idx2\n" ) ;
    exit(1) ;
  }
  INPARAMS -> NDATA[ idx2 ] = INPARAMS -> NDATA[ idx1 ] ;
  INPARAMS -> NRAW[ idx2 ] = INPARAMS -> NRAW[ idx1 ] ;
  INPARAMS -> traj_begin[ idx2 ] = INPARAMS -> traj_begin[ idx1 ] ;
  INPARAMS -> traj_end[ idx2 ] = INPARAMS -> traj_end[ idx1 ] ;
  INPARAMS -> traj_increment[ idx2 ] = INPARAMS -> traj_increment[ idx1 ] ;
  INPARAMS -> quarks[ idx2 ] = INPARAMS -> quarks[ idx1 ] ;
  int j ;
  for( j = 0 ; j < MAX_NFILES ; j++ ) {
    INPARAMS -> dimensions[ idx2 ][j] = INPARAMS -> dimensions[ idx1 ][j] ;
  }
}

// if we have multiple channels we expand our input params
void
expand_inparams( struct input_params *INPARAMS ,
		 const int NCHANNELS )
{
  int nfiles , idx = 0 ;
  for( nfiles = INPARAMS->NFILES-1 ; nfiles > -1 ; nfiles-- ) {
    for( idx = 0 ; idx < NCHANNELS ; idx++ ) {
      copy_inparams( INPARAMS , nfiles , idx + NCHANNELS*nfiles ) ;
    }
  }
  return ;
}

static int
read_corrfile( struct resampled *sample ,
	       struct mom_info **mominfo ,
	       FILE *file ,
	       const int meas ,
	       const int src , 
	       const int snk ,
	       const foldselection fold )
{
  struct mom_info *mom ;
  uint32_t NGSRC[1] , NGSNK[1] , LT[1] ;
  if( read_magic_gammas( &mom , file , NGSRC , NGSNK , LT ) == FAILURE ) {
    return FAILURE ;
  }

  const int gseek = (int)( snk + NGSNK[0] * src ) ;
  const int cseek = gseek * (int)( LT[0] * sizeof( double complex ) ) ;
  const int Lseek = ( gseek - 1 ) * (int)( sizeof( uint32_t ) ) ;
  fseek( file , Lseek + cseek , SEEK_CUR ) ;

  if( fread( LT , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    printf( "LT :: %d \n" , LT[0] ) ;
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , LT ) ;
  if( (int)LT[0] != L0 ) { 
    printf( "LT Read failure %d %d \n" , (int)LT[0] , L0 ) ; 
    return FAILURE ; 
  }

  double complex C[ L0 ] ;
  if( fread( C , sizeof( double complex ) , L0 , file ) != L0 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_64( 2 * L0 , C ) ;
  int t ;
  switch( fold ) {
  case NOFOLD :
    for( t = 0 ; t < L0 ; t++ ) {
      sample[t].resampled[meas] = creal( C[t] ) ;
      //sample[t].resampled[meas] = cimag( C[t] ) ;
    }
    break ;
  case PlPl :
    sample[0].resampled[meas] = fabs( creal( C[0] ) ) ;
    for( t = 1 ; t < L0/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( C[t] + C[L0-t] ) ;
    }
    break ;
  case PlMi :
    sample[0].resampled[meas] = creal( C[0] ) ;
    for( t = 1 ; t < L0/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( C[t] - C[L0-t] ) ;
    }
    break ;
  case MiPl :
    sample[0].resampled[meas] = -creal( C[0] ) ;
    for( t = 1 ; t < L0/2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( -C[t] + C[L0-t] ) ;
    }
    break ;
  case MiMi : 
    sample[0].resampled[meas] = -creal( C[0] ) ;
    for( t = 1 ; t < L0/2 ; t++ ) {
      sample[t].resampled[meas] = -0.5 * creal( C[t] + C[L0-t] ) ;
    }
    break ;
  }
  return SUCCESS ;
}

// read a raw correlator
struct resampled *
read_rawcorr( struct input_params *INPARAMS ,
	      struct mom_info **mom ,
	      const char *filename ,
	      const int nfile ,
	      const foldselection fold ,
	      const int src , 
	      const int snk ) 
{
  const int NMEAS = ( INPARAMS -> traj_end[nfile] - INPARAMS -> traj_begin[nfile] ) / 
    INPARAMS -> traj_increment[nfile] ;

  printf( "BUTT :: %d %d %d \n" , INPARAMS -> traj_end[nfile] , INPARAMS -> traj_begin[nfile] , INPARAMS -> traj_increment[nfile] ) ;
  printf( "IN DIS :: NMEAS :: %d nfiles :: %d \n" , NMEAS , nfile ) ;

  char filestr[ 256 ] ;
  sprintf( filestr , filename , INPARAMS -> traj_begin[nfile] ) ;

  // file open
  FILE *file = fopen( filestr , "rb" ) ;
  if( file == NULL ) {
    printf( "Cannot open %s \n" , filestr ) ;
    return NULL ;
  }
  print_time() ;

  uint32_t NGSRC[1] , NGSNK[1] , LT[1] = { 0 } ;
  if( read_magic_gammas( mom , file , NGSRC , NGSNK , LT ) == FAILURE ) {
    free( mom ) ;
    return NULL ;
  }
  printf( "Magic\n" ) ;

  // set this
  if( fold == NOFOLD ) {
    INPARAMS -> NDATA[ nfile ] = (int)LT[0] ;
  } else {
    INPARAMS -> NDATA[ nfile ] = (int)LT[0]/2 ;
  }

  fclose( file ) ;

  // allocate the resampled struct
  struct resampled *sample = malloc( INPARAMS -> NDATA[ nfile ]  * sizeof( struct resampled ) ) ;
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
    FILE *loc_file = fopen( loc_filestr , "rb" ) ;
    if( loc_file == NULL ) {
      printf( "Cannot open %s \n" , loc_filestr ) ;
      failure = true ;
    }

    // set the raw data
    struct mom_info *tmp = NULL ;
    read_corrfile( sample , &tmp , loc_file , meas , src , snk , fold ) ;
    free( tmp ) ;

    fclose( loc_file ) ;
  }

  return ( failure == false ) ? sample : NULL ;
}

// quick corr access
struct resampled *
corrs( double ***X , 
       struct input_params *INPARAMS ,
       struct mom_info **mom ,
       const char *filename ,
       const foldselection fold ,
       const int src , 
       const int snk ,
       const int idx )
{
  struct resampled *RAW ;
  RAW = read_rawcorr( INPARAMS , mom , filename , idx , fold , src , snk ) ;
  (*X)[idx] = (double*)malloc( INPARAMS -> NDATA[ idx ] * sizeof( double ) ) ;
  int j ;
  for( j = 0 ; j < INPARAMS -> NDATA[ idx ] ; j++ ) { (*X)[idx][j] = (double)j ; }
  return RAW ;
}
