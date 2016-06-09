/**
   @file  read_distribution.c
   @brief reader for the distribution type in UKhadron
 */

#include "fitfunc.h"
#include "GLU_bswap.h" // for the byte swaps
#include "stats.h"     // for the errors
#include "GLUdata_glueprop.h"

// read in an array of our resampled structure RES
static int
read_doubles( FILE *file , struct resampled *RES , const int NSAMPLES )
{
  if( fread( RES -> resampled , sizeof(double) , NSAMPLES , file ) != NSAMPLES ) return FAILURE ;
  if( !BigEndian ) bswap_64( NSAMPLES , RES-> resampled ) ;
  double average[1] ;
  if( fread( average , sizeof(double) , 1 , file ) != 1 ) {
    printf( "Average read failure \n" ) ;
    return FAILURE ;
  }
  if( !BigEndian ) bswap_64( 1 , average ) ;
  const double compar = average[ 0 ] ;
  compute_err( RES ) ;
  if( fabs( compar - RES -> avg ) > 1E-15 ) {
    //printf( "AVGs not same (Read Computed) !! %1.12e %1.12e %1.12e \n" , compar , RES -> avg , compar - RES -> avg ) ;
  } else {
    //printf( "AVG :: %1.12e \n" , RES -> avg ) ;
  }
  RES -> avg = compar ;
  printf( "AVG :: %1.12e \n" , RES -> avg ) ;
  return SUCCESS ;
}

// parse the two indicators for this format ...
static int
parse_sizes( FILE *file , struct resampled *RES , const int NSAMPLES , const resample_type res_check )
{
  int vals[ 2 ] ;
  if( fread( vals , sizeof(int) , 2 , file ) != 2 ) return FAILURE ;
  if( !BigEndian ) bswap_32( 2 , vals ) ;
  if( vals[0] != res_check ) { 
    printf( "Different sample type %d, expected %d Leaving \n" ,
	    vals[1] , res_check ) ;
    return FAILURE ;
  } else if( vals[1] != NSAMPLES ) {
    printf( "Different NSAMPLES %d, expected %d , Leaving \n" ,
	    vals[1] , NSAMPLES ) ;
    return FAILURE ;
  }
  RES -> restype = vals[0] ;
  RES -> NSAMPLES = vals[1] ;
  return SUCCESS ;
}

int
read_single_distribution( FILE *file , 
			  struct resampled *RES ,
			  struct mom_info *moms ,
			  const int NSAMPLES ,
			  const resample_type res_check ,
			  const int *dimensions ,
			  const momtype type )
{
  // read momenta
  uint32_t num_mom[1] ;
  if( fread( num_mom , sizeof(uint32_t) , 1 , file ) != 1 ) {
    printf( "[IO] num_mom reading failure %u\n" , num_mom[0] ) ;
    return FAILURE ;
  }
  if( !BigEndian ) bswap_32( 1 , num_mom ) ;
  // read a single momentum
  uint32_t mom[4] ;
  if( fread( mom , sizeof(uint32_t) , 4 , file ) != 4 ) { 
    printf( "[IO] mom reading failure \n" ) ;
    return FAILURE ;
  }
  if( !BigEndian ) bswap_32( 4 , mom ) ;
  printf( "MOM :: %d %d %d %d \n" , mom[0] , mom[1] , mom[2] , mom[3] ) ;
  moms -> n[0] = mom[1] ; moms -> n[1] = mom[2] ; moms -> n[2] = mom[3] ; 
  moms -> n[3] = 0 ;

  *moms = fill_mominfo( 4, moms->n , dimensions , type ) ;

  printf( "MOMS filled \n" ) ;

  if( parse_sizes( file , RES , NSAMPLES , res_check ) == FAILURE ) return FAILURE ;
  if( read_doubles( file , RES , NSAMPLES ) == FAILURE ) return FAILURE ;
  return SUCCESS ;
}

struct resampled **
read_dispersions( double ***X ,
		  struct input_params *INPARAMS ,
		  struct mom_info ***mominfo ,
		  int *NSLICES )
{
  *NSLICES = 1 ;
  struct resampled **RAW = malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = malloc( *NSLICES * sizeof( double* ) ) ;
  *mominfo = malloc( *NSLICES * sizeof( struct mom_info* ) ) ;


  RAW[0] = malloc( INPARAMS -> NFILES * sizeof( struct resampled ) ) ;
  (*X)[0] = malloc( INPARAMS -> NFILES * sizeof( double* ) ) ; 
  (*mominfo)[0] = malloc( INPARAMS -> NFILES * sizeof( struct mom_info ) ) ; 

  // set this
  INPARAMS -> NDATA[0] = INPARAMS -> NFILES ;

  // loop files
  int nfiles ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {

    // open one of the perhaps many files
    FILE *file = fopen( INPARAMS -> traj_file[ nfiles ] , "rb" ) ;
    if( file == NULL ) { 
      printf( "%s (file %d) not available ... leaving \n" , 
	      INPARAMS -> traj_file[nfiles] , nfiles ) ;
      exit(-1) ;
    }

    RAW[0][nfiles].resampled = malloc( INPARAMS -> NBOOTS * sizeof( double ) ) ;

    if( read_single_distribution( file , &RAW[0][nfiles] , 
				  &(*mominfo)[0][nfiles] , 
				  INPARAMS -> NBOOTS , 
				  INPARAMS -> resample ,
				  INPARAMS -> dimensions[nfiles] ,
				  INPARAMS -> mom_type ) == FAILURE ) {
      return NULL ;
    }

    // could do with a call to set_mominfo here ...
  }

  // set the xs
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    printf( "In this :: %e\n" , (*mominfo)[0][nfiles].p2 ) ;
    (*X)[0][nfiles] = (*mominfo)[0][nfiles].p2 ;
  }

  return RAW ;
}

/*
  OK, so UKhadron's output is ..

  Size of the distribution vector.size() 
  Distribution Type :: I think their enum is the same as mine
  Size of the number of resampled values Distribution.Nmeas() == vector.size()
  And then the average is at the end of this
  
  And then, back to the beginning

  expect two distributions back to back
*/
struct resampled *
read_distribution( int *NDATA ,
		   FILE *file )
{
  // look up the details ...
  /*
  int deets[ 3 ] ;
  if( fread( deets , sizeof(int) , 3 , file ) != 3 ) exit(-1) ;
  if( !BigEndian ) bswap_32( 3 , deets ) ;
  printf( "Array length %d\n" , deets[0] ) ;
  printf( "Resample type %d\n" , deets[1] ) ;
  printf( "Nsamples %d \n" , deets[2] ) ;
  // allocations
  struct resampled *res = malloc( deets[0] * sizeof( struct resampled ) ) ;
  int i ;
  for( i = 0 ; i < deets[0] ; i++ ) {
    res[i].resampled = malloc( deets[2] * sizeof( double ) ) ;
    res[i].restype   = deets[1] ;
    res[i].NSAMPLES  = deets[2] ;
  }
  *NDATA = deets[0] ;
  // and read in the data first one is special
  read_doubles( file , &res[0] , deets[2] ) ;
  int j ;
  for( j = 1 ; j < deets[0] ; j++ ) {
    parse_sizes( file , &res[j] , deets[2] , deets[1] ) ;
    read_doubles( file , &res[j] , deets[2] ) ;
  }
  */

  int deets[ 2 ] ;
  if( fread( deets , sizeof(int) , 2 , file ) != 2 ) exit(-1) ;
  if( !BigEndian ) bswap_32( 2 , deets ) ;
  //printf( "Resample type %d\n" , deets[0] ) ;
  //printf( "Nsamples %d \n" , deets[1] ) ;

  // allocations
  struct resampled *res = malloc( sizeof( struct resampled ) ) ;
  int i ;
  for( i = 0 ; i < 1 ; i++ ) {
    res[i].resampled = malloc( deets[1] * sizeof( double ) ) ;
    res[i].restype   = deets[0] ;
    res[i].NSAMPLES  = deets[1] ;
  }
  *NDATA = 1 ;

  // and read in the data first one is special
  read_doubles( file , &res[0] , deets[1] ) ;

  return res ;
}

// have a function that says whether we can flatten
// a distribution into an array of double
int
can_flatten( double *array , struct resampled *res , const int N )
{
  int i ;
  for( i = 0 ; i < N ; i++ ) {
    array[i] = res[i].avg ;
    int j ;
    for( j = 0 ; j < res[i].NSAMPLES ; j++ ) {
      if( fabs( array[i] - res[i].resampled[j] ) > 1E-15 ) {
	printf( "Cannot be flattened %e %e %1.8e\n" , array[i] , res[i].resampled[j] , 
		array[i] - res[i].resampled[j] ) ;
	return FAILURE ;
      }
    }
  }
  return SUCCESS ;
}

// read in the distribution
struct resampled **
read_UKhadron( double ***X ,
	       struct input_params *INPARAMS ,
	       int *NSLICES )
{
  struct resampled **RAW = (struct resampled**)malloc( INPARAMS -> NFILES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( INPARAMS -> NFILES * sizeof( double* ) ) ;

  int nfiles ;
  for( nfiles = 0 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {

    // open one of the perhaps many files
    FILE *file = fopen( INPARAMS -> traj_file[nfiles] , "rb" ) ;
    if( file == NULL ) { 
      printf( "%s (file %d) not available ... leaving \n" , 
	      INPARAMS -> traj_file[nfiles] , nfiles ) ;
      exit(-1) ;
    }
    // try and read any x-information
    struct resampled *tmp = read_distribution( &INPARAMS->NDATA[nfiles] , file ) ;
    (*X)[ nfiles ] = malloc( INPARAMS->NDATA[nfiles] * sizeof( double ) ) ;

    // can we pack it into a double array?
    if( can_flatten( (*X)[ nfiles ] , tmp , INPARAMS->NDATA[nfiles] ) == FAILURE ) {
      exit(-1) ;
    }

    // read in the distribution
    int testNDATA ;
    RAW[nfiles] = read_distribution( &testNDATA , file ) ;
    if( testNDATA != INPARAMS->NDATA[nfiles] ) { 
      exit(-1) ;
    }

    // data frees
    int k ;
    for( k = 0 ; k < INPARAMS->NDATA[nfiles] ; k++ ) {
      free( tmp[k].resampled ) ;
    }
    free( tmp ) ;
    fclose( file ) ;
  }
  *NSLICES = INPARAMS->NFILES ;
  return RAW ;
}


// read in the distribution
struct resampled **
read_SUSY( double ***X ,
	   struct input_params *INPARAMS ,
	   int *NSLICES )
{
  // look at a^2's
  int a2counter[ INPARAMS -> NFILES ] ;
  double ainvref = INPARAMS -> quarks[0].ainverse ;
  int nfiles , diff = 0 , n = 1 ;
  for( nfiles = 1 ; nfiles < INPARAMS -> NFILES ; nfiles++ ) {
    printf( "%f \n" , INPARAMS -> quarks[nfiles].ainverse ) ;
    if( INPARAMS -> quarks[nfiles].ainverse != ainvref ) {
      ainvref = INPARAMS -> quarks[nfiles].ainverse ;
      a2counter[ diff ] = n ;
      printf( "in here %d %d \n" , diff , a2counter[diff] ) ;
      n = 1 ;
      diff++ ;
    } else {
      n++ ;
    }
  }
  if( diff == 0 ) {
    a2counter[ diff ] = INPARAMS -> NFILES ;
  } else {
    a2counter[ diff ] = n ; //INPARAMS -> NFILES - a2counter[ diff - 1 ] ;
  }
  diff++ ;
  *NSLICES = diff ;

  struct resampled **RAW = (struct resampled**)malloc( *NSLICES * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( *NSLICES * sizeof( double* ) ) ;

  n = 0 ;
  for( nfiles = 0 ; nfiles < diff ; nfiles++ ) {

    RAW[ nfiles ] = malloc( a2counter[ nfiles ] * sizeof( struct resampled ) ) ;
    (*X)[ nfiles ] = malloc( a2counter[ nfiles ] * sizeof( double ) ) ;

    size_t i ;
    for( i = 0 ; i < a2counter[ nfiles ] ; i++ ) {

      // read the file
      FILE *file = fopen( INPARAMS -> traj_file[ n ] , "rb" ) ;
      if( file == NULL ) { 
	printf( "%s (file %d) not available ... leaving \n" , 
		INPARAMS -> traj_file[n] , n ) ;
	exit(-1) ;
      }

      // try and read any x-information
      struct resampled *tmp = read_distribution( &INPARAMS->NDATA[ nfiles ] , file ) ;

      // can we pack it into a double array?
      switch( INPARAMS->type ) {
      case KK_SS_EXTRAP :
	(*X)[nfiles][i] = ( INPARAMS -> quarks[ n ].ms - INPARAMS -> quarks[ n ].ZV ) ;
	printf( "tmp1 ---> " ) ;
	int testNDATA ;
	tmp = read_distribution( &testNDATA , file ) ;
	printf( "tmp2 \n" ) ;
	break ;
#ifdef EPSILON_FIT
      case KK_ML_EXTRAP :
	(*X)[nfiles][i] = INPARAMS -> quarks[ n ].m_pi * INPARAMS -> quarks[ n ].m_pi / 
	  ( INPARAMS -> quarks[ n ].f_pi * INPARAMS -> quarks[ n ].f_pi ) ;
	break ;
#else
      case KK_ML_EXTRAP :
	(*X)[nfiles][i] = INPARAMS -> quarks[ n ].m_pi * INPARAMS -> quarks[ n ].m_pi / ( 0.12229 * 0.12229 ) ;
	break ;
#endif
      default :
	return NULL ;
      }

      // set the raw data
      RAW[nfiles][i].resampled = malloc( tmp[0].NSAMPLES * sizeof( double ) ) ;
      equate( &RAW[nfiles][i] , tmp[0] ) ;

      // data frees
      int k ;
      for( k = 0 ; k < INPARAMS->NDATA[nfiles] ; k++ ) {
	free( tmp[k].resampled ) ;
      }
      free( tmp ) ;
      fclose( file ) ;

      n++ ;
    }
    INPARAMS -> NDATA[ nfiles ] = a2counter[ nfiles ] ;
    printf( "CHECK :: %d %d \n" , INPARAMS->NDATA[ nfiles ] , n ) ;
  }
  return RAW ;
}

