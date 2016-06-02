/**
   @file readflow_couple.c
   @brief read a plain text file for computing the gradient flow coupling
 */
#include "fitfunc.h"
#include <stdio.h>

static int
get_ndata( FILE *file )
{
  int i = 0 ;
  double t , Gt ;
  while( fscanf( file , "%lf %lf \n" , &t , &Gt ) != EOF ) {
    //printf( "%d %f %f \n" , i , t , Gt ) ;
    i++ ;
  }
  rewind( file ) ;
  return i ;
}

// read a single file
static int
read_data( struct resampled *sample ,
	   double *x ,
	   FILE *file ,
	   const int NDATA ,
	   const int meas ,
	   const double ainverse )
{
  const double alpha_fact = 4.0 * M_PI / 3.0 ;
  // expects to be 
  // time G(t)
  size_t i = 0 ;
  double t , Gt ;
  while( fscanf( file , "%lf %lf \n" , &t , &Gt ) != EOF ) {
    sample[i].resampled[meas] = alpha_fact * Gt ;
    x[i] = ainverse / ( sqrt( 8. * t ) ) ;
    //printf( "%f %f \n" , t , Gt ) ;
    i++ ;
  }
  //printf( "%zu %d \n" , i , NDATA ) ; 
  return i == NDATA ? SUCCESS : FAILURE ;
}

// read a raw GLU propagator file
struct resampled*
read_rawflow( double **x ,
	      struct input_params *INPARAMS ,
	      const int nfile )
{
  const int NMEAS = ( INPARAMS -> traj_end[nfile] - INPARAMS -> traj_begin[nfile] ) / 
    INPARAMS -> traj_increment[nfile] ;

  char filestr[ 512 ] ;
  sprintf( filestr , INPARAMS -> traj_file[nfile] , INPARAMS -> traj_begin[nfile] ) ;

  // read initial momlist
  FILE *file = fopen( filestr , "r" ) ;
  if( file == NULL ) {
    printf( "Cannot open %s \n" , filestr ) ;
    return NULL ;
  }

  // get the data
  const int NDATA = get_ndata( file ) ;
  printf( "NDATA :: %d \n" , NDATA ) ;
  INPARAMS -> NDATA[ nfile ] = NDATA ;

  *x = malloc( NDATA * sizeof( double ) ) ;

  // allocate the resampled struct
  struct resampled *sample = malloc( NDATA * sizeof( struct resampled ) ) ;
  int i ;
  for( i = 0 ; i < NDATA ; i++ ) {
    sample[i].resampled = malloc( NMEAS * sizeof( double ) ) ;
    sample[i].restype = RAWDATA ;
    sample[i].NSAMPLES = NMEAS ;
  }

  // loop trajectory files
  bool failure = false ;
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
      failure = true ;
    }

    // set the raw data
    if( read_data( sample , *x , loc_file , NDATA , meas ,
		   INPARAMS -> quarks[ nfile ].ainverse ) == FAILURE ) {
      printf( "Read data failure %d \n" , i ) ;
      failure = true ;
    }

    // close the local file
    fclose( loc_file ) ;
  }
  if( failure == true ) return NULL ;

  return sample ;
}

struct resampled**
read_flowcouple( double ***X ,
		 struct input_params *INPARAMS ,
		 int *NSLICES )
{
  const size_t NDATA = INPARAMS -> NFILES ;
  struct resampled **RAW = malloc( NDATA * sizeof( struct resampled* ) ) ;
  *X = (double**)malloc( NDATA * sizeof( double ) ) ;
  int nfile = 0 ;
  for( nfile = 0 ; nfile < INPARAMS -> NFILES ; nfile++ ) {
    RAW[nfile] = read_rawflow( *X+nfile , INPARAMS , nfile ) ;
  }
  *NSLICES = NDATA ;
  return RAW ;
}
