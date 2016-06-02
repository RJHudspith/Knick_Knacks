#include "fitfunc.h"
#include "Utils.h"

void
print_moms( const int ND , const int n[ ND ] ) 
{
  printf( "( %d %d %d %d ) \n" ,
	  n[0] , n[1] , n[2] , n[3] ) ;
  return ;
}

// standard insertion sort
static void
insertion_sort( const int N , int list[ N ] )
{
  int i , hole ;
  // standard insertion sort
  for( i = 1 ; i < N ; i++ ) {
    const int insert = list[i] ;
    hole = i ;
    while( hole > 0 && insert < list[hole - 1] ) {
      list[hole] = list[hole-1] ;
      hole-- ;
    }
    list[hole] = insert ;
  }
  return ;
}

// test for Bose-symmetry within a triplet
static bool
are_ZN_equivalent( const int *n , 
		   const int *np ,
		   const int DIMS )
{
  int mu ;
  // excess elements have to be equivalent
  for( mu = DIMS ; mu < 4 ; mu++ ) {
    if( abs( n[ mu ] ) != abs( np[ mu ] ) ) {
      //if( n[ mu ] != np[ mu ] ) {
      #ifdef verbose
      printf( "NOT INLIST %d %d %d %d == %d %d %d %d\n" , 
	      n[0] , n[1] , n[2] , n[3] ,
	      np[0] , np[1] , np[2] , np[3] ) ;
      #endif
      return false ;
    }
  }

  // set up temporaries
  int tn[ DIMS ] ;
  int tnp[ DIMS ] ;
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    tn[ mu ] = abs( n[ mu ] ) ;
    tnp[ mu ] = abs( np[ mu ] ) ;
  }

  // sort the temporaries
  insertion_sort( DIMS , tn ) ;
  insertion_sort( DIMS , tnp ) ;

  for( mu = 0 ; mu < DIMS ; mu++ ) {
    if( tn[ mu ] != tnp[ mu ] ) {
      #ifdef verbose
      printf( "NOT INLIST %d %d %d %d == %d %d %d %d\n" , 
	      n[0] , n[1] , n[2] , n[3] ,
	      np[0] , np[1] , np[2] , np[3] ) ;
      #endif
      return false ;
    }
  }

  #ifdef verbose
  printf( "INLIST %d %d %d %d == %d %d %d %d\n" , 
	  n[0] , n[1] , n[2] , n[3] ,
	  np[0] , np[1] , np[2] , np[3] ) ;
  #endif

  return true ;
}

// count dat, look for file
static int
count_equivalents( const struct mom_info *mominfo , 
		   const int N ,
		   const int *dims ,
		   const int ND ) 
{
  int *count = calloc( N , sizeof(int) ) ;
  char str[ 256 ] ;
  sprintf( str , "/home/jamie/PHYSICS/FITTER/LOCAL/ZND_%d_%dx%dx%dx%d" , N ,
	   dims[0] , dims[1] , dims[2] , dims[3] ) ;
  FILE *file = fopen( str , "rw" ) ;
  int i , j , NEW[1] = { 0 } ;
  if( file == NULL ) {
    // because this is c we have to do a two pass-variant
    // first pass is to count the number of equivalent X's
    for( i = 0 ; i < N ; i++ ) {
      if( count[i] == 0 ){
	// NEW is like the "equivalents index"
	NEW[0]++ ;
	// loop list to see if we have a hit
	for( j = i + 1 ; j < N ; j++ ) {
	  if( fabs( mominfo[i].p2 - mominfo[j].p2 ) > 1E-12 ) break ;
	  if( are_ZN_equivalent( mominfo[i].n , mominfo[j].n , ND-1 ) == true ) {
	    count[j] = NEW[0] ;
	  }
	}
      }
    }
    FILE *file2 = fopen( str , "wb" ) ;
    fwrite( NEW , sizeof(int) , 1 , file2 ) ;
    fclose( file2 ) ;
  } else {
    if( fread( NEW , sizeof(int) , 1 , file ) != 1 ) {
      return FAILURE ;
    }
    fclose( file ) ;
  }
  printf( "[ZN AVERAGE] FROM :: %d to %d \n" , N , NEW[0] ) ;
  free(count) ;
  return NEW[0] ;
}

// this one is for the doubles
static void
average_equivalents( struct mom_info *mominfoavg ,
		     const struct mom_info *mominfo ,
		     struct resampled *RAWavg ,
		     const struct resampled *RAW ,
		     const int N ,
		     const int ND , 
		     const int NEW )
{
  int *count = calloc( N , sizeof(int) ) ;
  int i , j , NN = 0 , hit = 0 ;
  for( i = 0 ; i < N ; i++ ) {

    if( count[i] == 0 ){

      struct resampled tmp ;
      tmp.resampled = malloc( RAW[i].NSAMPLES * sizeof( double ) ) ;
      equate( &tmp , RAW[i] ) ;

      // copy the index over
      memcpy( &mominfoavg[ hit ] , &mominfo[ i ] , sizeof( struct mom_info ) ) ;

      NN = 1 ;
      // loop list to see if we have a hit
      for( j = i + 1 ; j < N ; j++ ) {
	if( fabs( mominfo[i].p2 - mominfo[j].p2 ) > 1E-12 ) break ;
	if( are_ZN_equivalent( mominfo[i].n , mominfo[j].n , ND-1 ) == true ) {

	  count[ j ] = i + 1 ;
	  add( &tmp , RAW[ j ] ) ;
	  NN++ ;
	}
      }
      
      divide_constant( &tmp , (double)NN ) ;
      equate( &RAWavg[hit] , tmp ) ;
      hit++ ; // increment the hit index

      free( tmp.resampled ) ;
    }
  }
  free(count) ;
  return ;
}

// overwrites INPARAMS, equates Z_N Bose-symmetry invariants
struct resampled **
ZN_average_slow( struct mom_info ***mominfoavg ,
		 struct input_params *INPARAMS ,
		 const struct mom_info **mominfo , 
		 const struct resampled **BOOTS ,
		 const int NSLICES ,
		 const bool actually_doit ,
		 const int Z_N )
{
  struct resampled **bootavg = ( struct resampled** )malloc( NSLICES * sizeof( struct resampled* ) ) ;
  *mominfoavg = ( struct mom_info** )malloc( NSLICES * sizeof( struct mom_info* ) ) ;
  if( actually_doit ) {
    // very slow looped version, uses the fact we are sorted a little
    // so we can stop looking in p^2 bins that are not the one we like
    
    int j ;
    //#pragma omp parallel for private(j)
    for( j = 0 ; j < NSLICES ; j++ ) {

      // does it by file
      const int NEW = count_equivalents( mominfo[j] , INPARAMS -> NDATA[j] , 
					 INPARAMS -> dimensions[j] , 4 ) ;
      if( NEW == FAILURE ) {
	return NULL ;
      }

      // allocate new momenta and bootstraps
      (*mominfoavg)[j] = (struct mom_info*)malloc( NEW * sizeof( struct mom_info ) ) ;
      bootavg[j] = (struct resampled*)malloc( NEW * sizeof( struct resampled ) ) ;
      int k ;
      for( k = 0 ; k < NEW ; k++ ) {
	bootavg[j][k].resampled = (double*)malloc( BOOTS[j][0].NSAMPLES * sizeof( double ) ) ;
      }

      // average them
      average_equivalents( (*mominfoavg)[j] , mominfo[j] ,
			   bootavg[j] , BOOTS[j] ,
			   INPARAMS -> NDATA[j] ,
			   4 , NEW ) ;

      INPARAMS -> NDATA[j] = NEW ;
    }
    // 
  } else {
    // do almost nothing
    int j ;
    for( j = 0 ; j < NSLICES ; j++ ) {
      (*mominfoavg)[j] = ( struct mom_info* )malloc( INPARAMS -> NDATA[j] * sizeof( struct mom_info ) ) ;
      bootavg[j] = ( struct resampled* )malloc( INPARAMS -> NDATA[j] * sizeof( struct resampled ) ) ;
      int k ;
      for( k = 0 ; k < INPARAMS -> NDATA[j] ; k++ ) {
	bootavg[j][k].resampled = (double*)malloc( BOOTS[j][0].NSAMPLES * sizeof( double ) ) ;
	equate( &bootavg[j][k] , BOOTS[j][k] ) ;

	memcpy( &( *mominfoavg )[ j ][ k ] ,
		&mominfo[ j ][ k ] , sizeof( struct mom_info ) ) ;
      }
    }
    // and free the x-axis
    free_mominfo( (struct mom_info**)mominfo , NSLICES ) ;

    // free BOOTS
    free_resampled_dp( (struct resampled**)BOOTS , NSLICES , 
		       INPARAMS -> NDATA ) ;
    ////
  }
  return bootavg ;
}

// average a single set of measurements
struct resampled *
momavg_single( struct mom_info **mominfoavg ,
	       int *NDATA ,
	       const struct mom_info *mominfo , 
	       const struct resampled *BOOTS )
{
  int m ;
  int Nequiv = 0 ;
  for( m = 0 ; m < *NDATA ; m++ ) {
    int k ;
    for( k = 1 ; k < *NDATA ; k++ ) { 
      if( ( m + k ) > ( *NDATA  - 1 ) ) break ;
      if( fabs( mominfo[m].p2 - mominfo[k+m].p2 ) > 1E-14 ) {
	break ;
      }
    }
    Nequiv ++ ;
    m += ( k - 1 ) ;
  }

  // allocate new momenta and bootstraps
  *mominfoavg = (struct mom_info*)malloc( Nequiv * sizeof( struct mom_info ) ) ;
    struct resampled *bootavg = malloc( Nequiv * sizeof( struct resampled ) ) ;
  int k ;
  for( k = 0 ; k < Nequiv ; k++ ) {
    bootavg[k].resampled = (double*)malloc( BOOTS[0].NSAMPLES * sizeof( double ) ) ;
  }

  // and average
  struct resampled tmp ;
  tmp.resampled = malloc( BOOTS[0].NSAMPLES * sizeof( double ) ) ;

  int idx = 0 ;
  for( m = 0 ; m < *NDATA ; m++ ) {

    equate( &tmp , BOOTS[m] ) ;
    compute_err( &tmp ) ;

    int k , nequiv = 1 ;
    memcpy( &( *mominfoavg )[ idx ] ,
	    &mominfo[ m ] , sizeof( struct mom_info ) ) ;

    for( k = 1 ; k < *NDATA ; k++ ) { 
      if( ( m + k ) > ( *NDATA  - 1 ) ) break ;
      if( fabs( mominfo[m].p2 - mominfo[k+m].p2 ) > 1E-14 ) {
	break ;
      } else {
	add( &tmp , BOOTS[k+m] ) ;
	nequiv++ ;
      }
    }
    divide_constant( &tmp , (double)nequiv ) ;
    equate( &bootavg[idx] , tmp ) ;

    idx ++ ;
    m += ( k - 1 ) ;
  }

  free( tmp.resampled ) ;

  *NDATA = Nequiv ;
  return bootavg ;
}  

// overwrites INPARAMS
struct resampled **
momavg( struct mom_info ***mominfoavg ,
	struct input_params *INPARAMS ,
	const struct mom_info **mominfo , 
	const struct resampled **BOOTS ,
	const int NSLICES ,
	const bool actually_doit )
{
  struct resampled **bootavg = (struct resampled**)malloc( NSLICES * sizeof( struct resampled* ) ) ;
  *mominfoavg = (struct mom_info**)malloc( NSLICES * sizeof( struct mom_info* ) ) ;

  if( actually_doit ) {
    // average equivalent momenta
    int j , Nequiv[ NSLICES ] ;

    for( j = 0 ; j < NSLICES ; j++ ) {
      int m ;
      Nequiv[j] = 0 ;
      for( m = 0 ; m < INPARAMS -> NDATA[j] ; m++ ) {
	int k ;
	for( k = 1 ; k < INPARAMS -> NDATA[j] ; k++ ) { 
	  if( ( m + k ) > ( INPARAMS -> NDATA[j] - 1 ) ) break ;
	  if( fabs( mominfo[j][m].p2 - mominfo[j][k+m].p2 ) > 1E-14 ) {
	    break ;
	  }
	}
	Nequiv[j] ++ ;
	m += ( k - 1 ) ;
      }
      printf( "[MOMAVG] NEQUIV :: %d \n" , Nequiv[j] ) ;
    }
    
    // allocate new momenta and bootstraps
    for( j = 0 ; j < NSLICES ; j++ ) {
      (*mominfoavg)[j] = (struct mom_info*)malloc( Nequiv[j] * sizeof( struct mom_info ) ) ;
      bootavg[j] = (struct resampled*)malloc( Nequiv[j] * sizeof( struct resampled ) ) ;
      int k ;
      for( k = 0 ; k < Nequiv[j] ; k++ ) {
	bootavg[j][k].resampled = (double*)malloc( BOOTS[j][0].NSAMPLES * sizeof( double ) ) ;
      }
    }
    
    // and repeat
    for( j = 0 ; j < NSLICES ; j++ ) {
      
      struct resampled tmp ;
      tmp.resampled = malloc( BOOTS[j][0].NSAMPLES * sizeof( double ) ) ;

      int m , idx = 0 ;
      for( m = 0 ; m < INPARAMS -> NDATA[j] ; m++ ) {

	equate( &tmp , BOOTS[j][m] ) ;
	compute_err( &tmp ) ;

	int k , nequiv = 1 ;

	memcpy( &( *mominfoavg )[ j ][ idx ] ,
		&mominfo[ j ][ m ] , sizeof( struct mom_info ) ) ;

	for( k = 1 ; k < INPARAMS -> NDATA[j] ; k++ ) { 
	  if( ( m + k ) > ( INPARAMS -> NDATA[j] - 1 ) ) break ;
	  if( fabs( mominfo[j][m].p2 - mominfo[j][k+m].p2 ) < 1E-14 ) {
	    add( &tmp , BOOTS[j][k+m] ) ;
	    nequiv++ ;
	  } else {
	    break ;
	  }
	}

	divide_constant( &tmp , (double)nequiv ) ;
	equate( &bootavg[j][idx] , tmp ) ;
	idx++ ;
	m += ( k - 1 ) ;
      }
      free( tmp.resampled ) ;
    }

    // and free the x-axis
    free_mominfo( (struct mom_info**)mominfo , NSLICES ) ;

    // free BOOTS
    free_resampled_dp( (struct resampled**)BOOTS , NSLICES , 
		       INPARAMS -> NDATA ) ;

    // set inparams to be the size of Nequiv
    for( j = 0 ; j < NSLICES ; j++ ) {
      INPARAMS -> NDATA[j] = Nequiv[j] ;
    }

  } else {
    // do almost nothing
    int j ;
    for( j = 0 ; j < NSLICES ; j++ ) {
      (*mominfoavg)[j] = ( struct mom_info* )malloc( INPARAMS -> NDATA[j] * sizeof( struct mom_info ) ) ;
      bootavg[j] = ( struct resampled* )malloc( INPARAMS -> NDATA[j] * sizeof( struct resampled ) ) ;
      int k ;
      for( k = 0 ; k < INPARAMS -> NDATA[j] ; k++ ) {
	bootavg[j][k].resampled = (double*)malloc( BOOTS[j][0].NSAMPLES * sizeof( double ) ) ;
	equate( &bootavg[j][k] , BOOTS[j][k] ) ;

	memcpy( &( *mominfoavg )[ j ][ k ] ,
		&mominfo[ j ][ k ] , sizeof( struct mom_info ) ) ;
      }
    }
    // and free the x-axis
    free_mominfo( (struct mom_info**)mominfo , NSLICES ) ;

    // free BOOTS
    free_resampled_dp( (struct resampled**)BOOTS , NSLICES , 
		       INPARAMS -> NDATA ) ;
    ////
  }

  return bootavg ;
}

struct resampled **
momentum_average( struct mom_info ***mominfoavg ,
		  struct input_params *INPARAMS ,
		  const struct mom_info **mominfo , 
		  const struct resampled **BOOTS ,
		  const int NSLICES ,
		  const momavg_type type )
{
  switch( type ) {
  case PSQ_AVERAGE :
    printf( "[MOMAVG] PSQ averaging \n" ) ;
    return momavg( mominfoavg , INPARAMS , 
		   mominfo , BOOTS , NSLICES , true ) ;
  case Z_Nm1_AVERAGE :
    printf( "[MOMAVG] Z_%d averaging \n" , 3 ) ;
    return ZN_average_slow( mominfoavg , INPARAMS , 
			    mominfo , BOOTS , 
			    NSLICES , true , 3 ) ;
    //case default :
  case NONE :
    printf( "[MOMAVG] Not momentum averaging \n" ) ;
    return momavg( mominfoavg , INPARAMS , 
		   mominfo , BOOTS , NSLICES , false ) ;
  }
  return NULL ;
}
