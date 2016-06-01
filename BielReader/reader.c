/**
   @file reader.c
   @brief little code to turn quenched SU3 (Bielefeld) to NERSC configs
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <complex.h>
#include <math.h>

// success or failure
#define FAILURE (-1)
#define SUCCESS (!FAILURE)

// 4 dimensions
#define ND (4)

// 3 colors
#define NC (3)
#define NCNC (NC*NC)

// precision tolerance
#define PREC_TOL (1.E-15)

// argv enum
enum { DIMENSIONS = 1 , INFILE = 2 , OUTFILE = 3 } ;

// NERSC output enum
typedef enum { OUTPUT_GAUGE , NCxNC } outtype ;

// site is a lattice link
struct site {
  double complex U[ ND ][ NCNC ] ;
  size_t neighbor[ ND ] ;
  size_t back[ ND ] ;
} ;

// lattice stuff
struct latt_info {
  size_t dims[ ND ] ;
  size_t LCU ;
  size_t LVOLUME ;
} ;

// global lattice information
struct latt_info Latt ;

// computes the determinant
static double complex
det( const double complex U[ NCNC ] )
{
  return U[0] * ( U[4] * U[8] - U[5] * U[7] ) -
    U[1] * ( U[3] * U[8] - U[5] * U[6] ) +
    U[2] * ( U[3] * U[7] - U[4] * U[6] ) ;
} 

// check for unitarity of links
static bool
is_unitary( const double complex U[ NCNC ] )
{
  size_t j , i ; 
  double vv = 0. ; 
  for( j = 0 ; j < NC ; j++ ) {
    for( i = 0 ; i < NC ; i++ ) {
      vv += creal( U[ NC * i ] * conj( U[ NC * i + j ] ) ) ; 
    }
    if( fabs( vv - 1.0 )/NC > PREC_TOL ) {
      printf( "[column %zu] not orthogonal!!! ->  MUST REUNITARIZE vv %1.8f\n" , j , vv ) ; 
      //return false ;
    }
  }
  vv = 0 ; 
  //check if rows are orthogonal
  for( j = 0 ; j < NC ; j++ ) {
    for( i = 0 ; i < NC ; i++ ) {
      vv += creal( U[ i ] * conj( U[ i + NC * j ] ) ) ; 
    }
    if( fabs( vv - 1.0 )/NC > PREC_TOL ) {
      printf( "[row %zu] not orthogonal!!! -> MUST REUNITARIZE vv %1.8f\n" , j , vv ) ; 
      //return false ;
    }
  }
  return true ;
}

// read a checkerboarded timeslice
static int
read_timeslice( struct site *lat ,
		const char *filename ,
		const size_t timeslice )
{
  char tempstr[ 512 ] , tstr[ 5 ] ;
  // conversion because of Fortran
  if( timeslice < 10 ) {
    sprintf( tstr , "0%zu" , timeslice ) ;
  } else {
    sprintf( tstr , "%zu" , timeslice ) ;
  }
  sprintf( tempstr , filename , tstr ) ;
  printf( "Reading :: %s \n" , tempstr ) ;

  // read the binary data
  FILE *infile = fopen( tempstr , "rb" ) ;
  if( infile == NULL ) {
    printf( "404 File not found %s\n" , tempstr ) ;
    return FAILURE ;
  }

  size_t i , mu ;
  for( i = 0 ; i < Latt.LCU ; i++ ) {

    // read in all the links on a site
    double complex uout[ ND*NC*(NC-1) ] ;
    if( fread( uout , sizeof( double complex ) , 
	       ND*NC*(NC-1) , infile ) != ND*NC*(NC-1) ) {
      printf( "File read failure\n" ) ;
      return FAILURE ;
    }

    // loop mu
    for( mu = 0 ; mu < ND ; mu++ ) {

      // poke into this link ...
      lat[i].U[mu][0] = uout[ 0 + 6*mu ] ;
      lat[i].U[mu][1] = uout[ 2 + 6*mu ] ;
      lat[i].U[mu][2] = uout[ 4 + 6*mu ] ;

      lat[i].U[mu][3] = uout[ 1 + 6*mu ] ;
      lat[i].U[mu][4] = uout[ 3 + 6*mu ] ;
      lat[i].U[mu][5] = uout[ 5 + 6*mu ] ;

      // complete by cross product
      lat[i].U[mu][6] = conj( lat[i].U[mu][1] * lat[i].U[mu][5] - lat[i].U[mu][2] * lat[i].U[mu][4] ) ;
      lat[i].U[mu][7] = -conj( lat[i].U[mu][0] * lat[i].U[mu][5] - lat[i].U[mu][2] * lat[i].U[mu][3] ) ;
      lat[i].U[mu][8] = conj( lat[i].U[mu][0] * lat[i].U[mu][4] - lat[i].U[mu][1] * lat[i].U[mu][3] ) ;

      dt = det( lat[i].U[mu] ) ;

      // make sure the resulting matrix is special unitary
      if( cabs( det( lat[i].U[mu] ) - 1. ) > PREC_TOL || 
	  is_unitary( lat[i].U[mu] ) == false ) {
	printf( "Matrix is not Special Unitary!\n" ) ;
	return FAILURE ;
      }
      //
    }
  }
  fclose( infile ) ;
  return SUCCESS ; 
}

// initialise the lattice geometry
static int
init_geometry( char *dimstr )
{
  char *tok1 = strtok( dimstr , "," ) ;
  Latt.dims[0] = (size_t)atoi( tok1 ) ;
  if( Latt.dims[ 0 ] < 1 ) {
    printf( "Non-sensical lattice dimension L_%zu :: %zu \n" , 
	    (size_t)0 , Latt.dims[ 0 ] ) ;
    return FAILURE ;
  }
  size_t mu = 1 ;
  while( ( tok1 = strtok( NULL , "," ) ) != NULL ) {
    Latt.dims[ mu ] = (size_t)atoi( tok1 ) ;
    if( Latt.dims[ mu ] < 1 ) {
      printf( "Non-sensical lattice dimension L_%zu :: %zu \n" , 
	      mu , Latt.dims[ mu ] ) ;
      return FAILURE ;
    }
    if( mu == (ND-1) ) break ;
    mu++ ;
  }
  // set LCU and LVOLUME
  printf( "\nDimensions (" ) ;
  Latt.LCU = 1 ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    printf( " %zu x" , Latt.dims[mu] ) ;
    Latt.LCU *= Latt.dims[mu] ;
  }
  printf( " %zu )\n\n" , Latt.dims[ ND-1 ] ) ; 
  Latt.LVOLUME = Latt.LCU * Latt.dims[mu] ;

  printf( "LCU :: %zu , LVOLUME :: %zu \n" , Latt.LCU , Latt.LVOLUME ) ;
  return SUCCESS ;
}

// computes the lexicographical site index from the position vector in x
static size_t
gen_site( const size_t x[ ND ] )
{
  size_t res = x[ ND - 1 ] ;
  size_t mu ;
  for( mu = ND - 1 ; mu > 0 ; mu-- ) {
    res = Latt.dims[ mu - 1 ] * res + x[ mu - 1 ] ;
  }
  return res ;
}

// generic version of the below
static void 
get_mom_2piBZ( size_t x[ ND ] , 
	       const size_t i , 
	       const size_t DIMS )
{
  size_t mu , subvol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu != DIMS ) {
      x[ mu ] = ( ( i - i % subvol ) / subvol ) % Latt.dims[ mu ] ;
      subvol *= Latt.dims[ mu ] ;
    } else {// set it to 0?
      x[ mu ] = 0 ;
    }
  }
  return ;
}

// This is the new bona-fide generic shifting code
static size_t
gen_shift( size_t x[ ND ] ,
	   const int dir )
{
  if( dir >= 0 ) {
    x[ dir ] = ( x[ dir ] + 1 ) % Latt.dims[ dir ] ;
  } else {
    register const size_t numu = -dir - 1 ;
    if( x[ numu ] == 0 ) {
      x[ numu ] = x[ numu ] - 1 + Latt.dims[ numu ] ;
    } else {
      x[ numu ] = ( x[ numu ] - 1 ) % Latt.dims[ numu ] ;
    }
  }
  return gen_site( x ) ; 
}

// initialise the navigation
static void
init_navig( struct site *__restrict lat )
{
  size_t i ; 
  #pragma omp parallel for private(i)
  for(  i = 0 ; i < Latt.LVOLUME ; i++  )  {
    size_t x[ ND ] ; 
    get_mom_2piBZ( x , i , ND ) ;
    size_t mu ;
    for(  mu = 0 ; mu < ND ; mu++  )	{
      lat[i].neighbor[mu] = gen_shift( x , mu ) ; 
      lat[i].back[mu] = gen_shift( x , -mu - 1 ) ;  
    }
  }
  return ;
}

// set the cartesian coordinates
static void
cartesian_coordinates( size_t cx[ ND ] ,
		       size_t *cb ,
		       const size_t idx )
{
  size_t mu , si = idx , sum = 0 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    cx[ mu ] = si % Latt.dims[ mu ] , si /= Latt.dims[mu] ;
    sum += cx[ mu ] ;
  }
  *cb = sum & 1 ;
  return ;
}		       

// convert checkerboard to lexicographical
static size_t
cb_to_lexi( const size_t i )
{
  size_t cx[ ND ] , cb ;
  cartesian_coordinates( cx , &cb , i ) ;
  const size_t ix = ( cx[0] + Latt.dims[0] * ( cx[1] + Latt.dims[1] * cx[2] ) ) / 2 ;

  return cb + 2*ix ;
}

// compute the plaquette
static double
complete_plaquette( const double complex *__restrict a , 
		    const double complex *__restrict b , 
		    const double complex *__restrict c , 
		    const double complex *__restrict d )
{
  register double complex tra1 = (a[0]*b[0]+a[1]*b[3]+a[2]*b[6]) ;
  register double complex tra2 = (a[0]*b[1]+a[1]*b[4]+a[2]*b[7]) ;
  register double complex tra3 = (a[0]*b[2]+a[1]*b[5]+a[2]*b[8]) ;
  register double tra = creal( ( tra1 * conj(c[0]) + tra2 * conj( c[1] ) + tra3 * conj( c[2] ) ) * conj( d[0] ) ) ;
  tra += creal( ( tra1 * conj(c[3]) + tra2 * conj( c[4] ) + tra3 * conj( c[5] ) ) * conj( d[1] ) ) ;
  tra += creal( ( tra1 * conj(c[6]) + tra2 * conj( c[7] ) + tra3 * conj( c[8] ) ) * conj( d[2] ) ) ;
  tra1 = (a[3]*b[0]+a[4]*b[3]+a[5]*b[6]) ;
  tra2 = (a[3]*b[1]+a[4]*b[4]+a[5]*b[7]) ;
  tra3 = (a[3]*b[2]+a[4]*b[5]+a[5]*b[8]) ;
  register double trb = creal( ( tra1 * conj(c[0]) + tra2 * conj( c[1] ) + tra3 * conj( c[2] ) ) * conj( d[3] ) ) ;
  trb += creal( ( tra1 * conj(c[3]) + tra2 * conj( c[4] ) + tra3 * conj( c[5] ) ) * conj( d[4] ) ) ;
  trb += creal( ( tra1 * conj(c[6]) + tra2 * conj( c[7] ) + tra3 * conj( c[8] ) ) * conj( d[5] ) ) ;
  tra1 = (a[6]*b[0]+a[7]*b[3]+a[8]*b[6]) ;
  tra2 = (a[6]*b[1]+a[7]*b[4]+a[8]*b[7]) ;
  tra3 = (a[6]*b[2]+a[7]*b[5]+a[8]*b[8]) ;
  register double trc = creal( ( tra1 * conj(c[0]) + tra2 * conj( c[1] ) + tra3 * conj( c[2] ) ) * conj( d[6] ) ) ;
  trc += creal( ( tra1 * conj(c[3]) + tra2 * conj( c[4] ) + tra3 * conj( c[5] ) ) * conj( d[7] ) ) ;
  trc += creal( ( tra1 * conj(c[6]) + tra2 * conj( c[7] ) + tra3 * conj( c[8] ) ) * conj( d[8] ) ) ;
  return tra + trb + trc ;
}

// all of the plaquettes
double
all_plaquettes( const struct site *__restrict lat ,
		double *__restrict sp_plaq ,
		double *__restrict t_plaq )
{
  double spplaq = 0. , tplaq = 0.0 ;
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:spplaq) reduction(+:tplaq) 
  for( i = 0 ; i < Latt.LVOLUME ; i++ ) {
    double p = 0. , face ;
    size_t mu , nu , s , t ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      t = lat[i].neighbor[mu] ; 
      for( nu = 0 ; nu < mu ; nu++ ) {
        s = lat[i].neighbor[nu] ;
	face = complete_plaquette( lat[ i ].U[mu] , lat[ t ].U[nu] , 
				   lat[ s ].U[mu] , lat[ i ].U[nu] ) ; 
	p = p + (double)face ;
      }
    }
    spplaq = spplaq + (double)p ;
    // reinitialise p for the temporal plaquette ...
    p = 0.0 ;
    t = lat[i].neighbor[mu] ; 
    for( nu = 0 ; nu < mu ; nu++ ) {
      s = lat[i].neighbor[nu] ;
      face = complete_plaquette( lat[ i ].U[mu] , lat[ t ].U[nu] , 
				 lat[ s ].U[mu] , lat[ i ].U[nu] ) ; 
      p = p + (double)face ;
    }
    tplaq = tplaq + (double)p ;
  }
  *sp_plaq = 2.0 * spplaq / (double)( ( ND - 1 ) * ( ND - 2 ) * NC * Latt.LVOLUME ) ;
  *t_plaq  = 2.0 * tplaq / (double)( ( ND - 1 ) * ( ND - 2 ) * NC * Latt.LVOLUME ) ;
  return 0.5 * ( *sp_plaq + *t_plaq ) ; 
}

// convert lat_in ordering to NERSC
static int
convert_to_NERSC( struct site *lat_out ,
		  struct site *lat_in , 
		  const size_t timeslice )
{
  // is the 4D index
  const size_t tidx = timeslice * Latt.LCU ;
  size_t i ;
  for( i = 0 ; i < Latt.LCU ; i++ ) {
    const size_t idx = cb_to_lexi( i + tidx ) ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      memcpy( lat_out[ i + tidx ].U[mu] , lat_in[ idx ].U[mu] , NCNC * sizeof( double complex ) ) ;
    }
  }
  return SUCCESS ;
}

// compute the trace
static double
link_trace( const struct site *__restrict lat )
{
  double tr = 0.0 ;
  size_t i ;
#pragma omp parallel for private(i) reduction(+:tr)
  for( i = 0 ; i < Latt.LVOLUME ; i++ ) {
    register double loc_tr = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      loc_tr += 
	creal( lat[i].U[mu][0] ) + 
	creal( lat[i].U[mu][4] ) +
	creal( lat[i].U[mu][8] ) ;
    }
    tr = tr + (double)loc_tr ;
  }
  return tr / (double)( ND * NC * Latt.LVOLUME ) ;
}

// union for type-punning
typedef union 
{
  double val ;
  uint32_t chk[2] ;
} U32 ;

// compute the NERSC checksum
static uint32_t 
NERSC_checksum( struct site *__restrict lat ,
		const size_t LOOP )
{
  size_t i ;
  uint32_t csum = 0 ;
#pragma omp parallel for private(i) reduction(+:csum)
  for( i = 0 ; i < Latt.LVOLUME ; i++ ) {
    register uint32_t loc_csum = 0 ;
    double *p ;
    U32 in ;
    size_t j , mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      p = (double*)lat[i].U[mu] ;
      for( j = 0 ; j < LOOP ; j++ ) {
	in.val = *p ;
	loc_csum += ( in.chk[0] + in.chk[1] ) ;
	p++ ;
      }
    }
    csum = csum + (uint32_t)loc_csum ;
  }
  return csum ;
}

// writes the usual nersc header ...
static void
write_header_NERSC( FILE *__restrict out ,
		    const double tr ,
		    const double plaq ,
		    const uint32_t chksum ,
		    const char *details ,
		    const outtype type )
{
 // fill in the normal header crap //
 fprintf( out , "BEGIN_HEADER\n" ) ; 
 fprintf( out , "HDR_VERSION = 1.0\n" ) ; 
 switch( type ) {
 case OUTPUT_GAUGE : 
   fprintf( out , "DATATYPE = %dD_SU%d_GAUGE\n" , ND , NC ) ; 
   break ;
 case NCxNC :
 default :
   fprintf( out , "DATATYPE = %dD_SU%d_GAUGE_%dx%d\n" , ND , NC , NC , NC ) ;  
   break ;
 }
 fprintf( out , "STORAGE_FORMAT = 1.0\n" ) ; 
 size_t mu ;
 for( mu = 0 ; mu < ND ;  mu++ ) {
   fprintf( out , "DIMENSION_%zu = %zu\n" , mu + 1 , Latt.dims[mu] ) ; 
 }
 fprintf( out , "CHECKSUM = %x\n" , chksum ) ; 
 fprintf( out , "LINK_TRACE = %1.15f\n" , tr ) ; 
 fprintf( out , "PLAQUETTE = %1.15f\n" , plaq ) ; 
 // only use PERIODIC
 for( mu = 0 ; mu < ND ; mu++ ) {
   fprintf( out , "BOUNDARY_%zu = PERIODIC\n" , mu + 1 ) ; 
 }
 fprintf( out , "ENSEMBLE_ID = canlat15\n" ) ; 
 fprintf( out , "SEQUENCE_NUMBER = %d\n" , 0 ) ; 
 fprintf( out , "ENSEMBLE_LABEL = %s\n" , details ) ; 
 fprintf( out , "CREATOR_MACHINE = %s\n" , getenv( "USER" ) ) ; 
 //this bit's a pain -> going to have to fight the machine endianness
 fprintf( out , "FLOATING_POINT = IEEE64LITTLE\n" ) ; 
 fprintf( out , "END_HEADER\n" ) ; 
 return ;
}

// write out the binary data
static void
write_binary_data( FILE *outfile , 
		   const struct site *__restrict lat ,
		   const size_t LOOP ) 
{
  size_t i , mu ;
  for( i = 0 ; i < Latt.LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      fwrite( (double*)lat[i].U[mu] , sizeof( double ) , LOOP , outfile ) ;
    }
  }
  return ;
}

// expects ./READ Lx,Ly,Lz,Lt {infile} {outfile}
int 
main( const int argc , 
      char *argv[] )
{
  // check for the correct number of arguments
  if( argc != 4 ) {
    return printf( "usage ./VPF lx,ly,lz,lt {infile_T$d.g} {outfile}\n" ) ;
  }

  // set up lattice geometry
  if( init_geometry( argv[ DIMENSIONS ] ) == FAILURE ) {
    return FAILURE ;
  }

  // lat_out is the output NERSC format lattice
  struct site *lat_out = NULL ;
  lat_out = malloc( Latt.LVOLUME * sizeof( struct site ) ) ;
  init_navig( lat_out ) ;

  // read the slices in
  bool error = false ;
  size_t t ;
#pragma omp parallel for private(t)
  for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
    // allocate the in lat
    struct site *lat_in = NULL ;
    lat_in = malloc( Latt.LCU * sizeof( struct site ) ) ;
    // read a timeslice into lat_in
    if( read_timeslice( lat_in , argv[ INFILE ] , 
			t+1 ) == FAILURE ) {
      printf( "Timeslice read failure %zu\n" , t ) ;
      error = true ;
    }
    // poke into NERSC order
    convert_to_NERSC( lat_out , lat_in , t ) ;
    free( lat_in ) ;
  }

  // if we errd we leave
  if( error == true ) goto memfree ;

  // print plaquettes
  double sp_plaq , t_plaq ;
  const double plaq = all_plaquettes( lat_out , &sp_plaq , 
				      &t_plaq ) ;
  printf( "\nPlaquettes ( xyz , t ) :: ( %f , %f )\n" ,
	  sp_plaq , t_plaq ) ;

  // compute the average link trace
  const double tr = link_trace( lat_out ) ;
  printf( "\nLink trace :: %e\n" , tr ) ;

  // set the loop value
  outtype type = OUTPUT_GAUGE ;
  size_t LOOP = NC * ( NC - 1 ) * 2 ;
  switch( type ) {
  case OUTPUT_GAUGE :
    LOOP = NC * ( NC - 1 ) * 2 ;
    break ;
  case NCxNC :
  default :
    LOOP = NCNC * 2 ;
    break ;
  }

  // compute the checksum
  const uint32_t chksum = NERSC_checksum( lat_out , LOOP ) ;
  printf( "NERSC Checksum :: %x\n" , chksum ) ;

  // write out a NERSC file
  FILE *outfile = fopen( argv[ OUTFILE ] , "wb" ) ;
  write_header_NERSC( outfile , tr , plaq , chksum , 
		      "wflow_runs" , OUTPUT_GAUGE ) ;

  printf( "\nWriting a NERSC gauge configuration to %s \n\n" , argv[ OUTFILE ] ) ;
  
  // write the binary data
  write_binary_data( outfile , lat_out , LOOP ) ;

  fclose( outfile ) ;
  
 memfree :

  // free the slices
  free( lat_out ) ;

  return SUCCESS ;
}
