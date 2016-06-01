/**
   Using GSL's SVD we attempt to compute the coefficients of a polynomial

   i.e. solve the (possibly) over-constrained VanderMonde equation


   | 1 x_0 x_0^2 x_0^3 .... x_0^N | | a_0 |   | y( x_0 ) |
   | 1 x_1 x_1^2 x_1^3 .... x_1^N | | a_1 | = | y( x_1 ) |
   | ............................ | | ... | = |  .....   |
   | 1 .................... x_M^N | | a_N |   | y( x_M ) |
 
   Where M >= N and the a's are the coefficients of the polynomial

   Compiles with

   gcc test.c -lm -lgsl -lglcblas
*/

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define VERBOSE

#define FAILURE (0)
#define SUCCESS (!FAILURE)

static gsl_rng *r ;
static long unsigned int SEED = -1 ;

double
rng_double( void ) { return gsl_rng_uniform( r ) ; }

// between 0 and some number "max_idx"
int
rng_int( const int max_idx ) { 
  return gsl_rng_uniform_int( r , max_idx ) ;
}

double
rng_gaussian( const double sigma ) { 
  return gsl_ran_gaussian( r , sigma ) ; 
}

void
init_rng( long unsigned int Seed )
{
  gsl_rng_env_setup( ) ;
  
  const gsl_rng_type *type = gsl_rng_default ;
  r = gsl_rng_alloc (type) ;

  if( Seed == 0 ) {
    FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
    if( urandom == NULL ) exit( 1 ) ;
    if( fread( &Seed , sizeof( Seed ) , 1 , urandom ) != 1 ) exit(1) ;
    fclose( urandom ) ;
  }

  SEED = Seed ;
  printf( "USING GSL seed %lu \n" , SEED ) ;

  gsl_rng_set( r , Seed ) ;

  return ;
}

void
rng_reseed( void ) 
{
  if( SEED != -1 ) {
    gsl_rng_set( r , SEED ) ;
  } else {
    printf( "GSL rng not seeded properly!\n" );
    exit(1) ;
  }
  return ;
}

void
free_rng( void )
{
  gsl_rng_free( r );
  return ;
}

// little print utility
static void
printmatrix( const double **mat ,
	     const int NROWS , 
	     const int NCOLS )
{
  printf( "\n" ) ;
  int i , j ;
  for( i = 0 ; i < NROWS ; i++ ) {
    printf( "|" ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      printf( " %f " , mat[i][j] ) ;
    }
    printf( "|\n" ) ;
  }
  printf( "\n" ) ;
  return ;
}

/**
   Computes A = U S V^T for a general M * N matrix where M > N
 */
static int 
svd_inverse( double *__restrict *__restrict Ainv , 
	     const double *__restrict *__restrict A ,
	     const int NCOLS ,
	     const int NROWS )
{
  // initial case dies
  if( NROWS < NCOLS ) {
    printf( "WHAT M :: %d vs %d \n" , NROWS , NCOLS ) ;
    return FAILURE ;
  }
  if( NROWS != NCOLS ) {
    printf( "Computing the Pseudo Inverse \n" ) ;
  }

  // allocations
  gsl_matrix *m  = gsl_matrix_alloc( NROWS , NCOLS ) ;
  gsl_matrix *Q  = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_vector *S  = gsl_vector_alloc( NCOLS ) ;
  gsl_vector *WORK  = gsl_vector_alloc( NCOLS ) ;
  double Diag[ NCOLS ] , InvDiag[ NCOLS ] ;
  double diff = 0.0 ;
  register double sum ;
  double tmp ;
  int i , j , k , FLAG = SUCCESS ;

  for( i = 0 ; i < NROWS ; i++ ){
   for( j = 0 ; j < NCOLS ; j++ ) {
     gsl_matrix_set( m , i , j , A[i][j] ) ;
   }
  }

  // do the decomposition
  if( gsl_linalg_SV_decomp( m , Q , S , WORK ) ) {
    printf( "GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  // loop these set 1.0 / Diag[i] to 0.0 if Diag[i] is crazy small
  for( i = 0 ; i < NCOLS ; i++ ) {
    Diag[ i ] = gsl_vector_get( S , i ) ; 
    tmp = 1.0 / Diag[i] ;
    InvDiag[ i ] = fabs( Diag[i] ) < 1E-15 ? 0.0 : tmp ;
#ifdef VERBOSE
    printf( "SVD %d %le\n" , i , Diag[i] ) ;
    printf( "SVD %d %le\n" , i , InvDiag[i] ) ;
#endif
  } 

  // test the solution to make sure it isn't too bad
  for( i = 0 ; i < NROWS ; i++ ) {
    for( k = 0 ; k < NCOLS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	sum += gsl_matrix_get( m , i , j ) * Diag[j] * gsl_matrix_get( Q , k , j ) ;
      }
      diff += fabs( sum - A[i][k] ) ;
    }
  }

  // tell us how good the solution is and cry if it is too bad
  printf( "Decomposition accuracy :: %le \n" , diff ) ;
  diff /= (double)( NCOLS * NROWS ) ;
  if( diff > 1E-8 ) {
    printf( "SVD accuracy considered too low %e \n" , diff ) ;
    FLAG = FAILURE ;
    goto FREE ;
  } 

  // compute the product Ainv = V * ( 1.0 / Diag ) * U^T
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( k = 0 ; k < NROWS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	if( InvDiag[ j ] != 0.0 ) {
	  sum += gsl_matrix_get( Q , i , j ) * InvDiag[j] * gsl_matrix_get( m , k , j ) ;
	}
      }
      Ainv[i][k] = sum ;
    }
  }

#ifdef VERBOSE
  // test our solution, should equal the NxN identity
  double **sol = malloc( NCOLS * sizeof( double* ) ) ;
  for( i = 0 ; i < NCOLS ; i++ ) {
    sol[i] = (double*)malloc( NCOLS * sizeof( double ) ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      sum = 0.0 ;
      for( k = 0 ; k < NROWS ; k++ ) {
	sum += Ainv[i][k] * A[k][j] ;
      }
      sol[i][j] = sum ;
    }
  }
  printmatrix( (const double**)A , NROWS , NCOLS ) ;
  printmatrix( (const double**)Ainv , NCOLS , NROWS ) ;
  printmatrix( (const double**)sol , NCOLS , NCOLS ) ;

  for( i = 0 ; i < NCOLS ; i++ ) {
    free( sol[i] ) ;
  }
  free( sol ) ;
#endif

 FREE :
  gsl_matrix_free( m ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( WORK ) ;

  return FLAG ;
}

// 
static int
compute_coefficients( double *__restrict coeffs ,
		      const double *__restrict y ,
		      const double *__restrict x ,
		      const double M ,
		      const double N )
{
  // allocate the matrices
  double **A = malloc( M * sizeof( double* ) ) ;
  double **Ainv = malloc( N * sizeof( double* ) ) ;
  double x0 ;
  register double sum ;
  int i , j , FLAG = SUCCESS ;

  // matrix allocations
  for( i = 0 ; i < M ; i++ ) {
    A[i] = (double*)malloc( N * sizeof( double ) ) ;
    if( i < N ) 
      Ainv[i] = (double*)malloc( M * sizeof( double ) ) ;
  }

  // allocate the matrix with increasing powers of x
  for( i = 0 ; i < M ; i++ ) {
    x0 = 1 ;
    for( j = 0 ; j < N ; j++ ) {
      A[i][j] = x0 ;
      x0 *= x[i] ;
    }
  }

  // if the SVD screws up we set failure and free allocations
  if( svd_inverse( Ainv , (const double**)A , N , M ) == FAILURE ) {
    FLAG = FAILURE ;
    goto FREE ;
  }

  // multiply the inverse by the data to obtain the coefficients
  for( i = 0 ; i < N ; i++ ) {
    sum = 0.0 ;
    for( j = 0 ; j < M ; j++ ) {
      sum += Ainv[i][j] * y[j] ;
    }
    coeffs[i] = sum ;
  }

 FREE :
  // free this stuff
  for( i = 0 ; i < M ; i++ ) {
    free( A[i] ) ;
    if( i < N )
      free( Ainv[i] ) ;
  }
  free( A ) ;
  free( Ainv ) ;
  return FLAG ;
}

// just a helpful little routine
static void
write_polynomial( const double *__restrict coeffs ,
		  const int POLY_ORDER )
{
  int j ;
  printf( "y = " ) ;
  for( j = 0 ; j < POLY_ORDER-1 ; j++ ) {
    printf( "%f * x^{%d} + " , coeffs[j] , j ) ;
  }
  printf( "%f * x^{%d) \n" , coeffs[j] , j ) ;
  printf( "\n" ) ;
}

int main( void )
{
  // POLY_ORDER^th polynomial with NDATA datapoints
  const int POLY_ORDER = 3 ;
  const int NDATA = 10 ;

  // allocations
  double *y = malloc( NDATA * sizeof( double ) ) ;
  double *x = malloc( NDATA * sizeof( double ) ) ;
  double *rcoeffs = malloc( POLY_ORDER * sizeof( double ) ) ;
  double *coeffs = malloc( POLY_ORDER * sizeof( double ) ) ;
  register double sum ;
  register double xn ;
  int i , j ;

  // get one from the entropy pool
  init_rng( 0 ) ;

  // set the coefficients to some random numbers
  for( j = 0 ; j < POLY_ORDER ; j++ ) {
    rcoeffs[ j ] = rng_double( ) ;
  }

  // set the matrix
  for( i = 0 ; i < NDATA ; i++ ) {
    x[i] = (double)i + 1 ;
    sum = 0.0 ;
    xn = 1.0 ;
    for( j = 0 ; j < POLY_ORDER ; j++ ) {
      sum += rcoeffs[ j ] * xn ;
      xn *= x[ i ] ;
    }
    y[ i ] = sum ;
  }

  // and compute
  if( compute_coefficients( coeffs , y , x , NDATA , POLY_ORDER ) == FAILURE ) {
    goto FREE ;
  }

  // and have a look at them
  write_polynomial( rcoeffs , POLY_ORDER ) ;
  write_polynomial( coeffs , POLY_ORDER ) ;

 FREE :
  free_rng( ) ;
  free( x ) ;
  free( y ) ;
  free( rcoeffs ) ;
  free( coeffs ) ;

  return 0 ;
}
