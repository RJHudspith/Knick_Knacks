/**
   @file pade_from_poly.c
   @brief computes pade coefficients from a polynomial

   For a given polynomial

   C(q^2) = C_0 + C_1 q^2 + C_2 q^4 + ...

   This code approximates the pade coefficients for a given

   ( n , m ) pade

   C(q^2) = ( P_0 + P_1 q^2 + P_2 q^4 + .. P_n q^{2n} ) 
                  / ( 1 + P_{n+1} q^2 + P_{n+2} q^4 + ... )


   This is solved by noticing that, C_0 = P_0 

   | C_1       0                 |  | 1       |   | P_0 |
   | C_2      C_1       0        |  | P_{n+1} |   | P_1 |
   | C_3      C_2                |  | .....   | = | ... |
   | ....     ...      ...    0  |  | P_{n+m} |   | P_n |
   | C_{n+m}  C_{n+m}  ...   C_1 |  |    0    |   |  0  |

   This can be used to compute pade representations for Taylor 
   series and whatnot
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

#define FAILURE 0
#define SUCCESS !FAILURE

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
#ifdef VERBOSE
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
#endif

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
    printf( "SVD %d %f\n" , i , Diag[i] ) ;
    printf( "SVD %d %f\n" , i , InvDiag[i] ) ;
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

// pade coefficient function
static int
pades( double *pade_coeffs ,
       const double *poly_coeffs ,
       const int n ,
       const int m )
{
  // the matrix order equation we must solve
  const int order = n + m + 1 ;

  // initialise solve matrix and its inverse
  double **sub_solve = malloc( m * sizeof( double* ) ) ;
  double **Ainv = malloc( m * sizeof( double* ) ) ;

  // vectors we use
  double *solves = malloc( m * sizeof( double ) ) ;
  double *left_side = malloc( ( order - 1 ) * sizeof( double ) ) ;
  double *right_side = malloc( ( order - 1 ) * sizeof( double ) ) ;

  // flag for whether it worked
  int FLAG = FAILURE ;

  // set up the submatrix we need to solve for the
  // left hand side
  int i , j ;
  for( i = 0 ; i < m ; i++ ) {
    sub_solve[ i ] = malloc( m * sizeof( double) ) ;
    Ainv[ i ] = malloc( m * sizeof( double) ) ;
    for( j = 0 ; j < m ; j++ ) {
      if( i - j + n > 0 ) {
	sub_solve[i][j] = poly_coeffs[ i - j + n ] ;
      } else {
	sub_solve[i][j] = 0.0 ;
      }
    }
    solves[i] = -poly_coeffs[ i + n + 1 ] ;
    #ifdef VERBOSE
    printf( "sols [ %d ] = %d :: %f \n" , i , m + i + 1 , solves[i] ) ;
    #endif
  }

  // printhimout
#ifdef VERBOSE
  printmatrix( (const double**)sub_solve , m , m ) ;
#endif

  // perform an inverse
  if( svd_inverse( Ainv , (const double**)sub_solve , m , m ) == FAILURE ) {
    FLAG = FAILURE ;
    goto FREE ;
  }

  // multiply ainverse by sols to get left hand side coefficients
  left_side[ 0 ] = 1.0 ;
  for( i = 1 ; i < order - 1 ; i++ ) {
    left_side[ i ] = 0.0 ;
    if( i <= m ) {
      for( j = 0 ; j < m ; j++ ) {
	left_side[ i ] += Ainv[i-1][j] * solves[j] ;
      }
    }
  }

  // multiply coefficient matrix by left hand side
  for( i = 0 ; i < order-1 ; i++ ) {
    register double sum = poly_coeffs[ i + 1 ] ;
    for( j = 1 ; j <= i ; j++ ) {
      sum += poly_coeffs[ i - j + 1 ] * left_side[ j ] ;
    }
    right_side[ i ] = sum ;
    #ifdef VERBOSE
    printf( "\n" ) ;
    printf( "right[ %d ] = %f \n" , i , right_side[ i ] ) ;
    #endif
  }

  // easiest one now
  pade_coeffs[ 0 ] = poly_coeffs[ 0 ] ;
  printf( "PADE( %d , %d ) :: %1.2f + ( " , n , m , pade_coeffs[ 0 ] ) ;

  // numerator is first in our scheme
  for( i = 0 ; i < n-1 ; i++ ) {
    pade_coeffs[ i + 1 ] = right_side[ i ] ;
    printf( " (%1.2f) q^{%d} +" , pade_coeffs[ i + 1 ] , 2*(i+1) ) ;
  }
  pade_coeffs[ i + 1 ] = right_side[ i ] ;
  printf( " (%1.2f) q^{%d} ) / ( 1 +" , pade_coeffs[ i + 1 ] , 2*(i+1) ) ;

  // denominator is last
  for( i = 0 ; i < m-1 ; i++ ) {
    pade_coeffs[ i + n + 1 ] = left_side[ i + 1 ] ;
    printf( " (%1.2f) q^{%d} +" , pade_coeffs[ i + n + 1 ] , 2*(i+1) ) ;
  } 
  pade_coeffs[ i + n + 1 ] = left_side[ i + 1 ] ;
  printf( " (%1.2f) q^{%d} ) \n" , pade_coeffs[ i + n + 1 ] , 2*(i+1) ) ;

  // free workspaces
 FREE :
  for( i = 0 ; i < m ; i++ ) {
    free( sub_solve[i] ) ;
    free( Ainv[i] ) ;
  }
  free( sub_solve ) ;
  free( Ainv ) ;
  free( solves ) ;
  free( left_side ) ;
  free( right_side ) ;

  return FLAG ;
}

//
int main( void )
{
  const int n = 2 ;
  const int m = 2 ;
  const int order = n + m + 1 ;

  double poly_coeffs[ order ] ; // inputs
  double pade_coeffs[ order ] ; // outputs

  // first nine coefficients
  poly_coeffs[0] = 1 ;
  poly_coeffs[1] = 1 ;
  poly_coeffs[2] = 1./2. ;
  poly_coeffs[3] = 1./6. ;
  poly_coeffs[4] = 1./24. ;

  pades( pade_coeffs , poly_coeffs , n , m ) ;

  int i ;
  for( i = 0 ; i < order ; i++ ) {
    printf( "%1.15f \n" , pade_coeffs[i] ) ;
  }

  return 0 ;
}
