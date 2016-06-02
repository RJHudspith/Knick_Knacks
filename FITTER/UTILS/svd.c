/**
   @file svd.c
   @brief code to compute the (psuedo)inverse of a matrix via svd using GSL

   Using GSL's SVD we attempt to compute the coefficients of a polynomial

   i.e. solve the (possibly) over-constrained VanderMonde equation


   | 1 x_0 x_0^2 x_0^3 .... x_0^N | | a_0 |   | y( x_0 ) |
   | 1 x_1 x_1^2 x_1^3 .... x_1^N | | a_1 | = | y( x_1 ) |
   | ............................ | | ... | = |  .....   |
   | 1 .................... x_M^N | | a_N |   | y( x_M ) |
 
   Where M >= N and the a's are the coefficients of the polynomial
 */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SVD_COL_BALANCE

//#define VERBOSE

#define FAILURE 0
#define SUCCESS !FAILURE

#ifdef VERBOSE
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
#endif

/**
   Computes A = U S V^T for a general M * N matrix where M > N
 */
int 
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
    //printf( "Computing the Pseudo Inverse \n" ) ;
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
#ifdef SVD_COL_BALANCE
  gsl_vector *D = gsl_vector_alloc( NCOLS ) ;
  gsl_linalg_balance_columns( m , D ) ;
#endif

  if( gsl_linalg_SV_decomp( m , Q , S , WORK ) ) {
    printf( "GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  // loop these set 1.0 / Diag[i] to 0.0 if Diag[i] is crazy small
  for( i = 0 ; i < NCOLS ; i++ ) {
    Diag[ i ] = gsl_vector_get( S , i ) ; 
    tmp = 1.0 / Diag[i] ;
    InvDiag[ i ] = fabs( Diag[i] ) < 1.E-32 ? 0.0 : tmp ;
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
      #ifdef SVD_COL_BALANCE
      diff += fabs( sum - A[i][k] / gsl_vector_get( D , k ) ) ;
      #else
      diff += fabs( sum - A[i][k] ) ;
      #endif
    }
  }

  // tell us how good the solution is and cry if it is too bad
  // printf( "Decomposition accuracy :: %le \n" , diff ) ;
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
      #ifdef SVD_COL_BALANCE
      Ainv[i][k] = sum / gsl_vector_get( D , i ) ;
      #else
      Ainv[i][k] = sum ;
      #endif
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

#ifdef SVD_COL_BALANCE
  gsl_vector_free( D ) ;
#endif

  return FLAG ;
}

// 
int
compute_coefficients( double *__restrict coeffs ,
		      double *__restrict chisq ,
		      const double *__restrict y ,
		      const double *__restrict sigma ,
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

  // compute residual r += ( y[ i ] - A * c )
  *chisq = 0.0 ;
  register double ri ;
  for( j = 0 ; j < M ; j++ ) {
    register double res = 0.0 ;
    for( i = 0 ; i < N ; i++ ) {
      res += A[j][i] * coeffs[i] ;
    }
    ri = ( y[ j ] - res ) / ( sigma[ j ] ) ;
    *chisq += ri * ri ;
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

// pade coefficient function
int
pades( double *pade_coeffs ,
       const double *poly_coeffs ,
       const int n ,
       const int m )
{
  // is guaranteed
  pade_coeffs[ 0 ] = poly_coeffs[ 0 ] ;

  // I do the simple cases here 
  switch( n ) {
  case 1 :
    switch( m ) {
    case 1 :     // ( 1 , 1 ) pade
      pade_coeffs[ 1 ] = poly_coeffs[ 1 ] ;
      pade_coeffs[ 2 ] = -poly_coeffs[ 2 ] / poly_coeffs[ 1 ] ;
      return SUCCESS ;
    case 2 :     // ( 1 , 2 ) pade
      pade_coeffs[ 1 ] = poly_coeffs[ 1 ] ;
      pade_coeffs[ 2 ] = -poly_coeffs[ 2 ] / poly_coeffs[ 1 ] ;
      pade_coeffs[ 3 ] = ( poly_coeffs[2] * poly_coeffs[2] - poly_coeffs[1]*poly_coeffs[3]) / 
	( poly_coeffs[1] * poly_coeffs[1] ) ;
      return SUCCESS ;
    }
    break ;
  case 2 :
    switch( m ) {
    case 1 :     // ( 2 , 1 ) pade
      pade_coeffs[ 1 ] = poly_coeffs[ 1 ] ;
      pade_coeffs[ 2 ] = poly_coeffs[2] - ( poly_coeffs[ 3 ] * poly_coeffs[ 1 ] ) / poly_coeffs[ 2 ] ;
      pade_coeffs[ 3 ] = -poly_coeffs[ 3 ] / poly_coeffs[ 2 ] ;
      return SUCCESS ;
    }
    break ;
  }

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

  // numerator is first in our scheme
  for( i = 0 ; i < n-1 ; i++ ) {
    pade_coeffs[ i + 1 ] = right_side[ i ] ;
  }
  pade_coeffs[ i + 1 ] = right_side[ i ] ;
  // denominator is last
  for( i = 0 ; i < m-1 ; i++ ) {
    pade_coeffs[ i + n + 1 ] = left_side[ i + 1 ] ;
  } 
  pade_coeffs[ i + n + 1 ] = left_side[ i + 1 ] ;

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

// just a helpful little routine
void
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

#ifdef SVD_COL_BALANCE
  #undef SVD_COL_BALANCE
#endif
