/**
   Simple fit to the gluon using a gribov-stingl sort of dealy

   Y = X.X * ( x[2] * ( 1.0 + x[3] * X.X ) ) / ( ( X.X + x[0] )^2 + x[1] * x[1] )

   derivatives
   dY/da = -X.X * ( x[2] * ( 1.0 + x[3] * X.X ) ) * 4.0 * x[0] * ( x[0]^2 + X.X ) / ( ( X.X + x[0] )^2 + x[1] )^2
   dY/db = -X.X * ( x[2] * ( 1.0 + x[3] * X.X ) ) * 2.0 * x[1] / ( ( X.X + x[0] )^2 + x[1] * x[1] )^2
   dY/dc = X.X * ( 1.0 + x[3] * X.X ) / ( ( X.X + x[0] )^2 + x[1] * x[1] )
   dY/dd = X.X * ( x[2] * X.X ) / ( ( X.X + x[0] )^2 + x[1] * x[1] )
 */
#include "fitfunc.h"

// evaluate the function at the point "X" using the fitparams
double
glue_eval( double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  return X.X * ( x[2] * ( 1.0 + x[3] * X.X ) ) / 
    ( pow( X.X + x[0] , 2 ) + x[1] ) ;
}

// compute the difference of the function from the data
int
glue_f( const gsl_vector *x , 
	void *data , 
	gsl_vector *f )
{
  //const size_t n  = ( (struct data *)data)->n ;
  const double *y = ( (struct data *)data)->y ;
  const double *X = ( (struct data *)data)->X ;
  const double *sigma = ( (struct data *) data)->sigma ;
  const int HI = ( (struct data *) data)->FIT_HI ;
  const int LO = ( (struct data *) data)->FIT_LO ;
  
  const double A = gsl_vector_get( x , 0 ) ;
  const double B = gsl_vector_get( x , 1 ) ;
  const double C = gsl_vector_get( x , 2 ) ;
  const double D = gsl_vector_get( x , 3 ) ;

  size_t i , count = 0 ;
  for ( i = LO ; i <= HI ; i++ ) {
    const double YI = X[i] * ( C * ( 1.0 + D * X[i] ) ) / 
      ( pow( X[i] + A*A , 2 ) + B*B ) ;
    gsl_vector_set ( f , count , 
		     ( YI - y[i] ) / sigma[i] ) ;
    count++ ;
  }

  return GSL_SUCCESS;
}

int
glue_df( const gsl_vector *x, 
	 void *data, 
	 gsl_matrix *J )
{
  const double *sigma = ((struct data *) data)->sigma ;
  const double *X = ((struct data *) data)->X ;

  const int HI = ( (struct data *) data)->FIT_HI ;
  const int LO = ( (struct data *) data)->FIT_LO ;

  const double A = gsl_vector_get( x , 0 ) ;
  const double B = gsl_vector_get( x , 1 ) ;
  const double C = gsl_vector_get( x , 2 ) ;
  const double D = gsl_vector_get( x , 3 ) ;

  size_t i , count = 0 ; 
  for (i = LO ; i <= HI ; i++) {
    const double s = sigma[i] ;
    const double denom = 1.0 / ( pow( A*A + X[i],2) + B*B ) ;
    const double numer = X[i] * ( C * ( 1.0 + D * X[i] ) ) ;
    gsl_matrix_set ( J , count , 0 , -4.0 * A * ( A*A + X[i] ) * numer * denom * denom / s ) ;
    gsl_matrix_set ( J , count , 1 , -2.0 * B *  numer * denom * denom / s ) ;
    gsl_matrix_set ( J , count , 2 , X[i] * ( 1.0 + D * X[i] )  * denom / s ) ;
    gsl_matrix_set ( J , count , 3 , X[i] * ( C * X[i] ) * denom / s ) ;
    count ++ ;
  }
  return GSL_SUCCESS;
}

int
glue_fdf (const gsl_vector *x , 
	  void *data ,
          gsl_vector *f , 
	  gsl_matrix *J )
{
  glue_f(x, data, f);
  glue_df(x, data, J);

  return GSL_SUCCESS;
}
