#ifndef POLY_H
#define POLY_H

void
poly_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
poly_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
// function evaluation
inline double
poly_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;

int
poly_f( const gsl_vector *params , 
	void *data , 
	gsl_vector *f ) ;
int
poly_df ( const gsl_vector *params , 
	  void *data , 
	  gsl_matrix * J ) ;
int
poly_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
