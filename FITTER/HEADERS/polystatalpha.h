#ifndef POLYSTATALPHA_H
#define POLYSTATALPHA_H

void
polystatalpha_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
polystatalpha_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
// function evaluation
inline double
polystatalpha_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;

int
polystatalpha_f( const gsl_vector *params , 
	void *data , 
	gsl_vector *f ) ;
int
polystatalpha_df ( const gsl_vector *params , 
	  void *data , 
	  gsl_matrix * J ) ;
int
polystatalpha_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
