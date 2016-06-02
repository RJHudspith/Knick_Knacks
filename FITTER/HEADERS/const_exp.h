#ifndef CONST_EXPFIT_H
#define CONST_EXPFIT_H

void
const_exp_description( const char *message ,
		       const double *__restrict params ,
		       const struct x_descriptor X ,
		       const int NPARAMS ) ;
void
const_exp_guess( double *params ,
		 void *data ,
		 const int NPARAMS ) ;
// function evaluation
inline double
const_exp_eval( const double *x ,
		const struct x_descriptor X ,
		const int NPARAMS ) ; 
int
const_exp_f( const gsl_vector *params , 
	     void *data , 
	     gsl_vector *f ) ;

int
const_exp_df ( const gsl_vector *params , 
	       void *data , 
	       gsl_matrix * J ) ;

int
const_exp_fdf( const gsl_vector *params , 
	       void *data , 
	       gsl_vector *f, 
	       gsl_matrix *J ) ;

#endif
