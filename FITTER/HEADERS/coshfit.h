#ifndef COSHFIT_H
#define COSHFIT_H

void
cosh_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
cosh_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
cosh_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
cosh_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f ) ;
int
cosh_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J ) ;
int
cosh_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J ) ;

#endif
