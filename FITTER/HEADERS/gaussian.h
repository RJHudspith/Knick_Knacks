#ifndef GAUSSIAN_H
#define GAUSSIAN_H

void
gaussian_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
gaussian_guess( double *__restrict params ,
		void *data ,
		const int NPARAMS ) ;
inline double
gaussian_eval( const double *__restrict x ,
	       const struct x_descriptor X ,
	       const int NPARAMS ) ;
int
gaussian_f( const gsl_vector *__restrict x , 
	    void *data , 
	    gsl_vector *__restrict f ) ;
int
gaussian_df( const gsl_vector *__restrict x, 
	     void *data, 
	     gsl_matrix *__restrict J ) ;
int
gaussian_fdf(const gsl_vector *__restrict x , 
	     void *data ,
	     gsl_vector *__restrict f , 
	     gsl_matrix *__restrict J ) ;
#endif
