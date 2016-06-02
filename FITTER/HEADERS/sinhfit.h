#ifndef SINHFIT_H
#define SINHFIT_H

void
sinh_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
sinh_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
sinh_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
sinh_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f ) ;
int
sinh_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J ) ;
int
sinh_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J ) ;
#endif
