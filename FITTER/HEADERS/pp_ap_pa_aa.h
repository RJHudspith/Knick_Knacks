#ifndef PPAAFIT_H
#define PPAAFIT_H

void
PPAA_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
PPAA_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
PPAAc_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
PPAA_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f ) ;
int
PPAA_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J ) ;
int
PPAA_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J ) ;

#endif
