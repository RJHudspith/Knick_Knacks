#ifndef PPAAWWFIT_H
#define PPAAWWFIT_H

void
PPAA_WW_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
PPAA_WW_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
PPAA_WW_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
PPAA_WW_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f ) ;
int
PPAA_WW_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J ) ;
int
PPAA_WW_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J ) ;

#endif
