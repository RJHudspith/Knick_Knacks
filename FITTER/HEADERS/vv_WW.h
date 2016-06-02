#ifndef VVWWFIT_H
#define VVWWFIT_H

void
VV_WW_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
VV_WW_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
VV_WW_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
VV_WW_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f ) ;
int
VV_WW_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J ) ;
int
VV_WW_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J ) ;

#endif
