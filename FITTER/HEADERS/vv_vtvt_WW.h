#ifndef VV_VTVT_WWFIT_H
#define VV_VTVT_WWFIT_H

void
VV_VTVT_WW_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
VV_VTVT_WW_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
VV_VTVT_WW_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
VV_VTVT_WW_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f ) ;
int
VV_VTVT_WW_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J ) ;
int
VV_VTVT_WW_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J ) ;

#endif
