#ifndef LOG_PLUS_C_H
#define LOG_PLUS_C_H

void
log_plus_c_description( const char *message ,
			const double *__restrict params ,
			const struct x_descriptor X ,
			const int NPARAMS ) ;
void
log_plus_c_guess( double *params ,
		  void *data ,
		  const int NPARAMS ) ;
inline double
log_plus_c_eval( const double *x ,
		 const struct x_descriptor X ,
		 const int NPARAMS ) ;
int
log_plus_c_f( const gsl_vector *x , 
	      void *data , 
	      gsl_vector *f ) ;
int
log_plus_c_df( const gsl_vector *x, 
	       void *data, 
	       gsl_matrix *J ) ;
int
log_plus_c_fdf (const gsl_vector *x , 
		void *data ,
		gsl_vector *f , 
		gsl_matrix *J ) ;
#endif
