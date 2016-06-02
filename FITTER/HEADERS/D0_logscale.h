#ifndef DO_LOGSCALE_H
#define DO_LOGSCALE_H

void set_mu_D0_logscale( const double munew );
void set_loops_D0_logscale( const int loopsnew );
void
D0_logscale_description( const char *message ,
		const double *__restrict params ,
		const struct x_descriptor X ,
		const int NPARAMS );
void
D0_logscale_guess( double *__restrict params ,
	  void *data ,
	  const int NPARAMS );
inline double
D0_logscale_eval( const double *__restrict x ,
	 const struct x_descriptor X ,
	 const int NPARAMS );
int
D0_logscale_f( const gsl_vector *__restrict x , 
      void *data , 
      gsl_vector *__restrict f );
int
D0_logscale_df( const gsl_vector *__restrict x, 
       void *data, 
       gsl_matrix *__restrict J );
int
D0_logscale_fdf (const gsl_vector *__restrict x , 
	void *data ,
	gsl_vector *__restrict f , 
	gsl_matrix *__restrict J );
#endif
