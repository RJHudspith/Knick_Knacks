#ifndef DO_DIFF_H
#define DO_DIFF_H

int D0_diff_n( void ) ;

void set_mu_D0_diff( const double munew );

void set_loops_D0_diff( const int loopsnew );

void set_q2q3_D0_diff( const double Q2 , 
		       const int NSIM ) ;

void
D0_diff_description( const char *message ,
		const double *__restrict params ,
		const struct x_descriptor X ,
		const int NPARAMS );
void
D0_diff_guess( double *__restrict params ,
	  void *data ,
	  const int NPARAMS );
inline double
D0_diff_eval( const double *__restrict x ,
	 const struct x_descriptor X ,
	 const int NPARAMS );
int
D0_diff_f( const gsl_vector *__restrict x , 
      void *data , 
      gsl_vector *__restrict f );
int
D0_diff_df( const gsl_vector *__restrict x, 
       void *data, 
       gsl_matrix *__restrict J );
int
D0_diff_fdf (const gsl_vector *__restrict x , 
	void *data ,
	gsl_vector *__restrict f , 
	gsl_matrix *__restrict J );
#endif
