#ifndef DO_DIFF_V1_H
#define DO_DIFF_V1_H

int D0_diff_v1_n( void ) ;

void set_mu_D0_diff_v1( const double munew );

void set_loops_D0_diff_v1( const int loopsnew );

void set_q2q3_D0_diff_v1( const double Q2 , 
		       const int NSIM ) ;

void
D0_diff_v1_description( const char *message ,
		const double *__restrict params ,
		const struct x_descriptor X ,
		const int NPARAMS );
void
D0_diff_v1_guess( double *__restrict params ,
	  void *data ,
	  const int NPARAMS );
inline double
D0_diff_v1_eval( const double *__restrict x ,
	 const struct x_descriptor X ,
	 const int NPARAMS );
int
D0_diff_v1_f( const gsl_vector *__restrict x , 
      void *data , 
      gsl_vector *__restrict f );
int
D0_diff_v1_df( const gsl_vector *__restrict x, 
       void *data, 
       gsl_matrix *__restrict J );
int
D0_diff_v1_fdf (const gsl_vector *__restrict x , 
	void *data ,
	gsl_vector *__restrict f , 
	gsl_matrix *__restrict J );
#endif
