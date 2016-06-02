#ifndef LAM_DO_H
#define LAM_DO_H

void set_mu_LAM_D0( const double munew );
void set_loops_LAM_D0( const int loopsnew );
void init_betas( const int nf ) ;
double
lambda_MS( const double alpha , const double mu , const int loops ) ;
double
alpha_Lambda( const double mu , const double Lambda , const int loops );

void
LAM_D0_description( const char *message ,
		    const double *__restrict params ,
		    const struct x_descriptor X ,
		    const int NPARAMS );
void
LAM_D0_guess( double *__restrict params ,
	      void *data ,
	      const int NPARAMS );
inline double
LAM_D0_eval( const double *__restrict x ,
	     const struct x_descriptor X ,
	     const int NPARAMS );
int
LAM_D0_f( const gsl_vector *__restrict x , 
	  void *data , 
	  gsl_vector *__restrict f );
int
LAM_D0_df( const gsl_vector *__restrict x, 
	   void *data, 
	   gsl_matrix *__restrict J );
int
LAM_D0_fdf (const gsl_vector *__restrict x , 
	    void *data ,
	    gsl_vector *__restrict f , 
	    gsl_matrix *__restrict J );
#endif
