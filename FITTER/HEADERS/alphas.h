#ifndef ALPHAS_H
#define ALPHAS_H

void set_mu( const double munew ) ;
void set_loops( const int loopsnew ) ;

void
alpha_description( const char *message ,
		   const double *__restrict params ,
		   const struct x_descriptor X ,
		   const int NPARAMS ) ;
void
alpha_guess( double *params ,
	     void *data ,
	     const int NPARAMS ) ;
inline double
alpha_eval( const double *x ,
	    const struct x_descriptor X ,
	    const int NPARAMS ) ;
int
alpha_f( const gsl_vector *x , 
	 void *data , 
	 gsl_vector *f ) ;
int
alpha_df( const gsl_vector *x, 
	  void *data, 
          gsl_matrix *J ) ;
int
alpha_fdf (const gsl_vector *x , 
	   void *data ,
	   gsl_vector *f , 
	   gsl_matrix *J ) ;

#endif
