#ifndef TWOPARAM_H
#define TWOPARAM_H

void alpha_twoparam_set_mu( const double munew ) ;
void alpha_twoparam_set_loops( const int loopsnew ) ;
void alpha_twoparam_set_subscale( const double munew ) ;

void
alpha_twoparam_description( const char *message ,
		   const double *__restrict params ,
		   const struct x_descriptor X ,
		   const int NPARAMS ) ;
void
alpha_twoparam_guess( double *params ,
	     void *data ,
	     const int NPARAMS ) ;
inline double
alpha_twoparam_eval( const double *x ,
		     const struct x_descriptor X ,
		     const int NPARAMS ) ;
int
alpha_twoparam_f( const gsl_vector *x , 
		  void *data , 
		  gsl_vector *f ) ;
int
alpha_twoparam_df( const gsl_vector *x, 
		   void *data, 
		   gsl_matrix *J ) ;
int
alpha_twoparam_fdf (const gsl_vector *x , 
		    void *data ,
		    gsl_vector *f , 
		    gsl_matrix *J ) ;

#endif
