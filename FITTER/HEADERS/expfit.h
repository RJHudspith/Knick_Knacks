#ifndef EXPFIT_H
#define EXPFIT_H

void
exp_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
exp_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
// function evaluation
inline double
exp_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ; 
int
exp_f( const gsl_vector *params , 
        void *data , 
        gsl_vector *f ) ;

int
exp_df ( const gsl_vector *params , 
	  void *data , 
          gsl_matrix * J ) ;

int
exp_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
