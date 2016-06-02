#ifndef MEXPFIT_H
#define MEXPFIT_H

int getm( void ) ;

void
mexp_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
mexp_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;

// function evaluation
inline double
mexp_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ; 

int
mexp_f( const gsl_vector *params , 
        void *data , 
        gsl_vector *f ) ;

int
mexp_df ( const gsl_vector *params , 
	  void *data , 
          gsl_matrix * J ) ;

int
mexp_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
