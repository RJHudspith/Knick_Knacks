#ifndef STATICPOT_H
#define STATICPOT_H

void
staticpot_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
staticpot_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
// function evaluation
inline double
staticpot_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;

int
staticpot_f( const gsl_vector *params , 
	void *data , 
	gsl_vector *f ) ;
int
staticpot_df ( const gsl_vector *params , 
	  void *data , 
	  gsl_matrix * J ) ;
int
staticpot_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
