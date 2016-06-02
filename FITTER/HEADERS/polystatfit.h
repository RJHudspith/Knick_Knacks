#ifndef POLYSTAT_H
#define POLYSTAT_H

void
polystat_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
polystat_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
// function evaluation
inline double
polystat_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;

int
polystat_f( const gsl_vector *params , 
	void *data , 
	gsl_vector *f ) ;
int
polystat_df ( const gsl_vector *params , 
	  void *data , 
	  gsl_matrix * J ) ;
int
polystat_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
