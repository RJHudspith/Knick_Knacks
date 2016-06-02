#ifndef PADE_H
#define PADE_H

void 
pade_set_nm( const int new_n ,
	     const int new_m ) ;

void 
pade_get_nm( int *new_n ,
	     int *new_m ) ;

void
pade_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
pade_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
// function evaluation
inline double
pade_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;

int
pade_f( const gsl_vector *params , 
	void *data , 
	gsl_vector *f ) ;
int
pade_df ( const gsl_vector *params , 
	  void *data , 
	  gsl_matrix * J ) ;
int
pade_fdf( const gsl_vector *params , 
	  void *data , 
	  gsl_vector *f, 
	  gsl_matrix *J ) ;

#endif
