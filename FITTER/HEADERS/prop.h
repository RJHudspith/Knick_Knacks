#ifndef PROPFIT_H
#define PROPFIT_H

void
prop_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS ) ;
void
prop_guess( double *params ,
	    void *data ,
	    const int NPARAMS ) ;
inline double
prop_eval( const double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;
int
prop_f( const gsl_vector *x , 
	void *data , 
	gsl_vector *f ) ;

int
prop_df( const gsl_vector *x, 
	 void *data, 
	 gsl_matrix *J ) ;

int
prop_fdf( const gsl_vector *x , 
	  void *data ,
          gsl_vector *f , 
	  gsl_matrix *J ) ;

#endif
