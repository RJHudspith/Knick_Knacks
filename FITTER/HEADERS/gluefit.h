#ifndef GLUEFIT_H
#define GLUEFIT_H

double
glue_eval( double *x ,
	   const struct x_descriptor X ,
	   const int NPARAMS ) ;

int
glue_f( const gsl_vector *x , 
	void *data , 
	gsl_vector *f ) ;

int
glue_df( const gsl_vector *x, 
	 void *data, 
	 gsl_matrix *J ) ;

int
glue_fdf (const gsl_vector *x , 
	  void *data ,
          gsl_vector *f , 
	  gsl_matrix *J ) ;

#endif
