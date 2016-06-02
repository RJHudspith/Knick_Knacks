#ifndef BKCHIRAL_H
#define BKCHIRAL_H

void
bkchiral_description( const char *message ,
			const double *__restrict params ,
			const struct x_descriptor X ,
			const int NPARAMS ) ;
void
bkchiral_guess( double *params ,
		  void *data ,
		  const int NPARAMS ) ;
inline double
bkchiral_eval( const double *x ,
		 const struct x_descriptor X ,
		 const int NPARAMS ) ;
int
bkchiral_f( const gsl_vector *x , 
	      void *data , 
	      gsl_vector *f ) ;
int
bkchiral_df( const gsl_vector *x, 
	       void *data, 
	       gsl_matrix *J ) ;
int
bkchiral_fdf (const gsl_vector *x , 
		void *data ,
		gsl_vector *f , 
		gsl_matrix *J ) ;

#endif
