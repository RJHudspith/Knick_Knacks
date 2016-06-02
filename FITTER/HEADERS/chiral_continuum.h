#ifndef CHIRAL_CONT_H
#define CHIRAL_CONT_H

void
chiral_cont_description( const char *message ,
			const double *__restrict params ,
			const struct x_descriptor X ,
			const int NPARAMS ) ;
void
chiral_cont_guess( double *params ,
		  void *data ,
		  const int NPARAMS ) ;
inline double
chiral_cont_eval( const double *x ,
		 const struct x_descriptor X ,
		 const int NPARAMS ) ;
int
chiral_cont_f( const gsl_vector *x , 
	      void *data , 
	      gsl_vector *f ) ;
int
chiral_cont_df( const gsl_vector *x, 
	       void *data, 
	       gsl_matrix *J ) ;
int
chiral_cont_fdf (const gsl_vector *x , 
		void *data ,
		gsl_vector *f , 
		gsl_matrix *J ) ;

#endif
