#ifndef CHIRAL_R23_H
#define CHIRAL_R23_H

void
chiral_R23_setfac( const double factor ) ;

void
chiral_R23_description( const char *message ,
			     const double *__restrict params ,
			     const struct x_descriptor X ,
			     const int NPARAMS ) ;
void
chiral_R23_guess( double *params ,
		       void *data ,
		       const int NPARAMS ) ;
inline double
chiral_R23_eval( const double *x ,
		      const struct x_descriptor X ,
		      const int NPARAMS ) ;
int
chiral_R23_f( const gsl_vector *x , 
		   void *data , 
		   gsl_vector *f ) ;
int
chiral_R23_df( const gsl_vector *x, 
		    void *data, 
		    gsl_matrix *J ) ;
int
chiral_R23_fdf (const gsl_vector *x , 
		     void *data ,
		     gsl_vector *f , 
		     gsl_matrix *J ) ;

#endif
