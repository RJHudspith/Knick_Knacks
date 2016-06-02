#ifndef DOPQ4_FIT_H
#define DOPQ4_FIT_H

void set_mu_D0pQ4( const double munew );
void set_loops_D0pQ4( const int loopsnew );
void
D0pQ4_description( const char *message ,
		   const double *__restrict params ,
		   const struct x_descriptor X ,
		   const int NPARAMS );
void
D0pQ4_guess( double *__restrict params ,
	     void *data ,
	     const int NPARAMS );
inline double
D0pQ4_eval( const double *__restrict x ,
	    const struct x_descriptor X ,
	    const int NPARAMS );
int
D0pQ4_f( const gsl_vector *__restrict x , 
	 void *data , 
	 gsl_vector *__restrict f );
int
D0pQ4_df( const gsl_vector *__restrict x, 
	  void *data, 
	  gsl_matrix *__restrict J );
int
D0pQ4_fdf (const gsl_vector *__restrict x , 
	   void *data ,
	   gsl_vector *__restrict f , 
	   gsl_matrix *__restrict J );
#endif
