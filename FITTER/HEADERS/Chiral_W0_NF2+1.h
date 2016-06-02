#ifndef CHIRAL_W0_NF2p1_H
#define CHIRAL_W0_NF2p1_H

inline double
chiralW0_nf2p1_eval( double *x ,
		     const struct x_descriptor X ,
		     const int NPARAMS ) ;

int
chiralW0_nf2p1_f( const gsl_vector *x , 
		  void *data , 
		  gsl_vector *f ) ;

int
chiralW0_nf2p1_df( const gsl_vector *x, 
		   void *data, 
		   gsl_matrix *J ) ;

int
chiralW0_nf2p1_fdf ( const gsl_vector *x , 
		     void *data ,
		     gsl_vector *f , 
		     gsl_matrix *J ) ;

#endif
