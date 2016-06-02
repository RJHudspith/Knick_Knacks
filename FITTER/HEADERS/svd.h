#ifndef GLUFIT_SVD_H
#define GLUFIT_SVD_H

/**
   @fn int svd_inverse( double *__restrict *__restrict Ainv , const double *__restrict *__restrict A , const int NCOLS , const int NROWS )
   @brief computes the (pseudo)inverse of the matrix A put into Ainv

   @return 0 (failure) or !0 (failure)
 */
int 
svd_inverse( double *__restrict *__restrict Ainv , 
	     const double *__restrict *__restrict A ,
	     const int NCOLS ,
	     const int NROWS ) ;

/**
   @fn int compute_coefficients( double *__restrict coeffs , double *__restrict chisq , const double *__restrict y , const double *__restrict sigma , const double *__restrict x , const double M , const double N )
   @brief compute the coefficients of a polynomial by svd
 */
int
compute_coefficients( double *__restrict coeffs ,
		      double *__restrict chisq ,
		      const double *__restrict y ,
		      const double *__restrict sigma , 
		      const double *__restrict x ,
		      const double M ,
		      const double N ) ;

/**
   @fn int pades( double *pade_coeffs , const double *poly_coeffs , const int n , const int m )
   @brief compute pade coefficients of an n,m pade from minimal polynomial coefficients
 */
int
pades( double *pade_coeffs ,
       const double *poly_coeffs ,
       const int n ,
       const int m ) ;

/**
   @fn void write_polynomial( const double *__restrict coeffs , const int POLY_ORDER )
   @brief write out a polynomial
 */
void
write_polynomial( const double *__restrict coeffs ,
		  const int POLY_ORDER ) ;

#endif
