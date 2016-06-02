/**
   @file gevp.h
   @brief prototype declaration for generalised eigensolver
 */
#ifndef GEVP_H
#define GEVP_H

/**
   @fn int solve_gevp( double *__restrict re_evalues , const double *__restrict A , const double *__restrict B , const int n )
   @brief generalised eigensolver
   @return SUCCESS(!FAILURE) or FAILURE(-1)
 */
int
solve_gevp( double *__restrict re_evalues , 
	    const double *__restrict A , 
	    const double *__restrict B ,
	    const int n ) ;

#endif
