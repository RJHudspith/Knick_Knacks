#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

double
D0_OPE( const double t , const double a , const int der , const int NLOOPS ) ;

double
D0_OPE_a2( const double t , const double a , const int der , const int NLOOPS ) ;

double
D2_V_OPE( const double t , const double a , const int der , const int NLOOPS , const double ml , const double ms ) ;

double
D4_GG_OPE( const double t , const double a , const int der , const int NLOOPS ) ;

double
D4_Vqq_OPE( const double t , const double a , const int der , const int NLOOPS , const double ml_ll , const double ms_ss ) ;

#endif
