/**
   @file 1dsolve.c
   @brief 1d numerical integrators
   Adaptive step size algorithm for solving 1d differential equations
   just in the variable x for the moment, extension is simple
 */
#include "common.h"
#include "integrators.h"

Integrator integrator ;

// does an integration step
// f is a function pointer
// h is the stepsize
static double
step( double (*f)( const point P ) ,
      double *h ,
      const point P )
{
  register double f1 = 0 , f2 = 0 ;
  size_t n = 0 ;

  // the answer
  integrator.error = 42 ;

  // usual shit. Do one full step at x and 
  // do 2 half steps
  while( ( integrator.error ) > 1.0 && 
	 ( n < integrator.maxsteps ) ) {

    // compute the difference between one step and 
    // two half steps to estimate our truncation error
    const double h_2 = *h/2. ;
    f1  = integrator.step_fwd( f , *h , P ) ;
    f2  = integrator.step_fwd( f , h_2 , P ) ;
    // step forward our point
    f2 += integrator.step_fwd( f , h_2 , (point){ P.x+h_2 , P.y+f2 } ) ;

    // compute the scaled difference, should be around1
    integrator.error = fabs( f2 - f1 ) / integrator.tolerance ;

    // if the difference is greater than the tolerance
    // shrink the step
    if( integrator. error > 1.0 ) {
      const double tmp = *h * integrator.adaptive_safety * 
	pow( integrator.error , 
	     integrator.adaptive_shrink ) ;
      const double tol = 0.1 * ( *h ) ;
      *h = ( 0. < tmp ) ?			\
	(tmp > tol ? tmp : tol) :\
	(tmp < tol ? tmp : tol) ;
      integrator.notoksteps++ ;
    }
    n++ ;
  }
  integrator.nsteps += n ;
  // return a NaN to ruin everything
  if( n == 25 ) {
    printf( "Ill convergence!" ) ;
    return sqrt(-1) ;
  }
  // always return the most accurate evaluation
  return f2 ;
}

// perform an integration with initial condition y=ystart
double
integrate( double (*f)( const point P ) ,
	   const double xstart , 
	   const double xend ,
	   const double ystart )
{
  register point P = { .x=xstart , .y = ystart } ;
  double h = 0.1 ; // initial guess
  while( P.x < xend ) {
    // should exit if step is nan
    P.y += step( f , &h , P ) ;
    P.x += h ;
    // grow the step
    if( (integrator.error) > ( integrator.adaptive_errcon ) ) {
      h *= pow( integrator.error , integrator.adaptive_growth ) ;
    } else {
      h *= integrator.adaptive_safety * 5.0 ;
    }
  }
  return P.y + integrator.step_fwd( f , xend-P.x , P ) ;
}

// some function with known solution
static inline double
fxy( const point P )
{
  return -2*P.y + P.x + 4 ;
}

void
query_integrator( const double actual ,
		  const double result )
{
  fprintf( stdout , "Finished in %zu steps\n" , integrator.nsteps ) ;
  fprintf( stdout , "Acceptance :: %f \n" , 
	  100*(integrator.nsteps-integrator.notoksteps)
	  /(double)integrator.nsteps ) ; 
  fprintf( stdout , "Result :: %g || Error :: %e \n" , result , result-actual ) ;
}

// main integrator function
int
main( void )
{
  initialise_integrator( RADAU5 , 1E-5 ) ;

  const double actual = -0.75*exp(-2*0.2)+0.5*0.2+1.75 ;
  const double result = integrate( fxy , 0 , 0.2 , 1 ) ;

  query_integrator( actual , result ) ;
    
  return 0 ;
}
