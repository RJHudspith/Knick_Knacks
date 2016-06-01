/**
   @brief a genetic algorithm
 */
#include <stdint.h> // for uint32_t
#include <stdlib.h> // for malloc
#include <string.h>
#include <stdio.h>

#define Ngen (150)
#define poly_degree (4)

// set up the table
static uint32_t *table ;

struct genes {
  float P[ poly_degree ] ;
} ;

static int
insertion_sort( void *list1 , void *list2 , 
		const size_t base1 , const size_t base2 ,
		const size_t lower , const size_t upper ,
		int(*compare)( const void *a , const void *b ) )
{
  char insert1[ base1 ] , insert2[ base2 ] ; // temporary storage
  int i , hole ;
  // standard insertion sort
  for( i = lower ; i < upper ; i++ ) {
    memcpy( insert1 , (char*)list1 + i*base1 , base1 ) ;
    memcpy( insert2 , (char*)list2 + i*base2 , base2 ) ;
    hole = i ;
    while( ( hole > 0 ) && ( compare( insert1 , (char*)list1 + base1*(hole - 1) ) ) ) {
      memcpy( (char*)list1 + base1*hole , (char*)list1 + base1*(hole - 1) , base1 ) ;
      memcpy( (char*)list2 + base2*hole , (char*)list2 + base2*(hole - 1) , base2 ) ;
      hole-- ;
    }
    memcpy( (char*)list1 + base1*hole , insert1 , base1 ) ;
    memcpy( (char*)list2 + base2*hole , insert2 , base2 ) ;
  }
  return 0 ;
}

int
lt_dbl( const void *a , 
	const void *b )
{
  return (*(const double*)a < *(const double*)b) ; 
}

// set the allocated table
void
GLU_set_KISS_table( const uint32_t seed )
{
  table = ( uint32_t* )malloc( 5 * sizeof( uint32_t ) ) ;
  table[0] = seed ;
  int i ;
  for( i = 1 ; i < 5 ; i++ ) {
    table[i] = ( 1812433253UL * ( table[i-1] ^ ( table[i-1] >> 30)) + i ) ;
  }
  return ;
}

// free the allocated table
void
GLU_free_KISS_table( void )
{
  free( table ) ;
}

// unsigned integer version
static uint32_t
JKISS32( void )
{
  table[0] ^= ( table[0] << 5 ) ;
  table[0] ^= ( table[0] >> 7 ) ;
  table[0] ^= ( table[0] << 22 ) ;
  const uint32_t t = table[1] + table[2] + table[3] ;
  table[1] = table[2] ;
  table[3] = 0 ; //t < 0 ; is always >= 0. Unsigned!
  table[4] += 1411392427 ;
  table[2] = t&2147483647 ;
  return table[0] + table[2] + table[4] ; //x + y + w ;
}

// generates a random double number
double
KISS_dbl( void ) 
{
  return JKISS32( ) * 2.3283064365386963e-10 ;
}

// give a random double within some upper and lower bound
double
KISS_range( const double upper , 
	    const double lower )
{
  return (upper - lower)*KISS_dbl() + lower ;
}

// fitness function is a chi-square?
static double
poly_eval( const float *params ,
	   const float x )
{
  size_t i ;
  register double sum ;
  switch( poly_degree ) {
  case 1 :
    return params[0] ;
  case 2 :
    return params[0] + params[1] * x ;
  case 3 : 
    return params[0] + x * ( params[1] + x * params[2] ) ;
  default :
    // horners rule evaluation
    sum = x * params[ poly_degree-1 ] ;
    for( i = poly_degree - 2 ; i > 0 ; i-- ) {
      sum = x * ( params[i] + sum ) ;
    }
    return sum + params[0] ;
  }
}

// ( y(x) - f(x) )^2 = 0
static double
poly_chisq( const float *params ,
	    const float *y ,
	    const float *x ,
	    const size_t Ndata )
{
  register double sumsq = 0.0 ;
  register float cache ;
  // evaluate function
  size_t i ;
  for( i = 0 ; i < Ndata ; i++ ) {
    cache = ( y[ i ] - poly_eval( params , x[i] ) ) ;
    sumsq += cache * cache ;
  }
  return sumsq ;
}

// compute the fitness of our whole population
// what we do is compare the fitness for each element of the 
// poly with the best values of the other coefficients
static double
fitness( struct genes *gene_pool , // is [Ngen][Npoly] matrix
	 const float *y ,
	 const float *x ,
	 const size_t Ndata ,
	 const double range_upper , 
	 const double range_lower )
{
  // collect a measure of the fitnesses for our various generations
  double fitnesses[ Ngen ] = {} ;

  struct genes new_pool[ Ngen ] ;

  // more efficient to compare all possible variants against one-another
  // and then take the best NGen? is this like a crossover anyway?
  size_t i , j , mu , nu ;
  for( i = 0 ; i < Ngen ; i++ ) { // fix the first parameter of the loop

    // so the trial solution is the first in the list
    float best[ poly_degree ] ;
    for( nu = 0 ; nu < poly_degree ; nu++ ) {
      best[ nu ] = gene_pool[ i ].P[ nu ] ;
    }

    // compute our trial best fitness
    fitnesses[i] = poly_chisq( best , y , x , Ndata ) ;

    // compare to all others getting the best fitness for this parameter
    // loop all possible variants for each coefficient
    float trial[ poly_degree ] , trial_chisq = 1000 ;
    memcpy( trial , best , poly_degree * sizeof( float ) ) ;
    for( nu = 1 ; nu < poly_degree ; nu++ ) {
      for( j = 0 ; j < Ngen ; j++ ) {
	trial[nu] = gene_pool[j].P[nu] ;
	trial_chisq = poly_chisq( trial , y , x , Ndata ) ;
	if( trial_chisq < fitnesses[i] ) {
	  best[nu] = trial[nu] ;
	  fitnesses[i] = trial_chisq ;
	}
      }
    }

    // store the best mixture of our fit coeffiecients in new pool
    for( nu = 0 ; nu < poly_degree ; nu++ ) {
      new_pool[ i ].P[ nu ] = best[ nu ] ;
    }
  }

  // sort the new pool by fitness
  insertion_sort( fitnesses , new_pool , 
		  sizeof( double ) , sizeof( struct genes ) ,
		  0 , Ngen , lt_dbl ) ;

  // top some amount persist
  {
    const size_t cross = (size_t)(Ngen*0.15) ;
    const size_t mutate = (size_t)(Ngen*0.75) ;
    for( i = 0 ; i < cross ; i++ ) {
      memcpy( gene_pool[i].P , new_pool[i].P , poly_degree * sizeof( float ) ) ;
    }
    // next set are crossed over as an average of two parents
    for( i = cross ; i < mutate ; i++ ) {
      // randomly select two parents
      size_t dad = (size_t)KISS_dbl()*cross ;
      size_t mum = (size_t)KISS_dbl()*cross ;
      for( mu = 0 ; mu < poly_degree ; mu++ ) {
	// crossover is an average of the two parameters
	gene_pool[i].P[mu] = 0.5 * ( new_pool[dad].P[mu] + 
				     new_pool[mum].P[mu] ) ;
      }
    }
    // final set are mutations of the remaining population?
    for( i = cross ; i < Ngen ; i++ ) {
      for( mu = 0 ; mu < poly_degree ; mu++ ) {
	gene_pool[i].P[mu] = KISS_range( range_upper , range_lower ) ;
      }
      //
    }
  }

  return fitnesses[0] ; // sorted so is the best fitness
}

// set up our pool with some random shit
static void
initialise_pool( struct genes *gene_pool ,
		 const double range_upper ,
		 const double range_lower )
{
  size_t i , mu ;
  for( i = 0 ; i < Ngen ; i++ ) {
    for( mu = 0 ; mu < poly_degree ; mu++ ) {
      gene_pool[i].P[mu] = KISS_range( range_upper , range_lower ) ;
    }
  }
  return ;
}

// main prog
int main( int argc , 
	  char *argv[] )
{
  GLU_set_KISS_table( 12345 ) ;

  // initialise our gene pool
  struct genes *gene_pool = malloc( Ngen * sizeof( struct genes ) ) ;
  size_t mu ;
  const double lower = -1.1 ;
  const double upper = +1.1 ;
  initialise_pool( gene_pool , upper , lower ) ;

  float xdata[ 6 ] = { 0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 } ;
  float ydata[ 6 ] = { 0.5 , 0.51 , 0.54 , 0.59 , 0.66 , 0.75 } ;

  float chisq = 1 ;
  size_t iters = 0 ;
  while( chisq > 1E-6 && iters < 5000 ) {
    chisq = fitness( gene_pool , ydata , xdata , 6 , upper , lower ) ;
    iters++ ;
  }
  printf( "Finished after %zu iterations \n" , iters ) ;

  // print out our most selfish gene
  printf( "GENE :: %f" , gene_pool[0].P[0] ) ; 
  for( mu = 1 ; mu < poly_degree ; mu++ ) {
    printf( " + %f.x^%zu" , gene_pool[0].P[mu] , mu ) ; 
  }
  printf( " :: fitness %e \n" , chisq ) ;

  GLU_free_KISS_table(  ) ;

  return 0 ;
}
