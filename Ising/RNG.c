#include "Ising.h"

// set up the table
static uint32_t *table ;

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
  table[3] = 0 ; //t < 0 ; is always >= 0 for unsigned int ...
  table[4] += 1411392427 ;
  table[2] = t&2147483647 ;
  return table[0] + table[2] + table[4] ; //x + y + w ;
}

// generates a random double number
double
KISS_dbl( void ) 
{
#ifdef FULL_DOUBLE
  // to remove the type-punning error gcc always gives about strict aliasing
  union {
    uint32_t theInts[2] ;
    uint64_t ulint ;
    double x ;
  } lu2dbl ;
  // pack the bottom and top half of an unsigned long with uints
  lu2dbl.theInts[0] = JKISS32( ) ;
  lu2dbl.theInts[1] = JKISS32( ) ;
  // shift down so that only the mantissa has non-zero entries
  lu2dbl.ulint = ( lu2dbl.ulint >> 12 ) | 0x3FF0000000000000UL ;
  // have to subtract one because of the implicit 2^0+mantissa
  return lu2dbl.x - 1.0 ;
#else
  return JKISS32( ) * 2.3283064365386963e-10 ;
#endif
}
