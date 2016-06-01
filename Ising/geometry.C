#include "Ising.h"

// compute the geometry
static void 
get_mom_2piBZ( int x[ ND ] , 
	       const int i , 
	       const int DIMS )
{
  int mu , subvol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu != DIMS ) {
      x[ mu ] = ( ( i - i % subvol ) / subvol ) % Latt.dims[ mu ] ;
      subvol *= Latt.dims[ mu ] ;
    } else {// set it to 0?
      x[ mu ] = 0 ;
    }
  }
  return ;
}

// generate a site
static int
gen_site( const int x[ ND ] )
{
  int res = x[ ND - 1 ] ;
  int mu ;
  for( mu = ND - 1 ; mu > 0 ; mu-- ) {
    res = Latt.dims[ mu - 1 ] * res + x[ mu - 1 ] ;
  }
  return res ;
}

// This is the new bona-fide generic shifting code
static int 
gen_shift( const int i , 
	   const int dir )
{
  int x[ ND ] ; 
  get_mom_2piBZ( x , i , ND ) ;
  if( dir >= 0 ) {
    x[ dir ] = ( x[ dir ] + 1 ) % Latt.dims[ dir ] ;
  } else {
    register const int numu = -dir - 1 ;
    if( x[ numu ] == 0 ) {
      x[ numu ] = x[ numu ] - 1 + Latt.dims[ numu ] ;
    } else {
      x[ numu ] = ( x[ numu ] - 1 ) % Latt.dims[ numu ] ;
    }
  }
  return gen_site( x ) ; 
}

// initialise navigation
void 
init_navig( struct site *__restrict lat )
{
  int i ; 
  for(  i = 0 ; i < LVOLUME ; i++  )  {
    int mu ;
    for(  mu = 0 ; mu < ND ; mu++  )	{
      lat[i].neighbor[mu] = gen_shift( i , mu ) ; 
      lat[i].neighbor[mu+ND] = gen_shift( i , -mu - 1 ) ;  
    }
  }
  return ;
}
