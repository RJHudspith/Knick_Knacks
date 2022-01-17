## Python interface to the 1d integrator library
import math
import os
import ctypes
from ctypes import *

## mimics the "point" struct that exists in the c-code
class Point(Structure):
    _fields_ = [("x", c_double),
                ("y", c_double)]

## function we plan on integrating
def callback(P):
    return -2*P.y + P.x + 4 ;

## callback function
MyCB = CFUNCTYPE(c_double,Point)
cb = MyCB(callback)

## open shared library in this directory
fun = ctypes.CDLL(os.getcwd()+"/Lib1D.so")
fun.initialise_integrator.argtypes = [c_int,c_double]

## this is the enum of all the integrators we have and an initialisation
( EULER , BACKWARD_EULER , MIDPOINT , IMPLICIT_MIDPOINT ,
  HEUN , RALSTON , RK3 , RK4 , RK4_38 , GAUSS4 , GAUSS5 ,
  RADAU3 , RADAU5 ) = (0,1,2,3,4,5,6,7,8,9,10,11,12)
fun.initialise_integrator(GAUSS5,1E-14)

## check against known result
fun.integrate.argtypes = [MyCB,c_double,c_double,c_double]
fun.integrate.restype = c_double
## integrate arguments are : callback,xstart,xend,ystart
result = fun.integrate( cb , 0 , 0.2 , 1 )

## print result and compare accuracy to known result
fun.query_integrator.argtypes = [c_double,c_double]
fun.query_integrator(result,-0.75*math.exp(-2*0.2)+0.5*0.2+1.75)
