Generic 1d differential solver
==============================

Here is a code that will generically solve equations of the form:

dy/dx = f(x,y)

numerically with initial conditions and final integration range in x. Of course this is the same as numerically integrating the function f(x,y) in x.

The idea is to perform many small steps between xmin and xmax using various integration schemes. If we consider the simplest case, the Euler step,

    f(x+h,y) = h*f(x,y)

For small step-size h. I adaptively integrate this guy by a double half-step method:

    f(x+h,y) = h*f(x,y)

    g(x+h/2,y)  = h/2*f(x,y)
    g(x+h,y)    = h/2*g(x+h/2,y+g(x,y))

And compare these two results for this single step, decreasing h until this step satisfies some desired tolerance.

    err = g(x+h,y)-f(x+h,y)

The full integral will always come from g as it is more accurate.

    Int = \sum g(x,y)

Various integrators supported
=============================

Following https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods we support various explicit methods, Euler, Heun's, RK3, RK4 ... etc, which can be derived from their corresponding Butcher Tableau.

We also support several implicit methods, Implicit Midpoint, Radau3 and Radau5, and Gauss-Legendre 4th and 5th order. These methods have an implicit dependence on previous steps that I solve by simple fixed-point iteration to some reasonably-high tolerance.

The adaptive growth and shrinking factors are determined accordingly by the assumption of the order of the integrator used.

Python interface
================

I also have a python interface to the library. For this to work the makefile must compile the .so library that the python code links to. In the python script the function to be integrated is defined:

  def callback(P):
      return -2*P.y + P.x + 4 ;

Which takes a "point" object P that has members P.x and P.y. Otherwise the code works the same as the c-code. This can be tested by running ./INT or the script, they do the same thing and give the same result.