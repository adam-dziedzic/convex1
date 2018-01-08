from math import exp
from math import log
from math import sin
from math import cos
import math
import numpy as np
import os

# bisection gets handle `fun' to the oracle of a function f, performs
# bisection search on interval [MIN,MAX] and returns an 
# eps-suboptimal solution x; i.e. f(x)-f(x^*) <= eps .
# You only need to add three lines to the code:
#   1) upper bound the suboptimality of MID
#   2,3) update the search interval 

def bisection( fun, MIN, MAX, eps=1e-5 ):

  # counting the number of iterations
  counter = 0
  while True:
    counter +=1

    MID = ( MAX + MIN ) / 2

    # Oracle access to the function value and gradient
    gradient, value = fun( MID )
    norm = np.linalg.norm(MAX-MIN)
    #print("gradient: ", gradient, " value: ", value, " norm: ", norm, " min: ", MIN, " max: ", MAX, "mid: ", MID)

    # provide an upper bound for the suboptimality of MID in terms of
    # the magnitude of the gradient and distance from the optimum
    ###############################
    # TODO: suboptimality = ???
    ###############################
    suboptimality = abs(gradient) * norm
    
    if suboptimality <= eps:
      break

    if gradient > 0:
      ###############################
      # TODO: Updating the interval #
      ###############################
      MAX = MID
    else:
      ###############################
      # TODO: Updating the interval #
      ###############################
      MIN = MID
  print( 'Number of Iterations: %d' %counter )
  print( 'Suboptimal point: %1.15f' %MID )
  print( 'Suboptimal value: %1.15f' %value )
  return MID

###############################################################
# weird_func gets scalar x as input and returns f(x) and f'(x)

def weird_func( x ):

  # f(x) = x^4 + 6x^2 + 12*(x-4)*e^(x-1)
  value = pow(x, 4) + 6 * pow(x, 2) + 12 * (x - 4) * exp(x - 1)

  # f'(x) = 4x^3 + 12*x + 12*(x-3)*e^(x-1)
  gradient = 4 * pow(x, 3) + 12 * x + 12 * (x - 3) * exp(x - 1)

  return (gradient,value)

def simple_func( x ):

  # f(x) = x^2 − 3x + 4
  value = pow(x, 2) - 3*x + 4

  # f'(x) = 2x - 3
  gradient = 2*x - 3

  return (gradient,value)

def example_func( x ):

  # f(x) = 100*(1 − x)*(1 − 3*log(x))
  value = 100*(1 - x)*(1 - 3*log(x))

  # f'(x) = 100*(log(x) - 3/x + 2)
  gradient = 100*(3*log(x) - 3/x + 2)

  return (gradient,value)

def complex_func( x ):

  # f(x) = x^4 + x^2 + 3 + (x - 3)*e^(x-2) + 2*sin(x) - x*cos(x)
  value = pow(x,4) + pow(x,2) + 3 + (x - 3)*pow(math.e,x-2) + 2*sin(x) - x*cos(x)

  # f'(x) = 4*x^3 + 2*x + (x - 2)*e^(x-2) + cos(x) + x*sin(x)
  gradient = 4*pow(x,3) + 2*x + (x - 2)*pow(math.e,x-2) + cos(x) + x*sin(x)

  return (gradient,value)

###############################################################

# You may need to define other functions below

print("BISECTIONS")
for funargs in [dict(fun=weird_func, MIN=-5, MAX=5),
                dict(fun=simple_func, MIN=0, MAX=10),
                dict(fun=example_func, MIN=1, MAX=1.3),
                dict(fun=weird_func, MIN=-5, MAX=5, eps=1e-15),
                dict(fun=complex_func, MIN=-1000, MAX=1000, eps=1e-15)]:
  fun = funargs["fun"]
  print()
  print(fun.__name__, ": ")
  optimal_point = bisection(**funargs)
  if "eps" in funargs:
    print("eps: ", funargs["eps"])
  else:
    print("eps: default(1e-5)")
                  
print("Finished the script: ", os.path.basename(__file__))
