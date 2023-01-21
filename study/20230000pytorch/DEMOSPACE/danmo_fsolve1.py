import numpy as np
from scipy.optimize import fsolve
print( 'Hello world!' )
def demofunc(x):
	return ( x - np.pi )**2
rootx = fsolve( demofunc, 1.0 )
print( 'rootx = ', rootx )
print( 'demofunc(rootx) = ', demofunc(rootx) )
