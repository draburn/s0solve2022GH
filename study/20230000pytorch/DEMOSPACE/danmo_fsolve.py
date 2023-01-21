import numpy as np
from scipy.optimize import fsolve
def demofunc(x):
	return [ x[0] * np.cos(x[1]) - 4.0,
	         x[1] * x[0] - x[1] - 5.0 ]
rootx = fsolve( demofunc, [1, 1])
print( 'rootx = ', rootx )
print( 'rootx[0] = ', (rootx[0]) )
print( 'rootx[1] = ', (rootx[1]) )
print( 'demofunc(rootx) = ', demofunc(rootx) )
