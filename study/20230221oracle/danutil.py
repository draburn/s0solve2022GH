import inspect
import numpy as np
from numpy.linalg import norm

def var( xAvg, xSqAvg ):
	assert (type(xAvg) == type(xSqAvg) )
	if ( type(xAvg) == np.ndarray):
		xVarSq = xSqAvg - (xAvg**2)
		xVarSq[xVarSq<0.0] = 0.0
		return np.sqrt(xVarSq)
	xVarSq = xSqAvg - (xAvg**2)
	if ( xVarSq > 0.0 ):
		return np.sqrt(xVarSq)
	return 0.0
# End def var().

def tfInt( tf ):
	if (tf):
		return 1
	else:
		return 0
# End def tfInt().

def reldiff( a, b ):
    sa = np.sum(np.abs(a))
    sb = np.sum(np.abs(b))
    if ( 0.0 == sa and 0.0 == sb ):
        return 0.0
    return ( np.sum(np.abs(a-b)) / ( sa + sb ) )
# End def reldiff().

def msg( *arguments, **keywords ):
	print(f'[{inspect.stack()[1].filename}.{inspect.stack()[1].lineno:05d}]', *arguments, **keywords)
# End def msg().
