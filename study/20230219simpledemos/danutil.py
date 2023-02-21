import numpy as np
import inspect

def var( xAvg, xSqAvg ):
	assert ( type(xAvg) == type(xSqAvg) )
	if ( type(xAvg) == float ):
		xVarSq = xSqAvg - (xAvg**2)
		if ( xVarSq > 0.0 ):
			return np.sqrt(xVarSq)
		return 0.0
	assert( type(xAvg) == np.ndarray )
	xVarSq = xSqAvg - (xAvg**2)
	xVarSq[xVarSq<0.0] = 0.0
	return np.linalg.norm(xVarSq)

def tfInt( tf ):
	if (tf):
		return 1
	else:
		return 0

def reldiff( a, b ):
    sa = np.sum(np.abs(a))
    sb = np.sum(np.abs(b))
    if ( 0.0 == sa and 0.0 == sb ):
        return 0.0
    return ( np.sum(np.abs(a-b)) / ( sa + sb ) )

def msg(*arguments, **keywords):
	print(f'[{inspect.stack()[1].filename}.{inspect.stack()[1].lineno:05d}]', *arguments, **keywords)
