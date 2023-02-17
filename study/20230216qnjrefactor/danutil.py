import numpy as np

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
