import inspect
import time
import numpy as np
from numpy.linalg import norm

danutil_import_time = time.time()

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
def msgtime():
	msg(f'It is {time.asctime()}; time since start is {time.time()-danutil_import_time:0.3f}s.')
# End msgtime().

# Upper-triangular orthonormalization, with drop.
# DRaburn 2023-02-21:
#  Ideally, we'd also calculate and return the matR of the QR factorization, but, POITROME.
def utorthdrop( matA, dropRelThresh, dropAbsThresh ):
	matV = matA.copy() # Superfluous?
	sizeK = matA.shape[1]
	vecKeep = np.zeros((sizeK), dtype='bool')
	vecKeep[:] = True
	for k in range (sizeK):
		vNorm = norm(matV[:,k])
		if ( vNorm <= dropAbsThresh ):
			vecKeep[k] = False
			matV[:,k] = 0.0
		else:
			matV[:,k] /= vNorm
	for k in range (1, sizeK):
		matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
		vNorm = norm( matV[:,k] )
		if ( vNorm <= dropRelThresh ):
			vecKeep[k] = False
			matV[:,k] = 0.0
		else:
			matV[:,k] /= vNorm
			matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
			vNorm = norm( matV[:,k] )
			if ( vNorm <= dropRelThresh ):
				vecKeep[k] = False
				matV[:,k] = 0.0
			else:
				matV[:,k] /= vNorm
	matV = matV[:,vecKeep]
	return matV, vecKeep
# End utorthdrop().

def symm( matA ):
	return (matA.T + matA)/2.0
# End symm().

def bye():
	print(f'[{inspect.stack()[1].filename}.{inspect.stack()[1].lineno:05d}]')
	print(f'[{inspect.stack()[1].filename}.{inspect.stack()[1].lineno:05d}] It is {time.asctime()}; time since start is {time.time()-danutil_import_time:0.3f}s.')
	print(f'[{inspect.stack()[1].filename}.{inspect.stack()[1].lineno:05d}] Goodbye.')
	exit()
# End bye().

def linishrootOfQuad( a, b, c, tol=1.0E-6 ):
	# Find "linish root" of "y = a*(x^2) + b*x + c";
	# this is the root that corresponds to the a = 0 solution;
	# if there is no root, this returns the extremum.
	if ( 0.0 == b ):
		# We have no "forward" direction.
		if ( a*c >= 0.0 ):
			# We're at the extremum and there's no root.
			return 0.0
		msg('WARNING: Ambiguous direction. Returning positive root.')
		return np.sqrt(abs(c/a))
	discrim = (b*b) - (4.0*a*c)
	if ( discrim < 0.0 ):
		# No root; return extremum.
		return -b/(2.0*a)
	if ( abs(a*c) < tol*(b*b) ):
		# Use near-linear model.
		return -c*( 1.0 - ((a*c)/(b*b)) )/b
	# Use general quadratic model.
	return b*( np.sqrt(1.0 - ((4.0*a*c)/(b*b))) - 1.0 )/(2.0*a)
# End linishrootOfQuad().
