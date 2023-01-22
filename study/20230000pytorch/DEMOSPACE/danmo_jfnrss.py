import numpy as np
from numpy.random import default_rng
from scipy import linalg

rngSeed = 0
sizeX = 3
sizeF = 3
numIter = 50
print( 'rngSeed = ', rngSeed )
print( 'sizeX = ', sizeX )
print( 'sizeF = ', sizeF )
print( 'numIter = ', numIter )

rng = default_rng(rngSeed)
matJCrit = rng.standard_normal( (sizeF,sizeX) )
vecXCrit = rng.standard_normal( sizeX )
def f( x ):
	return matJCrit @ ( x - vecXCrit )
#print( 'matJCrit = ', matJCrit )
#print( 'vecXCrit = ', vecXCrit )

vecX0 = np.zeros( sizeX )
vecF0 = f( vecX0 )
#print( 'vecX0 = ', vecX0 )
#print( 'vecF0 = ', vecF0 )
print( '||f0|| = ', np.sqrt( vecF0 @ vecF0 ) )

vecX = vecX0;
vecF = vecF0;

for iterCount in range( 0, numIter ):
	# Construct a random subspace.
	sizeK = 2
	#print( 'sizeK = ', sizeK )
	matV = linalg.orth( rng.standard_normal( (sizeX,sizeK) ) )
	matW = np.zeros( (sizeF,sizeK) )
	epsF = 1.0E-4
	for k in range( 0, sizeK ):
		vecFP = f( vecX + (epsF * matV[:,k]) )
		vecFM = f( vecX - (epsF * matV[:,k]) )
		matW[:,k] = ( vecFP - vecFM ) / ( 2.0*epsF )
	vecDelta = -matV @ linalg.lstsq( matW, vecF )[0]
	vecX += vecDelta
	vecF = f( vecX )
	if 0 == iterCount % 5:
		print( iterCount, np.sqrt( vecDelta @ vecDelta ), np.sqrt( vecF @ vecF ) )
print( iterCount, np.sqrt( vecDelta @ vecDelta ), np.sqrt( vecF @ vecF ) )
