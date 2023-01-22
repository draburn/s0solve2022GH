import numpy as np
from numpy.random import default_rng
from scipy import linalg

rngSeed = 0
sizeX = 2
sizeF = 2
print( 'rngSeed = ', rngSeed )
print( 'sizeX = ', sizeX )
print( 'sizeF = ', sizeF )

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

# We could loop over this next part.
matJ = np.zeros( (sizeF,sizeX) )
epsF = 1.0E-4
for n in range( 0, sizeX ):
	vecXP = vecX.copy()
	vecXM = vecX.copy()
	vecXP[n] += epsF
	vecXM[n] -= epsF
	vecFP = f( vecXP )
	vecFM = f( vecXM )
	#print( 'vecXP = ', vecXP )
	#print( 'vecXM = ', vecXM )
	#print( 'vecFP = ', vecFP )
	#print( 'vecFM = ', vecFM )
	matJ[:,n] = ( vecFP - vecFM ) / ( 2.0 * epsF )
print( 'matJ = ', matJ )
if sizeX == sizeF:
	vecDelta = -linalg.solve( matJ, vecF )
else:
	vecDelta = -linalg.lstsq( matJ, vecF )
vecXNext = vecX + vecDelta
vecFNext = f( vecXNext )
print( '||fNext|| = ', np.sqrt( vecFNext @ vecFNext ) )
