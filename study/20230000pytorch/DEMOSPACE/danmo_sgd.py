import numpy as np
from numpy.random import default_rng
#from scipy import linalg

rngSeed = 0
sizeX = 3
sizeF = sizeX
print( 'rngSeed = ', rngSeed )
print( 'sizeX = ', sizeX )
print( 'sizeF = ', sizeF )
rng = default_rng(rngSeed)
matA = rng.standard_normal(( sizeF, sizeX ))
matHCrit = matA.T @ matA
vecXCrit = rng.standard_normal(( sizeX ))
fCrit = 10.0
def funcFG( x ):
	d = x - vecXCrit
	g = matHCrit @ d
	f = fCrit + (( d @ g )/2.0)
	return ( f, g )
vecX0 = np.zeros(( sizeX ))
f0, vecG0 = funcFG( vecX0 )
print( 'f0 = ', f0 )
print( 'vecG0 = ', vecG0 )

iterLimit = 10000
learningRate = 0.01
momentumFactor = 0.9
print( 'iterLimit = ', iterLimit )
print( 'learningRate = ', learningRate )
print( 'momentumFactor = ', momentumFactor )
vecX = vecX0
vecP = np.zeros(( sizeX ))
for iterCount in range( 0, iterLimit ):
	f, vecG = funcFG( vecX )
	vecP = ( momentumFactor * vecP ) - ( learningRate * vecG )
	vecX += vecP

vecXF = vecX
fF, vecGF = funcFG( vecXF )
print( 'fF = ', fF )
print( 'vecGF = ', vecGF )
