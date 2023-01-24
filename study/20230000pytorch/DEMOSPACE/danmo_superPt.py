import inspect
import numpy as np
from numpy.random import default_rng
#from scipy import linalg

# Init logging.
frame = inspect.currentframe()
def msg( *arguments, **keywords ):
	#print( f"[", __file__, ".", frame.f_lineno, "] ", *arguments, **keywords )
	print( f"[{__file__}.{frame.f_lineno:05d}]", *arguments, **keywords )
def norm( v ):
	return np.sqrt( v @ v )

# Init problem.
rngSeed = 0
sizeX = 3
sizeF = sizeX
msg( 'rngSeed = ', rngSeed )
msg( 'sizeX = ', sizeX )
msg( 'sizeF = ', sizeF )
rng = default_rng(rngSeed)
matA = rng.standard_normal(( sizeF, sizeX ))
matHCrit = np.zeros(( sizeX, sizeX ))
matHCrit[:,:] = matA.T @ matA
vecXCrit = rng.standard_normal(( sizeX ))
fCrit = 10.0
def funcFG( x ):
	d = x - vecXCrit
	g = matHCrit @ d
	f = fCrit + (( d @ g )/2.0)
	return ( f, g )
vecX0 = np.zeros(( sizeX ))
f0, vecG0 = funcFG( vecX0 )
msg( 'f0 = ', f0 )
msg( '||vecG0|| = ', norm(vecG0) )

# Init SGD solver.
fBail = f0 * 1E8
fevalLimit = 100000
learningRate = 0.01
momentumFactor = 0.9
msg( 'fevalLimit = ', fevalLimit )
msg( 'learningRate = ', learningRate )
msg( 'momentumFactor = ', momentumFactor )
fevalCount = 0
vecX = vecX0.copy()
vecP = np.zeros(( sizeX ))

# Init superPt.
numFevalPerSuperPt = 100
superPtLimit = 1000
fTol = f0*1.0E-12
gTol = norm(vecG0)*1.0E-12
msg( 'numFevalPerSuperPt = ', numFevalPerSuperPt)
msg( 'superPtLimit = ', superPtLimit )
msg( 'fTol = ', fTol)
msg( 'gTol = ', gTol )
running_fevalCount = 0
running_fTot = 0.0
running_xtgTot = 0.0
running_vecGTot = np.zeros(( sizeX ))
running_vecXTot = np.zeros(( sizeX ))
superPtCount = 0
vecXSeed = vecX.copy()
vecPSeed = vecP.copy()
vecXHarvest = np.zeros( sizeX )
vecPHarvest = np.zeros( sizeX )
superPt_f = 0.0
superPt_vecG = np.zeros( sizeX )
superPt_vecX = np.zeros( sizeX )

# Main loop.
doMainLoop = True
while doMainLoop:
	# Perform feval
	f, vecG = funcFG( vecX )
	fevalCount += 1
	
	# Update superPt running totals before updating SGD.
	xtg = vecX @ vecG
	running_fevalCount += 1
	running_fTot += f
	running_xtgTot += xtg
	running_vecGTot[:] += vecG[:]
	running_vecXTot[:] += vecX[:]
	
	# Updage SGD.
	vecP[:] = ( momentumFactor * vecP[:] ) - ( learningRate * vecG[:] )
	vecX[:] += vecP[:]
	
	# Check per-feval stop crit.
	if ( f > fBail ):
		msg( "IMPOSED STOP: f > fBail. This strongly indicates divergence." )
		doMainLoop = False
	elif ( fevalLimit > 0 and fevalCount >= fevalLimit ):
		msg( "IMPOSED STOP: fevalCount >= fevalLimit." )
		doMainLoop = False
	# Check elapsed time?.
	# Check for "stop signal on disk"?
	if ( not doMainLoop ):
		break
	
	# Have we finished the super point?
	if ( running_fevalCount < numFevalPerSuperPt ):
		continue
	vecXHarvest[:] = vecX[:]
	vecPHarvest[:] = vecP[:]
	
	# Do super-point analysis.
	superPtCount += 1
	#
	superPt_vecG[:] = running_vecGTot[:] / running_fevalCount
	superPt_vecX[:] = running_vecXTot[:] / running_fevalCount
	superPt_fAvg = running_fTot / running_fevalCount
	superPt_xtgAvg = running_xtgTot / running_fevalCount
	superPt_f = superPt_fAvg - (( superPt_xtgAvg - ( superPt_vecX @ superPt_vecG ) )/2.0)
	#
	running_fevalCount = 0
	running_fTot = 0.0
	running_xtgTot = 0.0
	running_vecGTot[:] = 0.0
	running_vecXTot[:] = 0.0
	
	# Check superPt stop crit.
	if ( norm(superPt_vecG) <= gTol ):
		msg( "SUCCESS: norm(superPt_vecG) <= gTol." )
		doMainLoop = False
	elif ( superPt_f <= fTol ):
		msg( "SUCCESS: superPt_f <= fTol." )
		doMainLoop = False
	elif ( superPtCount > 0 and superPtCount >= superPtLimit ):
		msg( "IMPOSED STOP: superPtCount >= superPtLimit." )
		doMainLoop = False
	if ( not doMainLoop ):
		break
	
	# Print progress log.
	msg( f"  {fevalCount:7d}, {superPtCount:5d};",
	  f"  {norm( superPt_vecX - vecX0 ):8.2E};",
	  f"  {norm( vecXHarvest - vecXSeed ):8.2E};",
	  f"  {superPt_f:8.2E};",
	  f"  {norm(superPt_vecG):8.2E}" )
	
	# Prepare for next iteration.
	# Record seed for posterity.
	vecXSeed[:] = vecX
	vecPSeed[:] = vecP
	

# Look at results.
msg( f"  {fevalCount:7d}, {superPtCount:5d};",
  f"  {norm( superPt_vecX - vecX0 ):8.2E};",
  f"  {norm( vecXHarvest - vecXSeed ):8.2E};",
  f"  {superPt_f:8.2E};",
  f"  {norm(superPt_vecG):8.2E}" )
vecXF = vecX
fF, vecGF = funcFG( vecXF )
msg( '||vecXF - vecX0|| = ', norm( vecXF - vecX0 ) )
msg( '||vecXF - vecXCrit|| = ', norm( vecXF - vecXCrit ) )
msg( 'fF = ', fF )
msg( 'fF - fCrit = ', fF - fCrit )
msg( '||vecGF|| = ', norm(vecGF) )
