import inspect
import numpy as np
from numpy.random import default_rng
#from scipy import linalg

# Init logging.
frame = inspect.currentframe()
def msg( *arguments, **keywords ):
	#print( f"[", __file__, ".", frame.f_lineno, "] ", *arguments, **keywords )
	print( f'[{__file__}.{frame.f_lineno:05d}]', *arguments, **keywords )
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
superPtLimit = 100
fTol = f0*1.0E-12
gTol = norm(vecG0)*1.0E-18
msg( 'numFevalPerSuperPt = ', numFevalPerSuperPt)
msg( 'superPtLimit = ', superPtLimit )
msg( 'fTol = ', fTol)
msg( 'gTol = ', gTol )
running_fevalCount = 0
running_fTot = 0.0
running_fSqTot = 0.0
running_xtgTot = 0.0
running_vecGTot = np.zeros(( sizeX ))
running_vecXTot = np.zeros(( sizeX ))
superPtCount = 0
vecXSeed = vecX.copy()
vecPSeed = vecP.copy()
vecXHarvest = np.zeros(( sizeX ))
vecPHarvest = np.zeros(( sizeX ))
superPt_f = 0.0
superPt_fVar = 0.0
superPt_vecG = np.zeros(( sizeX ))
superPt_vecX = np.zeros(( sizeX ))

# Init minf and best...
#  "best" is superPt with min ||vecG|| sbjt f not too much larger than minf_f.
coeff_best_minf = 1.0
coeff_best_best = 1.0
coeff_best_curr = 1.0
msg( 'coeff_best_minf = ', coeff_best_minf )
msg( 'coeff_best_best = ', coeff_best_best )
msg( 'coeff_best_curr = ', coeff_best_curr )
minf_present = False
minf_f = 0.0
minf_fVar = 0.0
minf_vecG = np.zeros(( sizeX ))
minf_vecX = np.zeros(( sizeX ))
best_present = False
best_f = 0.0
best_fVar = 0.0
best_vecG = np.zeros(( sizeX ))
best_vecX = np.zeros(( sizeX ))
best_vecXHarvest = np.zeros(( sizeX ))
best_vecPHarvest = np.zeros(( sizeX ))
badCount = 0

# Init records.
maxNumRecords = 20
msg( 'maxNumRecords = ', maxNumRecords )
record_matX = np.zeros(( sizeX, maxNumRecords ))
record_matG = np.zeros(( sizeX, maxNumRecords ))
record_rvcF = np.zeros(( 1, maxNumRecords ))
numRecords = 0

# Init QNJ.
useQNJ = True
maxSubspaceSize = maxNumRecords
msg( 'useQNJ = ', useQNJ )
msg( 'maxSubspaceSize = ', maxSubspaceSize )
# Pre-alloc workspaces.
matD = np.zeros(( sizeX, maxNumRecords ))
matV = np.zeros(( sizeX, maxSubspaceSize ))



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
	running_fSqTot += f*f
	running_xtgTot += xtg
	running_vecGTot[:] += vecG[:]
	running_vecXTot[:] += vecX[:]
	
	# Updage SGD.
	vecP[:] = ( momentumFactor * vecP[:] ) - ( learningRate * vecG[:] )
	vecX[:] += vecP[:]
	
	# Check per-feval stop crit.
	if ( f > fBail ):
		msg( 'IMPOSED STOP: f > fBail. This strongly indicates divergence.' )
		doMainLoop = False
	elif ( fevalLimit > 0 and fevalCount >= fevalLimit ):
		msg( 'IMPOSED STOP: fevalCount >= fevalLimit.' )
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
	superPt_fSqVar = (running_fSqTot/running_fevalCount) - (superPt_fAvg**2)
	if ( 0.0 < superPt_fSqVar ):
		superPt_fVar = np.sqrt( superPt_fSqVar )
	else:
		superPt_fVar = 0.0
	# Note that part of fVar is due to f actually varying along the path;
	#  this is not desirable, but is probably acceptable.
	#
	running_fevalCount = 0
	running_fTot = 0.0
	running_xtgTot = 0.0
	running_vecGTot[:] = 0.0
	running_vecXTot[:] = 0.0
	
	# Do minf & best analysis.
	newIsMinf = False # Unless...
	newIsBest = False # Unless...
	if ( not minf_present ):
		newIsMinf = True
		newIsBest = True
	elif ( superPt_f < minf_f ):
		newIsMinf = True
		newIsBest = True
	else:
		fBestThresh = (  minf_f
		  + ( coeff_best_minf * minf_fVar )
		  + ( coeff_best_best * best_fVar )
		  + ( coeff_best_curr * superPt_fVar )  )
		if ( superPt_f <= fBestThresh and norm(superPt_vecG) < norm(best_vecG) ):
			newIsBest = True
	if ( newIsMinf ):
		minf_f = superPt_f
		minf_fVar = superPt_fVar
		minf_vecX[:] = superPt_vecX[:]
		minf_vecG[:] = superPt_vecG[:]
		minf_present = True
	if ( newIsBest ):
		best_f = superPt_f
		best_fVar = superPt_fVar
		best_vecX[:] = superPt_vecX[:]
		best_vecG[:] = superPt_vecG[:]
		best_vecXHarvest[:] = vecXHarvest[:]
		best_vecPHarvest[:] = vecPHarvest[:]
		best_present = True
	else:
		badCount += 1
	
	# Print progress log.
	if ( newIsMinf ):
		progLogSymbol = '*'
	elif ( newIsBest ):
		progLogSymbol = '.'
	else:
		progLogSymbol = 'X'
	msg(
	  f'  {fevalCount:7d}, {superPtCount:5d}:',
	  f'  (X{badCount:4d}), {numRecords:3d};',
	  f'  {norm( best_vecX - vecX0 ):8.2E};',
	  f'  {norm( vecXHarvest - vecXSeed ):8.2E};',
	  f'  {best_f:8.2E};',
	  f'  {norm(best_vecG):8.2E}',
	  progLogSymbol )
	
	# Check superPt stop crit.
	if ( norm(superPt_vecG) <= gTol ):
		msg( 'SUCCESS: norm(superPt_vecG) <= gTol.' )
		doMainLoop = False
	elif ( superPt_f <= fTol ):
		msg( 'SUCCESS: superPt_f <= fTol.' )
		doMainLoop = False
	elif ( superPtCount > 0 and superPtCount >= superPtLimit ):
		msg( 'IMPOSED STOP: superPtCount >= superPtLimit.' )
		doMainLoop = False
	# Check superPt_vecX vs prev?
	# Check vecXHarvest vs vecXSeed?
	# Check superPt_f vs prev?
	if ( not doMainLoop ):
		break
	
	# Prepare for next iteration.
	# Record seed for posterity.
	# This will almost certainly be modified if use a quasi-newton jump.
	vecXSeed[:] = vecX
	vecPSeed[:] = vecP
	if ( not useQNJ ):
		continue
	
	# Add information to records.
	if ( numRecords == maxNumRecords ):
		record_matX = np.roll( record_matX, 1 )
		record_matG = np.roll( record_matG, 1 )
		record_rvcF = np.roll( record_rvcF, 1 )
		# Does this not require unnecessary mem alloc and copy?
	else:
		numRecords += 1
	record_matX[:,0] = superPt_vecX[:]
	record_matG[:,0] = superPt_vecG[:]
	record_rvcF[0,0] = superPt_f
	
	# Finally, QNJ!
	if ( numRecords < 2 ):
		continue
	
	# Generate basis.
	# DRaburn 2023-01-24: This should be as good as anything.
	vecXAnchor = best_vecX # Shallow copy / reference only / DO NOT MODIFY!
	vecGAnchor = best_vecG # Shallow copy / reference only / DO NOT MODIFY!
	fAnchor = best_f
	matD[:,0:maxNumRecords] = record_matX[:,0:maxNumRecords] - np.reshape( vecXAnchor, (sizeX,1) ) # Autobroadcast.
	
	msg( 'HACK!' )
	msg( 'WE NEED SOME EQUIVALENT TO UTORTHDROP HERE.' )
	doMainLoop = False
	break
	
	

# Look at results.
progLogSymbol = 'F'
msg(
  f'  {fevalCount:7d}, {superPtCount:5d}:',
  f'  (X{badCount:4d}), {numRecords:3d};',
  f'  {norm( superPt_vecX - vecX0 ):8.2E};',
  f'  {-1.0:8.1E};',
  f'  {superPt_f:8.2E};',
  f'  {norm(superPt_vecG):8.2E}',
  progLogSymbol )
vecXF = best_vecX
vecGF = best_vecG
fF = best_f
#fF, vecGF = funcFG( vecXF )
msg( '||vecXF - vecX0|| = ', norm( vecXF - vecX0 ) )
msg( '||vecXF - vecXCrit|| = ', norm( vecXF - vecXCrit ) )
msg( 'fF = ', fF )
msg( 'fF - fCrit = ', fF - fCrit )
msg( '||vecGF|| = ', norm(vecGF) )
