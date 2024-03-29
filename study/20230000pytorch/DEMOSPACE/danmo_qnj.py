import inspect
import numpy as np
from numpy.random import default_rng
from scipy import linalg
from scipy import optimize

# Init logging.
frame = inspect.currentframe()
def msg( *arguments, **keywords ):
	#print( f"[", __file__, ".", frame.f_lineno, "] ", *arguments, **keywords )
	print( f'[{__file__}.{frame.f_lineno:05d}]', *arguments, **keywords )
def reldiff( a, b ):
	sa = np.sum(np.abs(a))
	sb = np.sum(np.abs(b))
	if ( 0.0 == sa and 0.0 == sb ):
		return 0.0
	return ( np.sum(np.abs(a-b)) / ( sa + sb ) )
def utorthdrop( matA, dropRelThresh, dropAbsThresh ):
	#msg( 'Hey hey hey!' )
	matV = matA.copy() # Superfluous?
	sizeK = matA.shape[1]
	#msg( 'sizeK = ', sizeK )
	rvcDrop = np.zeros( (sizeK), dtype='bool' )
	rvcDrop[:] = False # Superfluous?
	#msg( 'matV = \n', matV )
	#msg( 'rvcDrop =', rvcDrop )
	for k in range ( 0, sizeK ):
		vNorm = linalg.norm( matV[:,k] )
		if ( vNorm <= dropAbsThresh ):
			rvcDrop[k] = True
			matV[:,k] = 0.0
		else:
			#msg( f'divding {k} by {vNorm}.' )
			matV[:,k] /= vNorm
	#msg( 'matV = \n', matV )
	#msg( 'rvcDrop =', rvcDrop )
	for k in range ( 1, sizeK ):
		###if ( rvcDrop[k] ):
		###	continue
		#msg( 'k = ', k )
		#msg( 'matV[:,k] = ', matV[:,k] )
		matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
		#msg( 'matV[:,k] = ', matV[:,k] )
		vNorm = linalg.norm( matV[:,k] )
		if ( vNorm <= dropRelThresh ):
			rvcDrop[k] = True
			matV[:,k] = 0.0
		else:
			matV[:,k] /= vNorm
			matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
			vNorm = linalg.norm( matV[:,k] )
			if ( vNorm <= dropRelThresh ):
				rvcDrop[k] = True
				matV[:,k] = 0.0
			else:
				matV[:,k] /= vNorm
				# Note: if dropThresh is too small, may end up keeping more than sizeX vectors.
	#msg( 'matV = \n', matV )
	#msg( 'rvcDrop =', rvcDrop )
	rvcKeep = ~rvcDrop
	matV = matV[:,rvcKeep]
	#msg( 'matV = \n', matV )
	return ( matV, rvcKeep )
def getLambdaFloor( f0, vecPhi, vecLambda, fFloor ):
	def fRes( lambdaF ):
		fRes = f0 - fFloor
		for k in range( 0, vecLambda.shape[0] ):
			fRes -= ((vecPhi[k])**2) / (2.0*max( vecLambda[k], lambdaF ))
		return fRes
	assert f0 > fFloor
	if ( min(vecLambda) > 0.0 ):
		if ( fRes( 0.0 ) >= 0.0 ):
			return 0.0
	lambdaHi = (vecPhi @ vecPhi) / ( 2.0 * ( f0 - fFloor )  )
	#msg( 'lambda: ', lambdaHi )
	if ( lambdaHi > max(vecLambda) ):
		return lambdaHi
	lambdaLo = MYEPS * max(vecLambda)
	#msg( 'lambda: ', lambdaHi, lambdaLo )
	if ( fRes( lambdaLo ) >= 0.0 ):
		return lambdaLo
	lambdaF = optimize.bisect( fRes, lambdaLo, lambdaHi )
	#msg( 'lambda: ', lambdaHi, lambdaLo, lambdaF )
	#msg( 'res = ', fRes(lambdaF) )
	#exit()
	return lambdaF
# Note:
#  "zeta" here was "phi" in 20230106stage\levsol0111.m;
#  "phi" here was "gamma" in 20230106stage\levsol0111.m.
def levsol( f0, vecPhi, matPsi, vecLambdaCurve, sMax, vecS, dMax, vecLambdaObjf, fMin ):
	def zetaOfP( p ):
		if ( p < MYEPS ):
			vecZeta = p * vecPhi / np.min( vecLambdaCurve )
		else:
			mu = np.min( vecLambdaCurve ) * ( (1.0/p) - 1.0 )
			vecZeta = vecPhi / ( vecLambdaCurve + mu )
		return vecZeta
	def deltaYOfP( p ):
		vecZeta = zetaOfP( p )
		vecDeltaY = ( matPsi @ vecZeta ) / vecS
		return vecDeltaY
	def fOfP( p ):
		vecZeta = zetaOfP( p )
		f = f0 - ( vecZeta @ vecPhi ) + (( vecZeta @ ( vecLambdaObjf * vecZeta ))/2.0)
		return f
	def sPastMaxOfP( p ):
		return linalg.norm( zetaOfP(p) ) - sMax
	def dPastMaxOfP( p ):
		return linalg.norm( deltaYOfP(p) ) - dMax
	def fTillMinOfP( p ):
		return fOfP( p ) - fMin
	assert np.min(vecS) > 0.0
	assert np.min(vecLambdaCurve) > 0.0
	assert linalg.norm(vecPhi) > 0.0
	p1 = 1.0
	if ( sMax > 0.0 ):
		assert sPastMaxOfP(0.0) < 0.0
		if ( sPastMaxOfP(p1) > 0.0 ):
			p1New = optimize.bisect( sPastMaxOfP, 0.0, p1 )
			p1 = p1New
	if ( dMax > 0.0 ):
		assert dPastMaxOfP(0.0) < 0.0
		if ( dPastMaxOfP(p1) > 0.0 ):
			p1New = optimize.bisect( dPastMaxOfP, 0.0, p1 )
			p1 = p1New
	# Ugh. Just require fMin.
	assert fTillMinOfP(0.0) > 0.0
	if ( fTillMinOfP(p1) < 0.0 ):
		p1New = optimize.bisect( fTillMinOfP, 0.0, p1 )
		p1New = p1
	return deltaYOfP( p1 )
# Init problem.
MYEPS = 1.0E-8
rngSeed = 0
sizeX = 50
sizeF = sizeX
msg( 'rngSeed = ', rngSeed )
msg( 'sizeX = ', sizeX )
msg( 'sizeF = ', sizeF )
rng = default_rng(rngSeed)
matA = rng.standard_normal(( sizeF, sizeX ))
#matA = np.diag(np.linspace(1.0,sizeX,sizeX))
matHCrit = np.zeros(( sizeX, sizeX ))
matHCrit[:,:] = matA.T @ matA
vecXCrit = rng.standard_normal(( sizeX ))
#vecXCrit = np.ones(( sizeX ))
fCrit = 10.0
noiseX = 1.0E-5
#noiseX = 0.0
def funcFG( x ):
	#d = x - vecXCrit
	d = x - vecXCrit + noiseX*rng.standard_normal(( sizeX ))
	g = matHCrit @ d
	f = fCrit + (( d @ g )/2.0)
	return ( f, g )
vecX0 = np.zeros(( sizeX ))
f0, vecG0 = funcFG( vecX0 )
msg( 'f0 = ', f0 )
msg( '||vecG0|| = ', linalg.norm(vecG0) )

# Init SGD solver.
fBail = f0 * 1E8
fevalLimit = 100000
learningRate = 0.01
learningRate = 0.001
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
gTol = linalg.norm(vecG0)*1.0E-6
xTol = sizeX * 1.0E-12
#fTol = 1.0E-6
#gTol = 1.0E-6
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
useQNJ = True # Unless...
#useQNJ = False
maxSubspaceSize = maxNumRecords
qnj_dropThresh = 0.1
msg( 'useQNJ = ', useQNJ )
msg( 'maxSubspaceSize = ', maxSubspaceSize )
msg( 'qnj_dropThresh = ', qnj_dropThresh )
# Pre-alloc workspaces.
###matD = np.zeros(( sizeX, maxNumRecords ))
###matV = np.zeros(( sizeX, maxSubspaceSize ))
sizeK = 0
qnj_havePrev = False
qnj_sPrev = -1.0
qnj_sMax = 3.0
qnj_sMax_btCoeff = 0.1
qnj_sMax_ftCoeff = 2.0
qnj_dPrev = -1.0
qnj_dMax = -1.0
qnj_dMax_btCoeff = 0.1
qnj_dMax_ftCoeff = 2.0

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
	#print( 'vecX =\n', vecX )
	
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
		if ( superPt_f <= fBestThresh and linalg.norm(superPt_vecG) < linalg.norm(best_vecG) ):
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
	if ( 0 == superPtCount % 10 ):
		if ( newIsMinf ):
			progLogSymbol = '*'
		elif ( newIsBest ):
			progLogSymbol = '.'
		else:
			progLogSymbol = 'X'
		msg(
		  f' {superPtCount:4d} ({badCount:4d}X), {fevalCount:7d}:',
		  f' {sizeK:3d} / {numRecords:3d}:'
		  f'  {linalg.norm( best_vecX - vecX0 ):8.2E};',
		  f'  {linalg.norm( vecXHarvest - vecXSeed ):8.2E}, {qnj_dPrev:8.2E} / {qnj_dMax:8.2E};',
		  f'  {best_f:8.2E};',
		  f'  {linalg.norm(best_vecG):8.2E}',
		  progLogSymbol )
	
	# Check superPt stop crit.
	if ( linalg.norm(superPt_vecG) <= gTol ):
		msg( 'SUCCESS: linalg.norm(superPt_vecG) <= gTol.' )
		doMainLoop = False
	elif ( superPt_f <= fTol ):
		msg( 'SUCCESS: superPt_f <= fTol.' )
		doMainLoop = False
	elif ( superPtCount > 0 and superPtCount >= superPtLimit ):
		msg( 'IMPOSED STOP: superPtCount >= superPtLimit.' )
		doMainLoop = False
	elif ( qnj_havePrev and ( qnj_dPrev < xTol ) and (not newIsBest) ):
		msg( 'IMPOSED STOP: Failed to improve with a QNJ step smaller than xTol.' )
		doMainLoop = False
	# Check superPt_vecX vs prev?
	# Check vecXHarvest vs vecXSeed?
	# Check superPt_f vs prev?
	if ( not doMainLoop ):
		break
	
	# Prepare for next iteration.
	# Record seed for posterity.
	# This will almost certainly be modified if use a quasi-newton jump.
	sizeK = 0
	vecXSeed[:] = vecX[:]
	vecPSeed[:] = vecP[:]
	
	forceBasisGen = True # For comparison to Octave code.
	if ( (not useQNJ) and (not forceBasisGen) ):
		continue
	
	# Add information to records.
	# Always rolling is wasteful. POITROME.
	record_matX = np.roll( record_matX, 1 )
	record_matG = np.roll( record_matG, 1 )
	record_rvcF = np.roll( record_rvcF, 1 )
	# Does this not require unnecessary mem alloc and copy?
	if ( numRecords < maxNumRecords ):
		numRecords += 1
	record_matX[:,0] = superPt_vecX[:]
	record_matG[:,0] = superPt_vecG[:]
	record_rvcF[0,0] = superPt_f
	
	#msg( 'record_matX =\n', record_matX )
	
	# Finally, QNJ!... or not.
	if ( numRecords < 2 ):
		continue
	elif ( not newIsBest ):
		vecX[:] = best_vecXHarvest[:]
		vecP[:] = best_vecPHarvest[:]
		sizeK = 0
		vecXSeed[:] = vecX[:]
		vecPSeed[:] = vecP[:]
		continue
		# This is 'grad-if-bad'.
	
	# Generate basis.
	# DRaburn 2023-01-24: This is a crude two-pass QR method,
	#  involving an unfortunate cap to the number of records before the QR calculation.
	# POITROME
	vecXAnchor = best_vecX # Shallow copy / reference only / DO NOT MODIFY!
	vecGAnchor = best_vecG # Shallow copy / reference only / DO NOT MODIFY!
	fAnchor = best_f
	###matD[:,0:maxNumRecords] = record_matX[:,0:maxNumRecords] - np.reshape( vecXAnchor, (sizeX,1) ) # Autobroadcast.
	#msg( 'numRecords = ', numRecords )
	matD = record_matX[:,0:numRecords].copy() - np.reshape( vecXAnchor, (sizeX,1) ) # Autobroadcast.
	matG = record_matG[:,0:numRecords].copy()
	#msg( 'vecXAnchor = ', vecXAnchor )
	#msg( 'vecGAnchor = ', vecGAnchor )
	#msg( 'fAnchor = ', fAnchor )
	#msg( 'matD =\n', matD )
	#msg( 'matG =\n', matG )
	# We want an equivalent of my Octave "utorthdrop":
	#  construct a basis upper-triangularly, dropping any vectors that are below some threshold in orthogonality.
	matQ, rvcKeep = utorthdrop( matD, qnj_dropThresh, 1.0E-16 )
	matD = matD[:,rvcKeep]
	matG = matG[:,rvcKeep]
	matR = np.triu( matQ.T @ matD )
	sizeK = matQ.shape[1]
	#msg( 'matQ.T @ matQ =\n', matQ.T @ matQ )
	#msg( 'matQ =\n', matQ )
	#msg( 'matR =\n', matR )
	#msg( 'sizeK = ', sizeK )
	#msg( 'D =\n', matD )
	#msg( 'Q*R =\n', matQ @ matR )
	if ( 0 == sizeK ):
		continue
	matGamma = matQ.T @ matG
	vecGammaAnchor = matQ.T @ vecGAnchor
	if ( not useQNJ ):
		continue
	
	# Generate fit.
	# 2023-02-24: This is simplisic but reasonable.
	#  However, see "hessfit.m".
	vecGammaFit = vecGammaAnchor
	fFit = fAnchor
	matA = linalg.solve( matR.T, (matGamma - np.reshape( vecGammaAnchor, (sizeK,1) ) ).T )
	matHFit = ( matA.T + matA )/2.0
	#
	#vecGammaTrue = matQ.T @ matHCrit @ ( vecXAnchor - vecXCrit )
	#matHTrue = matQ.T @ matHCrit @ matQ
	#msg( 'vecGammaFit = ', vecGammaFit )
	#msg( 'vecGammaTrue = ', vecGammaTrue )
	#msg( 'matHFit =\n', matHFit )
	#msg( 'matHTrue = \n', matHTrue )
	#if ( np.sum(np.abs(matHFit-matHTrue)) >= 1.0e-8*( np.sum(np.abs(matHFit)) + np.sum(np.abs(matHTrue)) ) ):
	#	msg( 'ERROR: fit is not true.' )
	#	doMainLoop = False
	#	break
	
	#vecDelta = matQ @ linalg.solve( matHFit, -vecGammaFit )
	#vecX[:] = vecXAnchor + vecDelta
	#vecXSeed[:] = vecX
	#vecPSeed[:] = vecP
	#continue
	
	# Update trust region and scaling.

	if ( qnj_havePrev ):
		if ( newIsBest ):
			qnj_sMax = qnj_sPrev * qnj_sMax_ftCoeff
			qnj_dMax = qnj_dPrev * qnj_dMax_ftCoeff
		else:
			qnj_sMax = qnj_sPrev * qnj_sMax_btCoeff
			qnj_dMax = qnj_dPrev * qnj_dMax_btCoeff
	#msg( 'caps: ', qnj_sMax, qnj_dMax )
	#msg( 'matR =\n', matR )
	vecCap = np.max( np.abs(matR), 1 )
	vecCap[:] += np.sqrt(MYEPS)*np.max(vecCap)
	#msg( 'vecCap =', vecCap )
	vecS = 1.0 / vecCap.copy()
	#msg( 'vecS = ', vecS )
	matS = np.diag(vecS)
	matSInv = np.diag(1.0/vecS)
	vecGammaScl = matSInv @ vecGammaFit
	matHScl = matSInv @ matHFit @ matSInv
	
	# Apply scaling and do eigenfactorization
	vecLambdaC, matPsi = linalg.eig( matHScl )
	for n in range( 0, vecLambdaC.shape[0]):
		assert np.isreal(vecLambdaC[n])
	vecLambdaOrig = np.real( vecLambdaC )
	vecPhi = matPsi.T @ (-vecGammaScl)
	# So, now:
	#  matM = matLambda + mu * matI
	#  vecZ = matSInv @ matPsi @ ( matM \ vecPhi )
	#  s = np.norm( matS * vecZ ) = np.norm( matM \ vecPhi )
	#  vecDelta = matV @ vecZ
	#  d = np.norm( vecDelta ) = np.norm( vecZ )
	#msg( "But... we want to do all of this from LAUNCH not ANCHOR?" )
	
	# Calculate lambdaMod,
	#  lambda perturbed so that Hessian is pos-def and fModMin >= 0.0
	# Note: we might consider doing this from our launch rather than anchor. Oh well.
	#msg( 'vecLambdaOrig = ', vecLambdaOrig )
	doLambdaFloorTest = False
	if (doLambdaFloorTest):
		fFit = 1.0
		vecPhi = np.array([1.0,1.0,1.0])
		#vecLambdaOrig = np.array([10.0,0.01,-1.0])
		vecLambdaOrig = np.array([0.0,0.0,-1.0])
		vecLambdaOrig = np.array([10.0,10.0,10.0])
		vecLambdaOrig = np.array([10.0,10.0,0.0])
		vecLambdaOrig = np.array([10.0,10.0,0.56])
		vecLambdaOrig = np.array([10.0,10.0,0.55])
		lambdaFloor = getLambdaFloor( fFit, vecPhi, vecLambdaOrig, 0.0 )
		msg( 'lambdaFloor =', lambdaFloor )
		vecLambdaMod = vecLambdaOrig.copy()
		for k in range ( 0, vecLambdaMod.shape[0] ):
			if ( vecLambdaMod[k] < lambdaFloor ):
				vecLambdaMod[k] = lambdaFloor
		fRes = fFit - (np.sum( vecPhi * vecPhi / vecLambdaMod )/2.0)
		msg( 'fRes = ', fRes )
		exit()
	lambdaFloor = getLambdaFloor( fFit, vecPhi, vecLambdaOrig, -0.01*fFit )
	vecLambdaMod = vecLambdaOrig.copy()
	for k in range ( 0, vecLambdaMod.shape[0] ):
		if ( vecLambdaMod[k] < lambdaFloor ):
			vecLambdaMod[k] = lambdaFloor
	matHMod = matS @ matPsi @ np.diag(vecLambdaMod) @ (matPsi.T) @ matS
	
	# Decompose "launch".
	vecXLaunch = best_vecXHarvest.copy()
	vecDLaunch = vecXLaunch - vecXAnchor
	vecYLaunch = matQ.T @ vecDLaunch
	vecXPerp = vecDLaunch - ( matQ @ vecYLaunch )
	vecGammaLaunch = vecGammaFit + ( matHMod @ vecYLaunch )
	fLaunch = fFit + ( vecYLaunch @ vecGammaFit ) + (( vecYLaunch @ vecGammaLaunch )/2.0)
	vecPLaunch = best_vecPHarvest.copy()
	vecT = matQ.T @ vecPLaunch
	vecPPerp = vecPLaunch - ( matQ @ vecT )
	assert linalg.norm( vecGammaLaunch ) > 0.0
	coeffPG = ( vecGammaLaunch @ vecT ) / ( vecGammaLaunch @ vecGammaLaunch )
	vecGammaPerp = vecT - ( coeffPG * vecGammaLaunch )
	testDecomp = True
	if ( testDecomp ):
		assert reldiff( vecXLaunch, vecXAnchor + (matQ @ vecYLaunch) + vecXPerp ) < MYEPS
		assert reldiff( vecPLaunch, (matQ @ ( (coeffPG*vecGammaLaunch) + vecGammaPerp )) + vecPPerp ) <= MYEPS
	if ( coeffPG > 0.0 ):
		coeffPG = 0.0
	
	# Calculate step.
	# TODO.
	# Placeholder: take a Newton step without any trust region.
	#vecZ = linalg.solve( matHMod, -vecGammaLaunch )
	vecLambdaCurve = vecLambdaMod  # Shallow copy / reference only / DO NOT MODIFY!
	vecLambdaObjf = vecLambdaOrig  # Shallow copy / reference only / DO NOT MODIFY!
	qnj_fMin = -0.01*fFit
	vecZ = levsol( fFit, vecPhi, matPsi, vecLambdaCurve, qnj_sMax, vecS, qnj_dMax, vecLambdaObjf, qnj_fMin )
	
	# Generate new seed.
	vecYNew = vecYLaunch + vecZ
	vecGammaNew = vecGammaLaunch + ( matHMod @ vecZ )
	fNew = fLaunch + ( vecZ @ vecGammaLaunch ) + (( vecZ @ vecGammaNew )/2.0)
	assert fLaunch > 0.0
	assert linalg.norm(vecGammaLaunch/vecS) > 0.0
	alphaF = fNew / fLaunch
	alphaG = linalg.norm( vecGammaNew / vecS ) / linalg.norm( vecGammaLaunch / vecS )
	vecDelta = matQ @ vecZ
	#
	vecX = vecXAnchor + ( matQ @ vecYNew ) + ( alphaG * vecXPerp )
	vecP = (matQ @ ( (coeffPG*vecGammaNew) + (alphaG*vecGammaPerp) )) + (alphaF*vecPPerp)
	
	vecXSeed[:] = vecX
	vecPSeed[:] = vecP
	qnj_havePrev = True
	qnj_dPrev = linalg.norm( vecDelta )
	qnj_sPrev = linalg.norm( matS @ vecZ )
	
	

# Look at results.
progLogSymbol = 'F'
msg(
  f' {superPtCount:4d} ({badCount:4d}X), {fevalCount:7d}:',
  f' {sizeK:3d} / {numRecords:3d}:'
  f'  {linalg.norm( superPt_vecX - vecX0 ):8.2E};',
  f'  {linalg.norm( vecXHarvest - vecXSeed ):8.2E}, {qnj_dPrev:8.2E} / {qnj_dMax:8.2E};',
  f'  {superPt_f:8.2E};',
  f'  {linalg.norm(superPt_vecG):8.2E}',
  progLogSymbol )
vecXF = best_vecX
vecGF = best_vecG
fF = best_f
#fF, vecGF = funcFG( vecXF )
msg( '||vecXF - vecX0|| = ', linalg.norm( vecXF - vecX0 ) )
msg( '||vecXF - vecXCrit|| = ', linalg.norm( vecXF - vecXCrit ) )
msg( 'fF = ', fF )
msg( 'fF - fCrit = ', fF - fCrit )
msg( '||vecGF|| = ', linalg.norm(vecGF) )
