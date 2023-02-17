import inspect
import numpy as np
import danutil
#import scipy
from scipy import optimize as scipy_optimize

MYEPS = 1.0E-8
main_frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{main_frame.f_lineno:05d}]', *arguments, **keywords)

# Upper-triangular orthonormalization, with drop
def utorthdrop( matA, dropRelThresh, dropAbsThresh ):
    matV = matA.copy() # Superfluous?
    sizeK = matA.shape[1]
    #msg( 'sizeK = ', sizeK )
    rvcDrop = np.zeros( (sizeK), dtype='bool' )
    rvcDrop[:] = False # Superfluous?
    #msg( 'matV = \n', matV )
    #msg( 'rvcDrop =', rvcDrop )
    for k in range ( 0, sizeK ):
        vNorm = np.linalg.norm( matV[:,k] )
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
        ###    continue
        #msg( 'k = ', k )
        #msg( 'matV[:,k] = ', matV[:,k] )
        matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
        #msg( 'matV[:,k] = ', matV[:,k] )
        vNorm = np.linalg.norm( matV[:,k] )
        if ( vNorm <= dropRelThresh ):
            rvcDrop[k] = True
            matV[:,k] = 0.0
        else:
            matV[:,k] /= vNorm
            matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
            vNorm = np.linalg.norm( matV[:,k] )
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
    lambdaF = scipy_optimize.bisect( fRes, lambdaLo, lambdaHi )
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
        return np.linalg.norm( zetaOfP(p) ) - sMax
    def dPastMaxOfP( p ):
        return np.linalg.norm( deltaYOfP(p) ) - dMax
    def fTillMinOfP( p ):
        return fOfP( p ) - fMin
    assert np.min(vecS) > 0.0
    assert np.min(vecLambdaCurve) > 0.0
    assert np.linalg.norm(vecPhi) > 0.0
    p1 = 1.0
    if ( sMax > 0.0 ):
        assert sPastMaxOfP(0.0) < 0.0
        if ( sPastMaxOfP(p1) > 0.0 ):
            p1New = scipy_optimize.bisect( sPastMaxOfP, 0.0, p1 )
            p1 = p1New
    elif ( 0.0 == sMax ):
        p1 = 0.0
    if ( dMax > 0.0 ):
        assert dPastMaxOfP(0.0) < 0.0
        if ( dPastMaxOfP(p1) > 0.0 ):
            p1New = scipy_optimize.bisect( dPastMaxOfP, 0.0, p1 )
            p1 = p1New
    elif ( 0.0 == dMax ):
        p1 = 0.0
    # Ugh. Just require fMin.
    assert fTillMinOfP(0.0) > 0.0
    if ( fTillMinOfP(p1) < 0.0 ):
        p1New = scipy_optimize.bisect( fTillMinOfP, 0.0, p1 )
        p1New = p1
    return deltaYOfP( p1 )

class prm:
	dropRelThresh = 0.1
	dropAbsThresh = 1.0E-12
	pBailThresh = 0.1
	#pBailThresh = 0.001

def calcJump( vecXLaunch, vecPLaunch, record_matX, record_matG, record_rvcF, trSize, prm ):
	# Set default return values
	vecXNext = None
	vecPNext = None
	sizeK = 0
	#
	# Startup.
	if ( trSize <= 0.0 ):
		return ( vecXNext, vecPNext, sizeK )
	numRecords = record_matX.shape[1]
	if ( numRecords <= 1 ):
		return ( vecXNext, vecPNext, sizeK )
	sizeX = record_matX.shape[0]
	#
	# Set anchor.
	vecXAnchor = record_matX[:,0].copy()
	vecGAnchor = record_matG[:,0].copy()
	fAnchor = record_rvcF[0,0]
	matD = record_matX[:,0:numRecords].copy() - np.reshape( vecXAnchor, (sizeX,1) ) # Autobroadcast
	matG = record_matG[:,0:numRecords].copy()
	#
	# Calculate basis.
	matQ, rvcKeep = utorthdrop( matD, prm.dropRelThresh, prm.dropAbsThresh )
	matD = matD[:,rvcKeep]
	matG = matG[:,rvcKeep]
	matR = np.triu( matQ.T @ matD ) # Should really be calc'd by utorthdrop, but, POITROME.
	sizeK = matQ.shape[1]
	if ( 0 == sizeK ):
		return ( vecXNext, vecPNext, sizeK )
	matGamma = matQ.T @ matG
	vecGammaAnchor = matQ.T @ vecGAnchor
	#
	# DRaburn 2023-02-16:
	#  This part is new!
	#  Bail if momentum is too outside of subspace!
	vecPIn = matQ.T @ vecPLaunch
	if ( np.linalg.norm(vecPIn) < prm.pBailThresh * np.linalg.norm(vecPLaunch) ):
		#msg(f'Bailing per pBailThresh: {np.linalg.norm(vecPIn)} < {prm.pBailThresh} * {np.linalg.norm(vecPLaunch)}.')
		sizeK = 0
		return ( vecXNext, vecPNext, sizeK )
	#
	# Generate "fit"/"model".
	# DRaburn 2023-02-16:
	#  'This is simplistic but reasonable.
	#   However, see "hessfit.m".'
	vecGammaFit = vecGammaAnchor
	fFit = fAnchor
	matA = np.linalg.solve( matR.T, (matGamma - np.reshape( vecGammaAnchor, (sizeK,1) ) ).T )
	matHFit = ( matA.T + matA )/2.0
	#
	# Apply scaling?
	vecS = np.ones(sizeK)
	matS = np.diag(vecS)
	matSInv = np.diag(1.0/vecS)
	vecGammaScl = matSInv @ vecGammaFit
	matHScl = matSInv @ matHFit @ matSInv
	#
	# Do eigenfactorization.
	vecLambdaC, matPsi = np.linalg.eig( matHScl )
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
	#
	# Apply floor
	lambdaFloor = getLambdaFloor( fFit, vecPhi, vecLambdaOrig, -0.01*fFit )
	vecLambdaMod = vecLambdaOrig.copy()
	for k in range ( 0, vecLambdaMod.shape[0] ):
		if ( vecLambdaMod[k] < lambdaFloor ):
			vecLambdaMod[k] = lambdaFloor
	matHMod = matS @ matPsi @ np.diag(vecLambdaMod) @ (matPsi.T) @ matS
	#
	vecDLaunch = vecXLaunch - vecXAnchor
	vecYLaunch = matQ.T @ vecDLaunch
	vecXPerp = vecDLaunch - ( matQ @ vecYLaunch )
	vecGammaLaunch = vecGammaFit + ( matHMod @ vecYLaunch )
	fLaunch = fFit + ( vecYLaunch @ vecGammaFit ) + (( vecYLaunch @ vecGammaLaunch )/2.0)
	vecT = matQ.T @ vecPLaunch
	vecPPerp = vecPLaunch - ( matQ @ vecT )
	assert np.linalg.norm( vecGammaLaunch ) > 0.0
	coeffPG = ( vecGammaLaunch @ vecT ) / ( vecGammaLaunch @ vecGammaLaunch )
	vecGammaPerp = vecT - ( coeffPG * vecGammaLaunch )
	testDecomp = True
	if ( testDecomp ):
		assert danutil.reldiff( vecXLaunch, vecXAnchor + (matQ @ vecYLaunch) + vecXPerp ) < MYEPS
		assert danutil.reldiff( vecPLaunch, (matQ @ ( (coeffPG*vecGammaLaunch) + vecGammaPerp )) + vecPPerp ) <= MYEPS
	if ( coeffPG > 0.0 ):
		coeffPG = 0.0

	# Calculate step.
	# TODO.
	# Placeholder: take a Newton step without any trust region.
	#vecZ = np.linalg.solve( matHMod, -vecGammaLaunch )
	vecLambdaCurve = vecLambdaMod  # Shallow copy / reference only / DO NOT MODIFY!
	vecLambdaObjf = vecLambdaOrig  # Shallow copy / reference only / DO NOT MODIFY!
	qnj_fMin = -0.01*fFit
	qnj_sMax = -1.0
	qnj_dMax = trSize
	#msg( 'caps = ', qnj_sMax, qnj_dMax )
	vecZ = levsol( fFit, vecPhi, matPsi, vecLambdaCurve, qnj_sMax, vecS, qnj_dMax, vecLambdaObjf, qnj_fMin )
	###vecZ = np.zeros( sizeK )

	# Generate new seed.
	vecYNew = vecYLaunch + vecZ
	vecGammaNew = vecGammaLaunch + ( matHMod @ vecZ )
	fNew = fLaunch + ( vecZ @ vecGammaLaunch ) + (( vecZ @ vecGammaNew )/2.0)
	assert fLaunch > 0.0
	assert np.linalg.norm(vecGammaLaunch/vecS) > 0.0
	alphaF = fNew / fLaunch
	alphaG = np.linalg.norm( vecGammaNew / vecS ) / np.linalg.norm( vecGammaLaunch / vecS )
	vecDelta = matQ @ vecZ
	#
	vecXNext = np.zeros( sizeX )
	vecPNext = np.zeros( sizeX )
	vecXNext[:] = vecXAnchor + ( matQ @ vecYNew ) + ( alphaG * vecXPerp )
	vecPNext[:] = (matQ @ ( (coeffPG*vecGammaNew) + (alphaG*vecGammaPerp) )) + (alphaF*vecPPerp)
	return ( vecXNext, vecPNext, sizeK )
