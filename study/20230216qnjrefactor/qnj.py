import inspect
import numpy as np
import danutil as du

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
        ###    continue
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
    elif ( 0.0 == sMax ):
        p1 = 0.0
    if ( dMax > 0.0 ):
        assert dPastMaxOfP(0.0) < 0.0
        if ( dPastMaxOfP(p1) > 0.0 ):
            p1New = optimize.bisect( dPastMaxOfP, 0.0, p1 )
            p1 = p1New
    elif ( 0.0 == dMax ):
        p1 = 0.0
    # Ugh. Just require fMin.
    assert fTillMinOfP(0.0) > 0.0
    if ( fTillMinOfP(p1) < 0.0 ):
        p1New = optimize.bisect( fTillMinOfP, 0.0, p1 )
        p1New = p1
    return deltaYOfP( p1 )

class prm:
	dummy_int = 137
	dummy_float = 0.123456789
	dummy_char = 'q'

def calcJump( vecXLaunch, vecPLaunch, record_matX, record_matG, record_rvcF, prm ):
	msg('Placeholder hack!')
	msg(f'prm.dummy_int = {prm.dummy_int}')
	return ( vecXLaunch, vecPLaunch )
