import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm
from scipy import optimize

class hessModelType():
	def __init(self):
		# Basis.
		self.vecXA = None
		self.matV = None
		# Core model.
		self.fA = None
		self.vecGammaA = None # "gamma" = V^T * g
		self.matH = None
		# Perpendicular compnent of gradient.
		self.vecGPerpA = None
		self.matW = None
		# A potential optimization is eliminating calculation of matW.
	def dump(self):
		msg(f'Begin hessModelType().dump...')
		msg(f'self = {self}')
		msg(f'  vecXA = {self.vecXA}')
		msg(f'  matV.shape[0] = {self.matV.shape[0]}')
		msg(f'  matV.shape[1] = {self.matV.shape[1]}')
		msg(f'  matV = ...\n{self.matV}')
		msg(f'  fA = {self.fA}')
		msg(f'  vecGammaA = {self.vecGammaA}')
		msg(f'  matH = ...\n{self.matH}')
		msg(f'  vecGPerpA = {self.vecGPerpA}')
		msg(f'  matW = ...\n{self.matW}')
		msg(f'End hessModelType.dump().')
	def copy(self):
		cp = hessModelType()
		cp.vecXA = self.vecXA.copy()
		cp.matV = self.matV.copy()
		cp.fA = self.fA
		cp.vecGammaA = self.vecGammaA.copy()
		cp.matH = self.matH.copy()
		cp.vecGPerpA = self.vecGPerpA.copy()
		cp.matW = self.matW.copy()
		return cp
	# Evaluation functions are largely for illustrative purpose.
	def evalYOfProjX(self, vecX):
		vecY = self.matV.T @ (vecX - self.vecXA)
		return vecY
	def evalXOfY(self, vecY):
		vecX = self.vecXA + (self.matV @ vecY)
		return vecX
	def evalFGammaGPerpOfY(self, vecY):
		vecHY = self.matH @ vecY
		f = self.fA + (vecY @ self.vecGammaA) + ((vecY @ vecHY)/2.0)
		vecGamma = self.vecGammaA + vecHY
		vecGPerp = self.vecGPerpA + (self.matW @ vecY)
		return f, vecGamma, vecGPerp
	def evalFGOfProjX(self, vecX):
		vecY = self.evalYOfProjX(vecX)
		f, vecGamma, vecGPerp = self.evalFGammaGPerpOfY(vecY)
		vecG = (self.matV @ vecGamma) + vecGPerp
		return f, vecG
# End class hessModelType().

class hessCurvesType():
	def __init__(self):
		# Basis.
		self.vecXA = None
		self.matV = None
		self.vecYLaunch = None
		self.vecS = None
		# Derived values.
		self.matPsi = None
		self.vecPhi = None # "phi" = PSI^T * (-gamma)
		self.vecLambdaLS = None # "LS" = "launch, scaled"
		self.vecLambdaPSD = None # "PSD" = "positive-semi-definite" (+LS)
		self.vecLambdaWB = None # "WB" = "well-behaved" = strictly pos-def and fMin >= 0.0 (+LS)
		self.tCauchyLS = None
		self.tCauchyPSD = None
		self.tCauchyWB = None
	def dump(self):
		msg(f'Begin hessCurvesType().dump...')
		msg(f'self = {self}')
		msg(f'  vecXA = {self.vecXA}')
		msg(f'  matV.shape[0] = {self.matV.shape[0]}')
		msg(f'  matV.shape[1] = {self.matV.shape[1]}')
		msg(f'  matV = ...\n{self.matV}')
		msg(f'  vecS = {self.vecS}')
		msg(f'  matPsi = ...\n{self.matPsi}')
		msg(f'  vecPhi = {self.vecPhi}')
		msg(f'  vecLambdaLS = {self.vecLambdaLS}')
		msg(f'  vecLambdaPSD = {self.vecLambdaPSD}')
		msg(f'  vecLambdaWB = {self.vecLambdaWB}')
		msg(f'  tCauchyLS = {self.tCauchyLS}')
		msg(f'  tCauchyPSD = {self.tCauchyPSD}')
		msg(f'  tCauchyWB = {self.tCauchyWB}')
		msg(f'End hessCurvesType.dump().')
	# DRaburn 2023-03-13:
	#  Re-defining calcYLevLSOfMu() and calcYLevPSDOfMu()
	#  so that an input mu of 0.0 is always end of curve.
	#  Previously, end would be larger of { 0.0, -mineig }.
	def calcYLevLSOfMu(self, mu):
		EPS_WB0 = 1.0e-6
		EPS_WB1 = 1.0e-6
		vecLambda = self.vecLambdaLS
		mu1 = mu
		muThresh = EPS_WB0 + EPS_WB1*max(abs(vecLambda))
		if ( min(vecLambda) <= muThresh ):
			mu1 += abs(min(vecLambda)) + muThresh
		vecZeta = self.vecPhi / (vecLambda + mu1)
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYLevPSDOfMu(self, mu):
		EPS_WB0 = 1.0e-6
		EPS_WB1 = 1.0e-6
		vecLambda = self.vecLambdaPSD
		mu1 = mu
		muThresh = EPS_WB0 + EPS_WB1*max(abs(vecLambda))
		if ( min(vecLambda) <= muThresh ):
			mu1 += abs(min(vecLambda)) + muThresh
		vecZeta = self.vecPhi / (vecLambda + mu1)
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYLevWBOfMu(self, mu):
		assert( mu + min(self.vecLambdaWB) > 0.0 )
		vecZeta = self.vecPhi / (self.vecLambdaWB + mu)
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYFloorOfMu_noZeroMin(self, mu):
		EPS_WB0 = 1.0e-6
		EPS_WB1 = 1.0e-6
		vecLambda = self.vecLambdaLS
		mu1 = mu
		muThresh = EPS_WB0 + EPS_WB1*max(abs(vecLambda))
		if ( min(vecLambda) <= muThresh ):
			mu1 += abs(min(vecLambda)) + muThresh
		vecMu = vecLambda.copy()
		vecMu[vecMu<mu1] = mu1
		vecZeta = self.vecPhi / vecMu
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYFloorOfMu(self, mu):
		mu1 = mu + min(self.vecLambdaWB)
		vecMu = self.vecLambdaWB.copy()
		vecMu[vecMu<mu1] = mu1
		vecZeta = self.vecPhi / vecMu
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYCauchyLSOfT(self, t):
		return self.vecYLaunch + ((t*self.tCauchyLS) * self.vecS * (self.matPsi @ self.vecPhi))
	def calcYCauchyPSDOfT(self, t):
		return self.vecYLaunch + ((t*self.tCauchyPSD) * self.vecS * (self.matPsi @ self.vecPhi))
	def calcYCauchyWBOfT(self, t):
		return self.vecYLaunch + ((t*self.tCauchyWB) * self.vecS * (self.matPsi @ self.vecPhi))
	def calcXLevLSOfMu(self, mu):
		return self.vecXA + (self.matV @ self.calcYLevLSOfMu(mu))
	def calcXLevPSDOfMu(self, mu):
		return self.vecXA + (self.matV @ self.calcYLevPSDOfMu(mu))
	def calcXLevWBOfMu(self, mu):
		return self.vecXA + (self.matV @ self.calcYLevWBOfMu(mu))
# End class class hessCurvesType().

class calcHessModel_prm():
	def __init__(self):
		self.dropRelThresh = 1.0e-1
		self.dropAbsThresh = 1.0e-12
		self.fdOrder = 2
		#self.epsRelFD = 1.0e-4 # DR 02-26: too small post "rollfix".
		self.epsRelFD = 1.0e-2
	def dump(self):
		msg(f'Begin calcHessModel_prm().dump...')
		msg(f'self = {self}')
		msg(f'  dropRelThresh = {self.dropRelThresh}')
		msg(f'  dropAbsThresh = {self.dropAbsThresh}')
		msg(f'  fdOrder = {self.fdOrder}')
		msg(f'  epsRelFD = {self.epsRelFD}')
		msg(f'End calcHessModel_prm().dump.')
def calcHessModel( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, prm=calcHessModel_prm() ):
	sizeX = vecXAnchor.shape[0]
	matD = record_matX - np.reshape(vecXAnchor, (sizeX,1)) # Autobroadcast
	matV, vecKeep = danutil.utorthdrop(matD, prm.dropRelThresh, prm.dropAbsThresh)
	sizeK = matV.shape[1]
	if ( 0 == sizeK ):
		#msg('ERROR: sizeK is zero.')
		return None
	#
	vecGammaA = matV.T @ vecGAnchor
	matY = np.triu( matV.T @ matD[:,vecKeep] ) # Ideally could be calc'd by utorthdrop.
	matGamma = matV.T @ record_matG[:,vecKeep]
	#matA = np.linalg.solve(matY.T, (  matGamma - np.reshape(vecGammaAnchor, (sizeK,1))  ).T) # Does not require matM.
	#matH = (matA.T + matA)/2.0  # Does not require matM.
	matM = np.linalg.solve(matY.T, (  record_matG[:,vecKeep] - np.reshape(vecGAnchor, (sizeX,1)) ).T).T
	matH = danutil.symm(matV.T @ matM)
	matW = matM - (matV @ matH)
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV.copy()
	hessModel.fA = fAnchor.copy()
	hessModel.vecGammaA = vecGammaA.copy()
	hessModel.matH = matH.copy()
	hessModel.vecGPerpA = vecGAnchor.copy() - (matV @ vecGammaA)
	hessModel.matW = matW.copy()
	return hessModel
def oracle_calcHessModel( funch_evalFG, vecXAnchor, record_matX, prm=calcHessModel_prm() ):
	sizeX = vecXAnchor.shape[0]
	matD = record_matX - np.reshape(vecXAnchor, (sizeX,1)) # Autobroadcast
	matV, vecKeep = danutil.utorthdrop(matD, prm.dropRelThresh, prm.dropAbsThresh)
	sizeK = matV.shape[1]
	if ( 0 == sizeK ):
		#msg('ERROR: sizeK is zero.')
		return None
	#
	fAnchor, vecGAnchor = funch_evalFG(vecXAnchor)
	vecGammaA = matV.T @ vecGAnchor
	#
	dScale = max(np.sqrt(np.sum(matD**2,0)))
	assert ( dScale > 0.0 )
	epsFD = prm.epsRelFD * dScale
	#matA = np.zeros((sizeK, sizeK)) # Does not require matM.
	matM = np.zeros((sizeX, sizeK))
	if (1 == prm.fdOrder):
		for k in range(sizeK):
			msg(f'{k} / {sizeK}')
			fP, vecGP = funch_evalFG(vecXAnchor + (epsFD*matV[:,k]))
			#matA[:,k] = (matV.T @ (vecGP-vecGA))/epsFD # Does note require matM.
			matM[:,k] = (vecGP - vecGAnchor)/epsFD
	elif (2 == prm.fdOrder):
		for k in range(sizeK):
			msg(f'{k} / {sizeK}')
			fP, vecGP = funch_evalFG(vecXAnchor + (epsFD*matV[:,k]))
			fM, vecGM = funch_evalFG(vecXAnchor - (epsFD*matV[:,k]))
			#matA[:,k] = (matV.T @ (vecGP-vecGM))/(2.0*epsFD) # Does note require matM.
			matM[:,k] = (vecGP - vecGM)/(2.0*epsFD)
			msg(f'  fP = {fP}')
			msg(f'  fM = {fM}')
			msg(f'  vecGP = {vecGP}')
			msg(f'  vecGM = {vecGM}')
	#matH = (matA.T + matA)/2.0
	matH = danutil.symm(matV.T @ matM)
	matW = matM - (matV @ matH)
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV.copy()
	hessModel.fA = fAnchor
	hessModel.vecGammaA = vecGammaA.copy()
	hessModel.matH = matH.copy()
	hessModel.vecGPerpA = vecGAnchor.copy() - (matV @ vecGammaA)
	hessModel.matW = matW.copy()
	return hessModel
# End def calcHessModel() etc.

class calcHessCurves_prm():
	def __init__(self):
		self.epsRelWB = 1.0e-6
		self.adjustFForLaunch = False
	def dump(self):
		msg(f'Begin calcHessCurves_prm().dump...')
		msg(f'self = {self}')
		msg(f'  epsRelWB = {self.epsRelWB}')
		msg(f'  adjustFForLaunch = {self.adjustFForLaunch}')
		msg(f'End calcHessCurves_prm.dump().')
def calcHessCurves__calcLambdaWBFloor( f0, vecPhi, vecLambda, fFloor, epsRelWB ):
	assert (f0 > fFloor)
	def fRes( lambdaF ):
		fRes = f0 - fFloor
		for k in range( 0, vecLambda.shape[0] ):
			fRes -= ((vecPhi[k])**2) / (2.0*max( vecLambda[k], lambdaF ))
		return fRes
	if (min(vecLambda) > 0.0):
		if (fRes(0.0) >= 0.0):
			return 0.0
	lambdaHi = (vecPhi @ vecPhi) / (2.0*(f0 - fFloor))
	if (max(vecLambda) < lambdaHi):
		return lambdaHi
	lambdaLo = epsRelWB * max(vecLambda)
	if (fRes(lambdaLo) >= 0.0):
		return lambdaLo
	lambdaF = optimize.bisect( fRes, lambdaLo, lambdaHi )
	return lambdaF
def calcHessCurves( hessModel, vecYLaunch=None, vecS=None, prm=calcHessCurves_prm() ):
	# Parse input.
	sizeX = hessModel.matV.shape[0]
	sizeK = hessModel.matV.shape[1]
	if (vecYLaunch is None):
		vecYLaunch = np.zeros(sizeK)
	else:
		pass
		#msg('WARNING: Use of non-zero vecYLaunch has not been tested.')
	if (vecS is None):
		vecS = np.ones(sizeK)
	else:
		msg('WARNING: Use of scaling matrix has not been tested.')
	assert (min(vecS) > 0.0)
	#
	# Move to launch point, scale, and get eigen-factorization.
	if (prm.adjustFForLaunch):
		fLS = hessModel.fA + (vecYLaunch @ hessModel.vecGammaA) + ((vecYLaunch @ hessModel.matH @ vecYLaunch)/2.0)
	else:
		fLS = hessModel.fA
	vecGammaLS = (hessModel.vecGammaA + (hessModel.matH @ vecYLaunch)) * vecS
	matHLS = danutil.symm((hessModel.matH * vecS).T * vecS)
	vecLambdaCmplx, matPsi = np.linalg.eig(matHLS)
	vecPhi = matPsi.T @ (-vecGammaLS)
	for n in range(len(vecLambdaCmplx)):
		assert ( np.isreal(vecLambdaCmplx[n]) )
	vecLambdaLS = np.real(vecLambdaCmplx)
	#
	vecLambdaPSD = vecLambdaLS.copy()
	vecLambdaPSD[vecLambdaPSD < 0.0] = 0.0
	#
	if (fLS > 0.0):
		lambdaWBFloor = calcHessCurves__calcLambdaWBFloor( fLS, vecPhi, vecLambdaLS, 0.0, prm.epsRelWB )
		vecLambdaWB = vecLambdaLS.copy()
		vecLambdaWB[vecLambdaWB < lambdaWBFloor] = lambdaWBFloor
	else:
		msg('\n\n\n')
		msg('********************************************************************************')
		msg(f'WARNING: fLS <= 0.0 (fLS = {fLS} while fA = {hessModel.fA}).')
		msg(f'  Please consider a different launch point.')
		msg(f'  Setting every element of vecLambdaWB to a large negative value.')
		msg('********************************************************************************')
		msg('\n\n\n')
		vecLambdaWB = vecLambdaLS.copy()
		vecLambdaWB[:] = -1.0E10
	#
	hessCurves = hessCurvesType()
	hessCurves.vecXA = hessModel.vecXA
	hessCurves.matV = hessModel.matV.copy()
	hessCurves.vecYLaunch = vecYLaunch.copy()
	hessCurves.vecS = vecS.copy()
	hessCurves.matPsi = matPsi.copy()
	hessCurves.vecPhi = vecPhi.copy()
	hessCurves.vecLambdaLS = vecLambdaLS.copy()
	hessCurves.vecLambdaPSD = vecLambdaPSD.copy()
	hessCurves.vecLambdaWB = vecLambdaWB.copy()
	#
	# Calcualte modified models...
	# We could re-use the original anchor and shift the "LS" quantities back,
	# but, it's probably easier to just shift the anchor.
	# Note that we are NOT doing the anchor shift for the hessCurves.
	# DRaburn 2023-02-25:
	#  Actually, it'd be nice to undo this.
	#msg('NOTE: vecXA (anchor) is shifted for modified Hessian models.')
	# Note "hessModelLS" is merely for debug.
	hessModelLS = hessModel.copy()
	hessModelLS.vecXA = hessCurves.vecXA.copy()
	hessModelLS.fA = fLS.copy()
	hessModelLS.vecGammaA = vecGammaLS / vecS
	hessModelLS.vecGPerpA = hessModel.vecGPerpA + (hessModel.matW @ vecYLaunch)
	hessModelLS.matH = danutil.symm( (( matPsi @ np.diag(vecLambdaLS) @ matPsi.T )/vecS).T / vecS )
	assert ( danutil.reldiff(hessModelLS.matH, hessModel.matH) < 1.0e-4 )
	#
	hessModelPSD = hessModelLS.copy()
	hessModelPSD.matH = danutil.symm( (( matPsi @ np.diag(vecLambdaPSD) @ matPsi.T )/vecS).T / vecS )
	#
	hessModelWB = hessModelLS.copy()
	hessModelWB.matH = danutil.symm( (( matPsi @ np.diag(vecLambdaWB) @ matPsi.T )/vecS).T / vecS )
	#
	# DRaburn 2023-03-14:
	#  Do calculations for Cauchy.
	hessCurves.tCauchyLS  = danutil.linishrootOfQuad( (vecPhi.T @ (vecLambdaLS  * vecPhi))/2.0, -vecPhi.T @ vecPhi, hessModel.fA )
	hessCurves.tCauchyPSD = danutil.linishrootOfQuad( (vecPhi.T @ (vecLambdaPSD * vecPhi))/2.0, -vecPhi.T @ vecPhi, hessModel.fA )
	hessCurves.tCauchyWB  = danutil.linishrootOfQuad( (vecPhi.T @ (vecLambdaWB  * vecPhi))/2.0, -vecPhi.T @ vecPhi, hessModel.fA )
	#
	# Note that hessModelLS should behave identically to hessModel, just with a different anchor.
	return hessCurves, hessModelPSD, hessModelWB, hessModelLS
# End def calcHessCurves() etc.

class searchHessCurve_prm():
	def __init__(self):
		self.curveSelector = "floor"
		self.tMin = 0.0
		self.tMax = 1.0
		self.dTTol = 1.0e-6
		self.iterLimit = 1000
		self.returnTAndMu = False
	def dump(self):
		msg(f'Begin searchHessCurve_prm().dump...')
		msg(f'self = {self}')
		msg(f'  curveSelector = {self.curveSelector}')
		msg(f'  tMin = {self.tMin} (depricated?)')
		msg(f'  tMax = {self.tMax} (depricated?)')
		msg(f'  dTTol = {self.dTTol} (depricated?)')
		msg(f'  iterLimit = {self.iterLimit} (depricated?)')
		msg(f'  returnTAndMu = {self.returnTAndMu}')
		msg(f'End searchHessCurve_prm.dump().')
def searchHessCurve( funch_evalFG, hessCurves, prm=searchHessCurve_prm() ):
	if (prm.curveSelector.lower() == "floor"):
		funchYOfMu = hessCurves.calcYFloorOfMu
	elif (prm.curveSelector.lower() == "ls"):
		funchYOfMu = hessCurves.calcYLevLSOfMu
	elif (prm.curveSelector.lower() == "psd"):
		funchYOfMu = hessCurves.calcYLevPSDOfMu
	elif (prm.curveSelector.lower() == "wb"):
		funchYOfMu = hessCurves.calcYLevWBOfMu
	else:
		msg(f'ERROR: unsupported value of prm.curveSelector.lower() ({prm.curveSelector.lower()}).')
		return None
	muScl = min(hessCurves.vecLambdaWB)
	assert( muScl > 0.0 )
	def xOfT(t):
		if (t == 0.0):
			return hessCurves.vecXA + (hessCurves.matV @ hessCurves.vecYLaunch )
		else:
			mu = muScl*((1.0/t) - 1.0)
			return hessCurves.vecXA + (hessCurves.matV @ funchYOfMu(mu))
	def fOfT( t ):
		f, _ = funch_evalFG(xOfT(t))
		# Optimization: make use of vecG. (We have to mux with the curve to get df/dt, though.)
		#scipy.optimize.fmin(func, x0, args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0, callback=None, initial_simplex=None)
		return f
	
	# DRaburn 2023-03-20.
	#  Let's end this.
	#zzz
	t1, f1 = findMin( fOfT, 0.0, 1.0 )
	msg(f'FOUND: {t1:15.9e}, {f1:15.9e}')
	#if ( t1 < 1.0e-7 or 1.0-t1 < 1.0e-7 ):
	if (False):
		msg('Generating plot...')
		numPts = 20
		tVals = np.linspace(0.00, 1.00, numPts)
		fVals = np.zeros(numPts)
		dVals = np.zeros(numPts)
		for n in range (numPts):
			fVals[n] = fOfT(tVals[n])
			dVals[n] = norm(xOfT(tVals[n])-xOfT(0.0))
			msg(f'  {n:3d}/{numPts:d} {tVals[n]:10.3e}, {dVals[n]:10.3e}, {fVals[n]:10.3e}')
		import matplotlib.pyplot as plt
		plt.semilogy(tVals, fVals, 'o-')
		plt.semilogy(t1, f1, 'p', markersize=20)
		plt.grid(True)
		plt.show()
	
	if ( t1 > 0.0 ):
		mu1 = muScl*((1.0/t1)-1.0)
	else:
		mu1 = 0.0
	if (prm.returnTAndMu):
		return xOfT(t1), t1, mu1
	return xOfT(t1)
	
	
	# DRaburn 2023-03-16.
	#  fminbound() does not do what I need. Enough of this.
	f0 = fOfT(0.0)
	t1 = 1.0
	#for n in range(30): # Reduce this.
	#for n in range(10):
	for n in range(20):
		f1 = fOfT(t1)
		if ( f1 <= f0 ):
			if (prm.returnTAndMu):
				msg(f'FOUND(ISH): {t1:9.3e}, {f1:9.3e}')
				return xOfT(t1), t1, muScl*((1.0/t1)-1.0)
			return xOfT(t1)
		t1 /= 2.0
		#t1 /= 5.0
	#msg('ERROR: Failed to reduce f')
	if (prm.returnTAndMu):
		return xOfT(0.0), 0.0, -1.0
	return xOfT(0.0)
	
	#msg('Calling fminbound...')
	# DRaburn 2023-03-16.
	#  Note that this search is inefficient, not making use of gradient information.
	#  But, this is merely for testing.
	#  Also, because I'm not returning f, redundant calculations are needed later.
	# DRaburn 2023-03-16.
	#  optimize.fminbound() can't handle NaN and other exceptions? Really?
	t0 = prm.tMin
	t1 = prm.tMax
	f0 = fOfT(t0)
	f1 = fOfT(t1)
	for n in range(20):
		if ( f1 <= 10.0*f0 ):
			break
		t1 = t0 + (t1-t0)/2.0
		f1 = fOfT(t1)
	assert ( f1 <= 100.0*f0 )
	t = optimize.fminbound( fOfT, t0, t1, xtol=prm.dTTol, maxfun=prm.iterLimit+2, disp=1 )
	f = fOfT(t)
	msg(f'{t0}, {t}, {t1}; {f0}, {f}, {f1}')
	return xOfT(t)
def multiSearchHessCurve( funch_evalFG, hessCurves ):
	shcPrm = searchHessCurve_prm()
	shcPrm.curveSelector = "floor"
	vecXFloor = searchHessCurve( funch_evalFG, hessCurves, shcPrm )
	fFloor = funch_evalFG(vecXFloor)
	shcPrm.curveSelector = "ls"
	vecXLS = searchHessCurve( funch_evalFG, hessCurves, shcPrm )
	fLS = funch_evalFG(vecXLS)
	#shcPrm.curveSelector = "psd"
	#vecXPSD = searchHessCurve( funch_evalFG, hessCurves, shcPrm )
	#fPSD = funch_evalFG(vecXPSD)
	#shcPrm.curveSelector = "wb"
	#vecXWB = searchHessCurve( funch_evalFG, hessCurves, shcPrm )
	#fWB = funch_evalFG(vecXWB)
	#
	vecX = vecXFloor
	f = fFloor
	if ( fLS < f ):
		vecX = vecXLS
		f = fLS
	#if ( fPSD < f ):
	#	vecX = vecXPSD
	#	f = fPSD
	#if ( fWB < f ):
	#	vecX = vecXWB
	#	f = fWB
	return vecX
# End searchHessCurve() etc.

class multiOracleMin_prm():
	def __init__(self):
		dummy = None
	def dump(self):
		msg(f'Begin multiOracleMin_prm.dump...')
		msg(f'self = {self}')
		msg(f'  dummy = {self.dummy}')
		msg(f'End multiOracleMin_prm.dump().')
def multiOracleMin(
  funch_evalFG,
  vecXAnchor,
  fAnchor,
  vecGAnchor,
  record_matX,
  record_vecF,
  record_matG,
  vecXLaunch = None,
  prm = multiOracleMin_prm() ):
	if ( vecXLaunch is None ):
		vecXLaunch = vecXAnchor.copy()
	mortal_chmPrm = calcHessModel_prm()
	oracle_chmPrm = calcHessModel_prm()
	oracle_chmPrm.dropRelThresh /= 10.0
	mortal_hm = calcHessModel( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, mortal_chmPrm )
	oracle_hm = oracle_calcHessModel( funch_evalFG, vecXAnchor, record_matX, oracle_chmPrm )
	mortal_vecYLaunch = mortal_hm.matV.T @ ( vecXLaunch - mortal_hm.vecXA )
	oracle_vecYLaunch = oracle_hm.matV.T @ ( vecXLaunch - oracle_hm.vecXA )
	vecS = None
	mortal_chcPrm = calcHessCurves_prm()
	oracle_chcPrm = calcHessCurves_prm()
	mortal_hc, _, _, _ = calcHessCurves( mortal_hm, mortal_vecYLaunch, vecS, mortal_chcPrm )
	oracle_hc, _, _, _ = calcHessCurves( oracle_hm, oracle_vecYLaunch, vecS, oracle_chcPrm )
	mortal_vecX = multiSearchHessCurve( funch_evalFG, mortal_hc )
	oracle_vecX = multiSearchHessCurve( funch_evalFG, oracle_hc )
	mortal_f = funch_evalFG(mortal_vecX)
	oracle_f = funch_evalFG(oracle_vecX)
	if ( oracle_f < mortal_f ):
		return oracle_vecX
	else:
		return mortal_vecX
# End multiOracleMin().

class multiSearchMin_prm():
	def __init__(self):
		dummy = None
	def dump(self):
		msg(f'Begin multiSearchMin_prm.dump...')
		msg(f'self = {self}')
		msg(f'  dummy = {self.dummy}')
		msg(f'End multiSearchMin_prm.dump().')
def multiSearchMin(
  funch_evalFG,
  vecXAnchor,
  fAnchor,
  vecGAnchor,
  record_matX,
  record_vecF,
  record_matG,
  vecXLaunch = None,
  prm = multiSearchMin_prm() ):
	if ( vecXLaunch is None ):
		vecXLaunch = vecXAnchor.copy()
	chmPrm = calcHessModel_prm()
	hm = calcHessModel( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, chmPrm )
	vecYLaunch = hm.matV.T @ ( vecXLaunch - hm.vecXA )
	vecS = None
	chcPrm = calcHessCurves_prm()
	hc, _, _, _ = calcHessCurves( hm, vecYLaunch, vecS, chcPrm )
	vecX = multiSearchHessCurve( funch_evalFG, hc )
	return vecX
# End multiSearchMin().

def searchMin(
  funch_evalFG,
  vecXAnchor,
  fAnchor,
  vecGAnchor,
  record_matX,
  record_vecF,
  record_matG,
  vecXLaunch = None,
  chmPrm = calcHessModel_prm(),
  chcPrm = calcHessCurves_prm(),
  shcPrm = searchHessCurve_prm() ):
	if ( vecXLaunch is None ):
		vecXLaunch = vecXAnchor.copy()
	hm = calcHessModel( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, chmPrm )
	vecYLaunch = hm.matV.T @ ( vecXLaunch - hm.vecXA )
	vecS = None
	hc, _, _, _ = calcHessCurves( hm, vecYLaunch, vecS, chcPrm )
	vecX = searchHessCurve( funch_evalFG, hc, shcPrm )
	return vecX
# End searchMin().

# _oracleP: oracle wrt calculating vecPLand.
# This is a bit beyond "hessmodel", in "qnj" territory.
class smop_dat():
	def __init(self):
		pass
def searchMin_sgd_oracleP(
  funch_evalFG,
  vecXAnchor,
  fAnchor,
  vecGAnchor,
  record_matX,
  record_vecF,
  record_matG,
  vecXLaunch,
  vecPLaunch,
  chmPrm = calcHessModel_prm(),
  chcPrm = calcHessCurves_prm(),
  shcPrm = searchHessCurve_prm() ):
	smopDat = smop_dat()
	smopDat.t = 0.0
	smopDat.mu = -1.0
	numRecords = record_matX.shape[1]
	if ( 0 == numRecords ):
		return vecXLaunch.copy(), vecPLaunch.copy()
	hm = calcHessModel( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, chmPrm )
	smopDat.hm = hm
	if ( hm is None ):
		sizeK = 0
	else:
		sizeK = hm.matV.shape[1]
	if ( 0 == sizeK ):
		return vecXLaunch.copy(), vecPLaunch.copy(), smopDat
	vecYLaunch = hm.matV.T @ ( vecXLaunch - hm.vecXA )
	smopDat.vecYLaunch = vecYLaunch
	vecS = None
	hc, hmPSD, hmWB, hmLS = calcHessCurves( hm, vecYLaunch, vecS, chcPrm )
	smopDat.hc = hc
	smopDat.hmPSD = hmPSD
	smopDat.hmWB = hmWB
	smopDat.hmLS = hmLS
	# Note: hmLS should functionally match hm.
	shcPrm.returnTAndMu = True
	#vecXLand = searchHessCurve( funch_evalFG, hc, shcPrm )
	###vecXLand, t, mu = searchHessCurve( funch_evalFG, hc, shcPrm )
	vecXLand, t, mu = searchHessCurve0321( funch_evalFG, hm, hc)
	
	# We'll use oracle/extra info for vecPLand...
	fLaunch, vecGLaunch = funch_evalFG( vecXLaunch )
	assert (vecGLaunch @ vecGLaunch > 0.0 )
	a = ( vecPLaunch @ vecGLaunch ) / (vecGLaunch @ vecGLaunch)
	vecB = vecPLaunch - (a*vecGLaunch)
	assert( norm(vecGLaunch @ vecB) <= 1.0e-6*(norm(vecGLaunch) * norm(vecB)) )
	fLand, vecGLand = funch_evalFG( vecXLand )
	alphaF = fLand/fLaunch
	vecPLand = (a*vecGLand) + (alphaF*vecB)
	#msg(f'{norm(vecXLand - vecXLaunch)}, {fLaunch}, {fLand}')
	return vecXLand, vecPLand, smopDat
def searchMin_sgd(
  funch_evalFG,
  vecXAnchor,
  fAnchor,
  vecGAnchor,
  record_matX,
  record_vecF,
  record_matG,
  vecXLaunch,
  vecPLaunch,
  chmPrm = calcHessModel_prm(),
  chcPrm = calcHessCurves_prm(),
  shcPrm = searchHessCurve_prm() ):
	smopDat = smop_dat()
	smopDat.t = 0.0
	smopDat.mu = -1.0
	numRecords = record_matX.shape[1]
	if ( 0 == numRecords ):
		return vecXLaunch.copy(), vecPLaunch.copy()
	hm = calcHessModel( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, chmPrm )
	smopDat.hm = hm
	if ( hm is None ):
		sizeK = 0
	else:
		sizeK = hm.matV.shape[1]
	if ( 0 == sizeK ):
		return vecXLaunch.copy(), vecPLaunch.copy(), smopDat
	vecYLaunch = hm.matV.T @ ( vecXLaunch - hm.vecXA )
	smopDat.vecYLaunch = vecYLaunch
	vecS = None
	hc, hmPSD, hmWB, hmLS = calcHessCurves( hm, vecYLaunch, vecS, chcPrm )
	smopDat.hc = hc
	smopDat.hmPSD = hmPSD
	smopDat.hmWB = hmWB
	smopDat.hmLS = hmLS
	# Note: hmLS should functionally match hm.
	shcPrm.returnTAndMu = True
	#vecXLand = searchHessCurve( funch_evalFG, hc, shcPrm )
	###vecXLand, t, mu = searchHessCurve( funch_evalFG, hc, shcPrm )
	vecXLand, t, mu = searchHessCurve0321( funch_evalFG, hm, hc)
	smopDat.t = t
	smopDat.mu = mu
	
	funch_evalFG = None # Don't oracle past here!
	# DRaburn 2023-03-17: Second attempt, only keep perp part.
	fLaunch, vecGLaunch = hm.evalFGOfProjX( vecXLaunch )
	fLand, vecGLand = hm.evalFGOfProjX( vecXLand )
	vecPLand = (fLand/fLaunch)*( vecPLaunch - (hm.matV @ (hm.matV.T @ vecPLaunch)) )
	return vecXLand, vecPLand, smopDat
	
	# DRaburn 2023-03-17: First attempt, below doesn't work...
	# We need to get a suitable vecPLand...
	fLaunch, vecGLaunch = hm.evalFGOfProjX( vecXLaunch )
	assert (vecGLaunch @ vecGLaunch > 0.0 )
	a = ( vecPLaunch @ vecGLaunch ) / (vecGLaunch @ vecGLaunch)
	vecB = vecPLaunch - (a*vecGLaunch)
	assert( norm(vecGLaunch @ vecB) <= 1.0e-6*(norm(vecGLaunch) * norm(vecB)) )
	fLand, vecGLand = hm.evalFGOfProjX( vecXLand )
	alphaF = fLand/fLaunch
	vecPLand = (a*vecGLand) + (alphaF*vecB)
	#msg(f'{norm(vecXLand - vecXLaunch)}, {fLaunch}, {fLand}')
	return vecXLand, vecPLand, smopDat

class findMin_prm():
	def __init__(self):
		self.xTol = 1.0e-5
def findMin( funch_evalFOfT, xLo, xHi, prm=findMin_prm() ):
	xL = xLo
	xR = xHi
	fL = funch_evalFOfT(xL)
	fR = funch_evalFOfT(xR)
	if ( fR < fL ):
		xC = (xL+xR)/2.0
		fC = funch_evalFOfT(xC)
	else:
		for n in range(10):
			xC = xL + (xR-xL)/10.0
			fC = funch_evalFOfT(xC)
			#msg(f'{n:2d}: {xL:9.3e}, {xC:9.3e}, {xR:9.3e};  {xL-xC:9.3e}, {xR-xC:9.3e};  {fL:9.3e}, {fC:9.3e}, {fR:9.3e};  {fL-fC:10.3e}, {fR-fC:10.3e}')
			msg(f'{n:2d}: {xL:15.9e}, {xC:15.9e}, {xR:15.9e};  {fL:15.9e}, {fC:15.9e}, {fR:15.9e}')
			if ( fC < fL ):
				break
			xR = xC
			fR = fC
		if ( fC >= fL ):
			#msg(f'{xL:9.3e}, {xC:9.3e}, {xR:9.3e};  {xL-xC:9.3e}, {xR-xC:9.3e};  {fL:9.3e}, {fC:9.3e}, {fR:9.3e};  {fL-fC:10.3e}, {fR-fC:10.3e}')
			msg('EXCEPTION: Failed to reduce f with small step.')
			return xL, fL
	for n in range(50):
		#msg(f'{n:2d}: {xL:9.3e}, {xC:9.3e}, {xR:9.3e};  {xL-xC:9.3e}, {xR-xC:9.3e};  {fL:9.3e}, {fC:9.3e}, {fR:9.3e};  {fL-fC:10.3e}, {fR-fC:10.3e}')
		msg(f'{n:2d}: {xL:15.9e}, {xC:15.9e}, {xR:15.9e};  {fL:15.9e}, {fC:15.9e}, {fR:15.9e}')
		assert( xL < xC )
		assert( xC < xR )
		if ( fC >= fL ):
			msg('EXCEPTION: Bad shape.')
			return xL, fL
		if ( abs(xR - xL) < prm.xTol ):
			if ( fR < fC ):
				return xR, fR
			return xC, fC
		if ( abs(xR - xC) > abs(xC - xL) ):
			#xTemp = (xR+xC)/2.0
			xTemp = xC + (xR-xC)/3.0
			fTemp = funch_evalFOfT(xTemp)
			#msg(f'    {xTemp:9.3e};  {xTemp-xC:9.3e};  {fTemp:9.3e};  {fTemp-fC:9.3e}')
			msg(f'    {xTemp:15.9e}; {fTemp:15.9e}')
			if ( fTemp < fC ):
				# C -> L, Temp -> C
				xL = xC
				fL = fC
				xC = xTemp
				fC = fTemp
			else:
				# Temp -> R
				xR = xTemp
				fR = fTemp
		else:
			#xTemp = (xL+xC)/2.0
			xTemp = xC + (xL-xC)/3.0
			fTemp = funch_evalFOfT(xTemp)
			#msg(f'    {xTemp:9.3e};  {xTemp-xC:9.3e};  {fTemp:9.3e};  {fTemp-fC:9.3e}')
			msg(f'    {xTemp:15.9e}; {fTemp:15.9e}')
			if ( fTemp < fC ):
				# C -> R, Temp ->C
				xR = xC
				fR = fC
				xC = xTemp
				fC = fTemp
			else:
				# Temp -> L
				xL = xTemp
				fL = fTemp
		xTemp = None
		fTemp = None
	msg('EXCEPTION: Failed to converge within allowed iterations.')
	return xL, fL
# End findMin().

class findZero_prm():
	def __init__(self):
		self.xTol = 1.0e-6
def findZero( funch_evalFOfT, xLo, xHi, prm=findZero_prm() ):
	msg(f'*** YAY, WE ARE CALLING findZero()! prm.xTol = {prm.xTol} ***')
	xL = xLo
	xR = xHi
	fL = funch_evalFOfT(xL)
	fR = funch_evalFOfT(xR)
	msg(f'{xL:15.9e}, {xR:15.9e};  {fL:15.9e}, {fR:15.9e}')
	if ( fL*fR > 0.0 ):
		msg('EXCEPTION: Initial guesses do not (properly) bracket a zero.')
		return xL
	if ( 0.0 == fL ):
		msg('EXCEPTION: Initial fL == 0.0.')
		return xL
	if ( 0.0 == fR ):
		msg('EXCEPTION: Initial fR == 0.0.')
		return xR
	for n in range (30):
		msg(f'{n:3d};  {xL:15.9e}, {xR:15.9e};  {fL:15.9e}, {fR:15.9e}')
		# We'll use basic secant-bisection.
		assert( fR > fL )
		xTemp = ((xL*fR) - (xR*fL))/(fR - fL)
		assert( xL <= xTemp )
		assert( xTemp <= xR )
		fTemp = funch_evalFOfT(xTemp)
		msg(f'      {xTemp:15.9e}; {fTemp:15.9e}')
		if ( 0.0 == fTemp ):
			return xTemp
		if ( fL * fTemp <= 0.0 ):
			xR = xTemp
			fR = fTemp
		else:
			xL = xTemp
			fL = fTemp
		if ( abs(xR-xL) < prm.xTol ):
			if ( abs(fR) < abs(fL) ):
				return xR
			else:
				return xL
		# Now, bisect.
		msg(f'{n:3d};  {xL:15.9e}, {xR:15.9e};  {fL:15.9e}, {fR:15.9e}')
		xTemp = (xL+xR)/2.0
		fTemp = funch_evalFOfT(xTemp)
		msg(f'      {xTemp:15.9e}; {fTemp:15.9e}')
		if ( 0.0 == fTemp ):
			return xTemp
		if ( fL * fTemp <= 0.0 ):
			xR = xTemp
			fR = fTemp
		else:
			xL = xTemp
			fL = fTemp
	msg('EXCEPTION: Failed to converge within allowed iterations.')
	return xL
# End findZero()

class searchHessCurve0321_prm():
	def __init__(self):
		self.curveSelector = "floor"
		#self.returnTAndMu = False
		#self.agreementCoeff = 0.5
		self.agreementCoeff = -1.0
	def dump(self):
		msg(f'Begin searchHessCurve0321_prm().dump...')
		msg(f'self = {self}')
		msg(f'  curveSelector = {self.curveSelector}')
		#msg(f'  returnTAndMu = {self.returnTAndMu}')
		msg(f'  agreementCoeff = {self.agreementCoeff}')
		msg(f'End searchHessCurve0321_prm.dump().')
def searchHessCurve0321( funch_evalFG, hessModel, hessCurves, prm=searchHessCurve0321_prm() ):
	if (prm.curveSelector.lower() == "floor"):
		funchYOfMu = hessCurves.calcYFloorOfMu
	elif (prm.curveSelector.lower() == "ls"):
		funchYOfMu = hessCurves.calcYLevLSOfMu
	elif (prm.curveSelector.lower() == "psd"):
		funchYOfMu = hessCurves.calcYLevPSDOfMu
	elif (prm.curveSelector.lower() == "wb"):
		funchYOfMu = hessCurves.calcYLevWBOfMu
	else:
		msg(f'ERROR: unsupported value of prm.curveSelector.lower() ({prm.curveSelector.lower()}).')
		return None
	muScl = min(hessCurves.vecLambdaWB)
	assert( muScl > 0.0 )
	def xOfT(t):
		if (t == 0.0):
			return hessCurves.vecXA + (hessCurves.matV @ hessCurves.vecYLaunch )
		else:
			mu = muScl*((1.0/t) - 1.0)
			return hessCurves.vecXA + (hessCurves.matV @ funchYOfMu(mu))
	def fOfT( t ):
		f, _ = funch_evalFG(xOfT(t))
		# Optimization: make use of vecG. (We have to mux with the curve to get df/dt, though.)
		#scipy.optimize.fmin(func, x0, args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0, callback=None, initial_simplex=None)
		return f
	
	t1, f1 = findMin( fOfT, 0.0, 1.0 )
	msg(f'FOUND: {t1:15.9e}, {f1:15.9e}')
	
	t1_min = t1
	f1_min = f1
	if (prm.agreementCoeff > 0.0):
		def zOfT( t ):
			vecX0 = xOfT(0.0)
			f0, _ = hessModel.evalFGOfProjX(vecX0)
			vecX = xOfT(t)
			fActual, _ = funch_evalFG(vecX)
			fModel, _ = hessModel.evalFGOfProjX(vecX)
			return fActual - fModel - (prm.agreementCoeff*f0)
		z1 = zOfT(t1)
		msg(f'   z1 = {z1:15.9e}')
		if (z1 > 0.0):
			t1 = findZero( zOfT, 0.0, t1 )
			z1 = zOfT(t1)
			msg(f'Found   t1 = {t1:15.9e}, z1 = {z1:15.9e}')
		t1_zer = t1
		#
	if (False):
		msg('Generating plot...')
		numPts = 20
		tVals = np.linspace(0.00, 1.00, numPts)
		fVals = np.zeros(numPts)
		mVals = np.zeros(numPts)
		dVals = np.zeros(numPts)
		for n in range (numPts):
			fVals[n] = fOfT(tVals[n])
			mVals[n], _ = hessModel.evalFGOfProjX(xOfT(tVals[n]))
			dVals[n] = norm(xOfT(tVals[n])-xOfT(0.0))
			msg(f'  {n:3d}/{numPts:d} {tVals[n]:10.3e}, {dVals[n]:10.3e}, {fVals[n]:10.3e}, {mVals[n]:10.3e}')
		f1_zer = fOfT(t1_zer)
		import matplotlib.pyplot as plt
		plt.semilogy(tVals, fVals, 'o-')
		plt.semilogy(tVals, mVals, 'x-')
		plt.semilogy(t1_min, f1_min, 'p', markersize=20)
		plt.semilogy(t1_zer, f1_zer, '*', markersize=15)
		plt.grid(True)
		plt.show()
		
	
	if ( t1 > 0.0 ):
		mu1 = muScl*((1.0/t1)-1.0)
	else:
		mu1 = 0.0
	#if (prm.returnTAndMu):
	#	return xOfT(t1), t1, mu1
	#return xOfT(t1)
	return xOfT(t1), t1, mu1
# End searchHessCurve0321().
