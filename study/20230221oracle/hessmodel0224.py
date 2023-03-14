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
		vecY = self.getYOfProjX(vecX)
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
		if ( mu1 + min(vecLambda) <= 0.0 ):
			mu1 = abs(min(vecLambda)) + EPS_WB0 + (EPS_WB1 * max(abs(vecLambda)))
		vecZeta = self.vecPhi / (vecLambda + mu1)
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYLevPSDOfMu(self, mu):
		EPS_WB0 = 1.0e-6
		EPS_WB1 = 1.0e-6
		vecLambda = self.vecLambdaPSD
		mu1 = mu
		if ( mu1 + min(vecLambda) <= 0.0 ):
			mu1 = abs(min(vecLambda)) + EPS_WB0 + (EPS_WB1 * max(abs(vecLambda)))
		vecZeta = self.vecPhi / (vecLambda + mu1)
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYLevWBOfMu(self, mu):
		assert( mu + min(self.vecLambdaWB) > 0.0 )
		vecZeta = self.vecPhi / (self.vecLambdaWB + mu)
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
	def calcYFloorOfMu(self, mu):
		EPS_WB0 = 1.0e-6
		EPS_WB1 = 1.0e-6
		vecLambda = self.vecLambdaLS
		mu1 = mu
		if ( max( mu1, min(vecLambda) ) <= 0.0 ):
			mu1 = EPS_WB0 + (EPS_WB1 * max(abs(vecLambda)))
		vecMu = vecLambda.copy()
		vecMu[vecMu<mu1] = mu1
		vecZeta = self.vecPhi / vecMu
		vecY = self.vecYLaunch + self.vecS * (self.matPsi @ vecZeta)
		return vecY
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
		msg('ERROR: sizeK is zero.')
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
		msg('ERROR: sizeK is zero.')
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
	def dump(self):
		msg(f'Begin calcHessCurves_prm().dump...')
		msg(f'self = {self}')
		msg(f'  epsRelWB = {self.epsRelWB}')
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
		msg('WARNING: Use of non-zero vecYLaunch has not been tested.')
	if (vecS is None):
		vecS = np.ones(sizeK)
	else:
		msg('WARNING: Use of scaling matrix has not been tested.')
	assert (min(vecS) > 0.0)
	#
	# Move to launch point, scale, and get eigen-factorization.
	fLS = hessModel.fA + (vecYLaunch @ hessModel.vecGammaA) + ((vecYLaunch @ hessModel.matH @ vecYLaunch)/2.0)
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
		msg(f'WARNING: fLS <= 0.0 (fLS = {fLS} while fA = {hessModel.fA}).')
		msg(f'  Please consider a different launch point.')
		msg(f'  Setting every element of vecLambdaWB to a large negative value.')
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
	msg('NOTE: vecXA (anchor) is shifted for modified Hessian models.')
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
	def dump(self):
		msg(f'Begin searchHessCurve_prm().dump...')
		msg(f'self = {self}')
		msg(f'  curveSelector = {self.curveSelector}')
		msg(f'  dTTol = {self.dTTol}')
		msg(f'  iterLimit = {self.iterLimit}')
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
			return hessCurves.vecXA
		else:
			mu = muScl*((1.0/t) - 1.0)
			return hessCurves.vecXA + (hessCurves.matV @ funchYOfMu(mu))
	def fOfT( t ):
		f, _ = funch_evalFG(xOfT(t))
		# Optimization: make use of vecG. (We have to mux with the curve to get df/dt, though.)
		#scipy.optimize.fmin(func, x0, args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0, callback=None, initial_simplex=None)
		return f
	msg('Calling fminbound...')
	t = optimize.fminbound( fOfT, prm.tMin, prm.tMax, xtol=prm.dTTol, maxfun=prm.iterLimit+2, disp=1 )
	msg(f'  t = {t}')
	return xOfT(t)
