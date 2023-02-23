import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm
from scipy import optimize

class hessModelType():
	def __init(self):
		self.vecXA = None
		self.matV = None
		self.fA = None
		self.vecGammaA = None
		self.matH = None
		self.havePerp = False
		self.vecGPerpA = None
		self.matW = None
	def dump(self):
		msg(f'Begin hessModelType().dump...')
		msg(f'self = {self}')
		msg(f'  vecXA = {self.vecXA}')
		msg(f'  matV = ...\n{self.matV}')
		msg(f'  fA = {self.fA}')
		msg(f'  vecGammaA = {self.vecGammaA}')
		msg(f'  matH = ...\n{self.matH}')
		msg(f'  vecGPerpA = {self.vecGPerpA}')
		msg(f'  matW = ...\n{self.matW}')
		msg(f'End hessModelType.dump().')
	def evalFGGamma(self, vecX):
		# This eval() is mostly for illustrative purpose.
		vecY = self.matV.T @ ( vecX - self.vecXA )
		vecHY = self.matH @ vecY
		vecGamma = self.vecGammaA + vecHY
		f = self.fA + (vecY @ self.vecGammaA) + ((vecY @ vecHY)/2.0)
		if (self.havePerp):
			vecG = (self.matV @ vecGamma) + self.vecGPerpA + (self.matW @ vecY)
		else:
			vecG = None
		return f, vecG, vecGamma
	def evalFGammaGPerpOfY(self, vecY):
		# This eval() is also mostly for illustrative purpose.
		vecHY = self.matH @ vecY
		vecGamma = self.vecGammaA + vecHY
		f = self.fA + (vecY @ self.vecGammaA) + ((vecY @ vecHY)/2.0)
		if (self.havePerp):
			vecGPerp = self.vecGPerpA + (self.matW @ vecY)
		else:
			vecGPerp = None
		return f, vecGamma, vecGPerp
# End class hessModelType().

class calcHessModel_prm():
	def __init__(self):
		self.dropRelThresh = 1.0e-1
		self.dropAbsThresh = 1.0e-12
		self.epsFD = 1.0e-4
		self.fdOrder = 2
	def dump(self):
		msg(f'Begin calcHessModel_prm().dump...')
		msg(f'self = {self}')
		msg(f'  dropRelThresh = {self.dropRelThresh}')
		msg(f'  dropAbsThresh = {self.dropAbsThresh}')
		msg(f'  epsFD = {self.epsFD}')
		msg(f'  fdOrder = {self.fdOrder}')
		msg(f'End calcHessModel_prm().dump.')
def calcHessModel_basic( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, prm=calcHessModel_prm() ):
	# Note: record_vecF is not actually used, but is included for consistency.
	sizeX = vecXAnchor.shape[0]
	matD = record_matX - np.reshape(vecXAnchor, (sizeX,1)) # Autobroadcast
	matV, vecKeep = danutil.utorthdrop(matD, prm.dropRelThresh, prm.dropAbsThresh)
	sizeK = matV.shape[1]
	#
	vecGammaA = matV.T @ vecGAnchor
	matY = np.triu( matV.T @ matD[:,vecKeep] ) # Ideally could be calc'd by utorthdrop.
	matGamma = matV.T @ record_matG[:,vecKeep]
	#matA = np.linalg.solve(matY.T, (  matGamma - np.reshape(vecGammaAnchor, (sizeK,1))  ).T) # Does not require matM.
	matM = np.linalg.solve(matY.T, (  record_matG[:,vecKeep] - np.reshape(vecGAnchor, (sizeX,1)) ).T).T
	matA = matV.T @ matM
	matH = ( matA.T + matA )/2.0
	matW = matM - ( matV @ matH )
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV
	hessModel.fA = fAnchor.copy()
	hessModel.vecGammaA = vecGammaA
	hessModel.matH = matH
	hessModel.havePerp = True
	hessModel.vecGPerpA = vecGAnchor - (matV @ vecGammaA)
	hessModel.matW = matW
	return hessModel
def calcHessModel_basicOracle( vecXAnchor, record_matX, funch_evalFG, prm=calcHessModel_prm() ):
	sizeX = vecXAnchor.shape[0]
	matD = record_matX - np.reshape(vecXAnchor, (sizeX,1)) # Autobroadcast
	matV, vecKeep = danutil.utorthdrop(matD, prm.dropRelThresh, prm.dropAbsThresh)
	sizeK = matV.shape[1]
	#
	fA, vecGA = funch_evalFG(vecXAnchor)
	vecGammaA = matV.T @ vecGA
	matM = np.zeros((sizeX, sizeK))
	#matA = np.zeros((sizeK, sizeK)) # Does not require matM.
	if (1 == prm.fdOrder):
		for k in range(sizeK):
			fP, vecGP = funch_evalFG(vecXAnchor + (prm.epsFD*matV[:,k]))
			#matA[:,k] = (matV.T @ (vecGP-vecGA))/prm.epsFD # Does note require matM.
			matM[:,k] = (vecGP - vecGA)/prm.epsFD
	elif (2 == prm.fdOrder):
		for k in range(sizeK):
			fP, vecGP = funch_evalFG(vecXAnchor + (prm.epsFD*matV[:,k]))
			fM, vecGM = funch_evalFG(vecXAnchor - (prm.epsFD*matV[:,k]))
			#matA[:,k] = (matV.T @ (vecGP-vecGM))/(2.0*prm.epsFD) # Does note require matM.
			matM[:,k] = (vecGP - vecGM)/(2.0*prm.epsFD)
	matA = matV.T @ matM
	matH = ( matA.T + matA )/2.0
	matW = matM - ( matV @ matH )
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV
	hessModel.fA = fA
	hessModel.vecGammaA = vecGammaA
	hessModel.matH = matH
	hessModel.havePerp = True
	hessModel.vecGPerpA = vecGA - (matV @ vecGammaA)
	hessModel.matW = matW
	return hessModel
def calcHessModel_fullspace( vecXAnchor, funch_evalFG, prm=calcHessModel_prm() ):
	sizeX = vecXAnchor.shape[0]
	matV = np.eye(sizeX)
	sizeK = sizeX
	#
	fA, vecGA = funch_evalFG(vecXAnchor)
	matA = np.zeros((sizeX, sizeX))
	if (1 == prm.fdOrder):
		for k in range(sizeK):
			fP, vecGP = funch_evalFG(vecXAnchor + (prm.epsFD*matV[:,k]))
			matA[:,k] = (vecGP-vecGA)/prm.epsFD
	elif (2 == prm.fdOrder):
		for k in range(sizeK):
			fP, vecGP = funch_evalFG(vecXAnchor + (prm.epsFD*matV[:,k]))
			fM, vecGM = funch_evalFG(vecXAnchor - (prm.epsFD*matV[:,k]))
			matA[:,k] = (vecGP-vecGM)/(2.0*prm.epsFD)
	matH = ( matA.T + matA )/2.0
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV
	hessModel.fA = fA
	hessModel.vecGammaA = vecGA
	hessModel.matH = matH
	hessModel.havePerp = True
	hessModel.vecGPerpA = np.zeros(sizeX)
	hessModel.matW = np.zeros((sizeX, sizeX))
	return hessModel
# End def calcHessModel() etc.

class calcCurves_prm():
	def __init__(self):
		self.fFloorC0 = 0.0
		self.fFloorC1 = -0.01
		self.epsEigWB = 1.0e-6
		self.approxMaxStep = -1.0
		self.numVals = 100
		self.curveExpA = 2.0
		self.curveExpB = 2.0
		self.epsCurve = 1.0e-6
	def dump(self):
		msg(f'Begin calcCurves_prm().dump...')
		msg(f'self = {self}')
		msg(f'  fFloorC0 = {self.fFloorC0}')
		msg(f'  fFloorC1 = {self.fFloorC1}')
		msg(f'  epsEigWB = {self.epsEigWB}')
		msg(f'  approxMaxStep = {self.approxMaxStep}')
		msg(f'  numVals = {self.numVals}')
		msg(f'  curveExpA = {self.curveExpA}')
		msg(f'  curveExpB = {self.curveExpB}')
		msg(f'  epsCurve = {self.epsCurve}')
		msg(f'End calcCurves_prm().dump.')
class calcCurves_datOut():
	def __init__(self):
		pass
	def dump(self):
		pass
def calcCurves__eigWB__floor( f0, vecPhi, vecLambda, prm ):
	fFloor = prm.fFloorC0 + (prm.fFloorC1 * f0)
	assert (fFloor < f0)
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
	lambdaLo = prm.epsEigWB * max(vecLambda)
	if (fRes(lambdaLo) >= 0.0):
		return lambdaLo
	lambdaF = optimize.bisect( fRes, lambdaLo, lambdaHi )
	return lambdaF
def calcCurves__eigWB( f0, vecGamma0, matH, prm ):
	vecLambdaC, matPsi = np.linalg.eig(matH)
	for n in range(len(vecLambdaC)):
		assert ( np.isreal(vecLambdaC[n]) )
	vecLambda = np.real(vecLambdaC)
	vecPhi = matPsi.T @ (-vecGamma0)
	lambdaFloor = calcCurves__eigWB__floor(f0, vecPhi, vecLambda, prm)
	vecLambdaWB = vecLambda.copy()
	vecLambdaWB[vecLambdaWB < lambdaFloor] = lambdaFloor
	return vecLambda, matPsi, vecPhi, vecLambdaWB
def calcCurves( hessModel, vecYLaunch=None, vecS=None, prm=calcCurves_prm() ):
	# Parse input.
	sizeX = hessModel.matV.shape[0]
	sizeK = hessModel.matV.shape[1]
	if (None == vecYLaunch):
		vecYLaunch = np.zeros(sizeK)
	if (None == vecS):
		vecS = np.ones(sizeK)
	vecXLaunch = hessModel.vecXA + (hessModel.matV @ vecYLaunch)
	# Move to launch point, scale, and get well-behaved eigen-factorization.
	f0 = hessModel.fA + (vecYLaunch @ hessModel.vecGammaA) + ((vecYLaunch @ hessModel.matH @ vecYLaunch)/2.0)
	vecGamma0 = (hessModel.vecGammaA + (hessModel.matH @ vecYLaunch)) * vecS
	matA = (hessModel.matH * vecS).T * vecS
	matH0 = (matA.T + matA)/2.0 # Just to be safe.
	vecLambda, matPsi, vecPhi, vecLambdaWB = calcCurves__eigWB( f0, vecGamma0, matH0, prm )
	#
	matHWB = matPsi @ np.diag(vecLambdaWB) @ matPsi.T # Not needed, except for datOut.
	vecZetaNewt = vecPhi / vecLambdaWB
	vecDeltaYNewt = vecS * (matPsi @ vecZetaNewt)
	vecDeltaXNewt = hessModel.matV @ vecDeltaYNewt
	vecYNewt = vecDeltaYNewt + vecYLaunch
	vecXNewt = vecDeltaXNewt + vecXLaunch
	#
	msg(f'f0 = {f0}')
	msg(f'vecGamma0 = {vecGamma0}')
	msg(f'matH0 = ...\n{matH0}')
	msg(f'(matPsi @ np.diag(vecLambda) @ matPsi.T) - matH0 = ...\n{(matPsi @ np.diag(vecLambda) @ matPsi.T) - matH0}')
	msg(f'vecLambda = {vecLambda}')
	msg(f'matPsi = ...\n{matPsi}')
	msg(f'vecPhi = {vecPhi}')
	msg(f'vecLambdaWB = {vecLambdaWB}')
	#
	numVals = prm.numVals
	mu0 = np.min(vecLambdaWB)
	tLo = 0.0
	tHi = 1.0 # Unless...
	if (0.0 < prm.approxMaxStep):
		if (norm(vecDeltaYNewt) > prm.approxMaxStep):
			tHi = prm.approxMaxStep / norm(vecDeltaYNewt)
	tVals = tLo + ((tHi-tLo)*( 1.0 - ((1.0-(np.linspace( 0.0, 1.0, numVals )**prm.curveExpA))**prm.curveExpB) ))
	muVals = np.zeros(numVals)
	vecZetaVals = np.zeros((sizeK, numVals))
	for n in range(numVals):
		t = tVals[n]
		if ( t < prm.epsCurve ):
			muVals[n] = 0.0
			vecZetaVals[:,n] = vecPhi*(t/mu0)
		else:
			mu = mu0*((1.0/t)-1.0)
			muVals[n] = mu
			vecZetaVals[:,n] = vecPhi / (vecLambdaWB + mu)
	#
	vecDeltaYVals = (np.reshape(vecS,(sizeK,1)) * (matPsi @ vecZetaVals))
	vecDeltaXVals = hessModel.matV @ vecDeltaYVals
	vecYVals = vecDeltaYVals + np.reshape(vecYLaunch, (sizeK,1)) # Autobroadcast.
	vecXVals = vecDeltaXVals + np.reshape(vecXLaunch, (sizeX,1)) # Autobroadcast.
	#
	datOut = calcCurves_datOut()
	datOut.vecYLaunch = vecYLaunch
	datOut.vecXLaunch = vecXLaunch
	datOut.vecS = vecS
	#
	datOut.matPsi = matPsi
	datOut.vecLambda = vecLambda
	datOut.vecLambdaWB = vecLambdaWB
	datOut.mu0 = mu0
	datOut.matHWB = matHWB
	datOut.vecZetaNewt = vecZetaNewt
	datOut.vecDeltaYNewt = vecDeltaYNewt
	datOut.vecDeltaXNewt = vecDeltaXNewt
	datOut.vecYNewt = vecYNewt
	datOut.vecXNewt = vecXNewt
	#
	datOut.tVals = tVals
	datOut.muVals = muVals
	datOut.vecZetaVals = vecZetaVals
	#
	datOut.vecDeltaYVals = vecDeltaYVals
	datOut.vecDeltaXVals = vecDeltaXVals
	datOut.vecYVals = vecYVals
	datOut.vecXVals = vecXVals
	#
	datOut.dVals = np.zeros(numVals)
	datOut.fVals = np.zeros(numVals)
	datOut.vecGammaVals = np.zeros((sizeK, numVals))
	datOut.vecGPerpVals = np.zeros((sizeX, numVals))
	datOut.fWBVals = np.zeros(numVals)
	datOut.vecGammaWBVals = np.zeros((sizeK, numVals))
	for n in range(numVals):
		temp_vecY = vecYVals[:,n]
		temp_f, temp_vecGamma, temp_vecGPerp = hessModel.evalFGammaGPerpOfY(temp_vecY)
		datOut.dVals[n] = norm(temp_vecY)
		datOut.fVals[n] = temp_f
		datOut.vecGammaVals[:,n] = temp_vecGamma
		datOut.vecGPerpVals[:,n] = temp_vecGPerp
		datOut.fWBVals[n] = f0 + (temp_vecY @ vecGamma0) + ((temp_vecY @ matHWB @ temp_vecY)/2.0)
		datOut.vecGammaWBVals[:,n] = vecGamma0 + (matHWB @ temp_vecY)
	return datOut
