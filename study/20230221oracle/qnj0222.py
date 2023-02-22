import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm

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
# End def calcHessModel().
