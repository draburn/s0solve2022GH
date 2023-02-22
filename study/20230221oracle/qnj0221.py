import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm
#from scipy import optimize as scipy_optimize

msg('*** WARNING: THIS FILE IS DEPRECATED. ***')
exit()

class hessModelType():
	def __init(self):
		self.matV = None
		self.vecXA = None
		self.fA = None
		self.vecGammaA = None
		self.matH = None
		self.vecGAPerp = None
		self.matW = None
	def dump(self):
		msg(f'self = {self}')
		msg(f'  matV = ...\n{self.matV}')
		msg(f'  vecXA = {self.vecXA}')
		msg(f'  fA = {self.fA}')
		msg(f'  vecGammaA = {self.vecGammaA}')
		msg(f'  matH = ...\n{self.matH}')
		msg(f'  matW = ...\n{self.matW}')
	def evalFGGamma(self, vecX):
		vecY = self.matV.T @ ( vecX - self.vecXA )
		vecHY = self.matH @ vecY
		vecGamma = self.vecGammaA + vecHY
		vecG = self.matV @ vecGamma + (self.matW @ vecY) + self.vecGAPerp
		f = self.fA + (vecY @ self.vecGammaA) + ((vecY @ vecHY)/2.0)
		return f, vecG, vecGamma
class calcHessModel_prm():
	def __init__(self):
		self.dropRelThresh = 1.0e-1
		self.dropAbsThresh = 1.0e-12
		self.epsFD = 1.0e-4
		self.fdOrder = 1
	def dump(self):
		msg(f'self = {self}')
		msg(f'  dropRelThresh = {self.dropRelThresh}')
		msg(f'  dropAbsThresh = {self.dropAbsThresh}')
		msg(f'  epsFD = {self.epsFD}')
		msg(f'  fdOrder = {self.fdOrder}')
#def calcHessModel( calcMethod, vecXAnchor, fAnchor, vecGAnchor, record_matX, record_vecF, record_matG, prm=calcHessModel_prm ):
#	if ( "basic" == calcMethod ):
#		return calcHessModel__basic(vecXAnchor, fAnchor, vecGAnchor, record_matX, record_matG, calcHessModel_prm)
#	if ( "fullspace" == calcMethod ):
#		return calcHessModel__fullspace(vecXAnchor, calcHessModel_prm)
#	if ( "basicOracle" == calcMethod ):
#		return calcHessModel__basicOracle(vecXAnchor, record_matX, calcHessModel_prm)
#	msg(f'ERROR: Invalid value of prm.calcMethod ({prm.calcMethod}).')
#	return None
def calcHessModel_basic( vecXAnchor, fAnchor, vecGAnchor, record_matX, record_matG, prm=calcHessModel_prm() ):
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
	hessModel.matW = matW
	return hessModel
# End calcHessModel_basic().
def calcHessModel_basicOracle( vecXAnchor, record_matX, funch_evalFG, prm=calcHessModel_prm() ):
	sizeX = vecXAnchor.shape[0]
	matD = record_matX - np.reshape(vecXAnchor, (sizeX,1)) # Autobroadcast
	matV, vecKeep = danutil.utorthdrop(matD, prm.dropRelThresh, prm.dropAbsThresh)
	sizeK = matV.shape[1]
	msg(f'matV.T @ matV = ...\n{matV.T @ matV}')
	#
	fA, vecGA = funch_evalFG(vecXAnchor)
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

	vecGammaA = matV.T @ vecGA
	msg(f'  vecGA = {vecGA}')
	msg(f'  vecGammaA = {vecGammaA}')
	
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV
	hessModel.fA = fA
	hessModel.vecGammaA = matV.T @ vecGA
	hessModel.matH = matH
	hessModel.matW = matW
	hessModel.vecGAPerp = vecGA - (matV @ vecGammaA)
	return hessModel
# End def calcHessModel_basicOracle().
def calcHessModel_fullspace( vecXAnchor, funch_evalFG, prm=calcHessModel_prm() ):
	sizeX = vecXAnchor.shape[0]
	matV = np.eye(sizeX)
	sizeK = sizeX
	#
	f0, vecG0 = funch_evalFG(vecX0)
	matA = np.zeros((sizeX, sizeX))
	if (1 == prm.fdOrder):
		for k in range(sizeK):
			fP, vecGP = funch_evalFG(vecX0 + (prm.epsFD*matV[:,k]))
			matA[:,k] = (vecGP-vecG0)/prm.epsFD
	elif (2 == prm.fdOrder):
		for k in range(sizeK):
			fP, vecGP = funch_evalFG(vecX0 + (prm.epsFD*matV[:,k]))
			fM, vecGM = funch_evalFG(vecX0 - (prm.epsFD*matV[:,k]))
			matA[:,k] = (vecGP-vecGM)/(2.0*prm.epsFD)
	matH = ( matA.T + matA )/2.0
	#
	hessModel = hessModelType()
	hessModel.vecXA = vecXAnchor.copy()
	hessModel.matV = matV
	hessModel.fA = f0
	hessModel.vecGammaA = matV.T @ vecG0
	hessModel.matH = matH
	hessModel.matW = np.zeros((sizeX, sizeX))
	return hessModel
