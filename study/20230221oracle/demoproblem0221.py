import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm

# Set base params; modify these as needed.
reportInit = True
prngSeed = 0
sizeX = 5
#
np.random.seed(prngSeed)
fC = float(abs(np.random.randn(1))) # Leave uncommented to preserve PRNG state.
fC = 0.0
vecXC = np.random.randn(sizeX)
matA = np.random.randn(sizeX, sizeX)
matH = matA.T @ matA
if (reportInit):
	msg(f'Initialized problem...')
	msg(f'  prngSeed = {prngSeed}')
	msg(f'  sizeX = {sizeX}')
	msg(f'  fC = {fC:0.18E}')

def genVecX0():
	return np.zeros(sizeX)

def evalF( vecX ):
	vecD = vecX - vecXC
	vecG = matH @ vecD
	f = fC + ((vecG @ vecD)/2.0)
	return f
# End def evalF().

def evalFG( vecX ):
	vecD = vecX - vecXC
	vecG = matH @ vecD
	f = fC + ((vecG @ vecD)/2.0)
	return f, vecG
# End def evalFG().

class evalSGD_prm():
	def __init__(self):
		self.learningRate = 1.0e-3
		self.momentumCoeff = 0.9
		self.numSteps = 10
		self.doStats = True
		self.doStore = False
	def dump(self):
		msg(f'Begin evalSGD_prm.dump()...')
		msg(f'self = {self}')
		msg(f'  learningRate = {self.learningRate}')
		msg(f'  momentumCoeff = {self.momentumCoeff}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  genStatsDat = {self.genStatsDat}')
		msg(f'  genRecordsDat = {self.genRecordsDat}')
		msg(f'End evalSGD_prm.dump().')
class evalSGD_statsDat():
	def __init__(self):
		self.numSteps = 0
		self.avg_vecX = np.zeros(sizeX)
		self.var_vecX = np.zeros(sizeX)
		self.avg_f = 0.0
		self.var_f = 0.0
		self.avg_xtg = 0.0
		self.var_xtg = 0.0
		self.avg_vecG = np.zeros(sizeX)
		self.var_vecG = np.zeros(sizeX)
	def dump(self):
		msg(f'Begin evalSGD_statsDat.dump()...')
		msg(f'self = {self}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  avg_vecX = {self.avg_vecX}')
		msg(f'  var_vecX = {self.var_vecX}')
		msg(f'  avg_f = {self.avg_f}')
		msg(f'  var_f = {self.var_f}')
		msg(f'  avg_xtg = {self.avg_xtg}')
		msg(f'  var_xtg = {self.var_xtg}')
		msg(f'  avg_vecG = {self.avg_vecG}')
		msg(f'  var_vecG = {self.var_vecG}')
		msg(f'End evalSGD_statsDat.dump().')
	def absorb(self, vecX, f, vecG):
		self.numSteps +=1
		self.avg_vecX[:] += vecX[:]
		self.var_vecX[:] += vecX[:]**2
		self.avg_f += f
		self.var_f += f**2
		xtg = vecX @ vecG
		self.avg_xtg += xtg
		self.var_xtg += xtg**2
		self.avg_vecG[:] += vecG[:]
		self.var_vecG[:] += vecG[:]**2
	def finalize(self):
		numSteps = self.numSteps
		self.avg_vecX[:] /= numSteps
		self.var_vecX[:] /= numSteps
		self.avg_f /= numSteps
		self.var_f /= numSteps
		self.avg_xtg /= numSteps
		self.var_xtg /= numSteps
		self.avg_vecG[:] /= numSteps
		self.var_vecG[:] /= numSteps
		self.var_vecX = danutil.var(self.avg_vecX, self.var_vecX)
		self.var_f = danutil.var(self.avg_f, self.var_f)
		self.var_xtg = danutil.var(self.avg_xtg, self.var_xtg)
		self.var_vecG = danutil.var(self.avg_vecG, self.var_vecG)
class evalSGD_storeDat():
	def __init__(self, storageSize):
		self.numSteps = 0
		self.matX = np.zeros((sizeX, storageSize))
		self.vecF = np.zeros(storageSize)
		self.matG = np.zeros((sizeX, storageSize))
	def dump(self):
		msg(f'Begin evalSGD_storeDat.dump()...')
		msg(f'self = {self}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  matX = ...\n{self.matX}')
		msg(f'  vecF = {self.vecF}')
		msg(f'  matG = ...\n{self.matG}')
		msg(f'End evalSGD_storeDat.dump().')
	def absorb(self, vecX, f, vecG):
		self.matX[:,self.numSteps] = vecX[:]
		self.vecF[self.numSteps] = f
		self.matG[:,self.numSteps] = vecG[:]
		self.numSteps += 1
	def finalize(self):
		pass
class evalSGD_datOut():
	def __init__(self, doStats, doStore, storageSize):
		self.statsDat = None # Unless...
		self.storeDat = None # Unless...
		if (doStats):
			self.statsDat = evalSGD_statsDat()
		if (doStore):
			self.storeDat = evalSGD_storeDat(storageSize)
	def dump(self):
		msg(f'Begin evalSGD_datOut.dump()...')
		msg(f'self = {self}')
		msg(f'  statsDat = {self.statsDat}')
		if (not(None == self.statsDat)):
			self.statsDat.dump()
		msg(f'  storeDat = {self.storeDat}')
		if (not(None == self.storeDat)):
			self.storeDat.dump()
		msg(f'End evalSGD_datOut.dump().')
	def absorb(self, vecX, f, vecG):
		if (not(None == self.statsDat)):
			self.statsDat.absorb(vecX, f, vecG)
		if (not(None == self.storeDat)):
			self.storeDat.absorb(vecX, f, vecG)
	def finalize(self):
		if (not(None == self.statsDat)):
			self.statsDat.finalize()
		if (not(None == self.storeDat)):
			self.storeDat.finalize()
def evalSGD( vecXSeed, vecPSeed=np.zeros(sizeX), prm=evalSGD_prm() ):
	vecX = vecXSeed.copy()
	vecP = vecPSeed.copy()
	stepCount = 0
	avg_f = 0.0
	datOut = evalSGD_datOut(prm.doStats, prm.doStore, prm.numSteps)
	for stepIndex in range(prm.numSteps):
		f, vecG = evalFG(vecX)
		datOut.absorb( vecX, f, vecG )
		stepCount += 1
		avg_f += f
		vecP[:] = (prm.momentumCoeff * vecP[:]) - (prm.learningRate * vecG[:])
		vecX[:] += vecP[:]
	datOut.finalize()
	avg_f /= stepCount
	return vecX, vecP, avg_f, datOut
# End evalSGD().
