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
		msg(f'self = {self}')
		msg(f'  learningRate = {self.learningRate}')
		msg(f'  momentumCoeff = {self.momentumCoeff}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  genStatsDat = {self.genStatsDat}')
		msg(f'  genRecordsDat = {self.genRecordsDat}')
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
class evalSGD_storeDat():
	def __init__(self, numSteps):
		self.numSteps = numSteps
		self.matX = np.zeros((sizeX, numSteps))
		self.vecF = np.zeros(numSteps)
		self.matG = np.zeros((sizeX, numSteps))
	def dump(self):
		msg(f'self = {self}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  matX = {self.matX}')
		msg(f'  vecF = {self.vecF}')
		msg(f'  matG = {self.matG}')
def evalSGD( vecXSeed, vecPSeed, prm=evalSGD_prm() ):
	vecX = vecXSeed.copy()
	vecP = vecPSeed.copy()
	numSteps = 0
	avg_f = 0.0
	if (prm.doStats):
		statsDat = evalSGD_statsDat()
	else:
		statsDat = None
	if (prm.doStore):
		msg(f'prm.numSteps = {prm.numSteps}')
		storeDat = evalSGD_storeDat(prm.numSteps)
	else:
		storeDat = None
	for stepIndex in range(prm.numSteps):
		f, vecG = evalFG(vecX)
		numSteps += 1
		avg_f += f
		if (prm.doStats):
			statsDat.avg_vecX[:] += vecX[:]
			statsDat.var_vecX[:] += vecX[:]**2
			statsDat.avg_f += f
			statsDat.var_f += f**2
			xtg = vecX @ vecG
			statsDat.avg_xtg += xtg
			statsDat.var_xtg += xtg**2
			statsDat.avg_vecG[:] += vecG[:]
			statsDat.var_vecG[:] += vecG[:]**2
		if (prm.doStore):
			storeDat.matX[:,stepIndex] = vecX[:]
			storeDat.vecF[stepIndex] = f
			storeDat.matG[:,stepIndex] = vecG[:]
		vecP[:] = (prm.momentumCoeff * vecP[:]) - (prm.learningRate * vecG[:])
		vecX[:] += vecP[:]
	# End step loop.
	if (prm.doStats):
		statsDat.numSteps = numSteps
		statsDat.avg_vecX[:] /= numSteps
		statsDat.var_vecX[:] /= numSteps
		statsDat.avg_f /= numSteps
		statsDat.var_f /= numSteps
		statsDat.avg_xtg /= numSteps
		statsDat.var_xtg /= numSteps
		statsDat.avg_vecG[:] /= numSteps
		statsDat.var_vecG[:] /= numSteps
		statsDat.var_vecX = danutil.var(statsDat.avg_vecX, statsDat.var_vecX)
		statsDat.var_f = danutil.var(statsDat.avg_f, statsDat.var_f)
		statsDat.var_xtg = danutil.var(statsDat.avg_xtg, statsDat.var_xtg)
		statsDat.var_vecG = danutil.var(statsDat.avg_vecG, statsDat.var_vecG)
	return vecX, vecP, avg_f, statsDat, storeDat
# End evalSGD().
