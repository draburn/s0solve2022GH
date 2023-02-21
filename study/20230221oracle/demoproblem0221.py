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

def evalF(vecX):
	vecD = vecX - vecXC
	vecG = matH @ vecD
	f = fC + ((vecG @ vecD)/2.0)
	return f
# End def evalF().

def evalFG(vecX):
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
		self.calcDatOut = True
class evalSGD_datOut():
	def __init__(self):
		self.numSteps = 0
		self.avg_vecX = np.zeros(sizeX)
		self.var_vecX = np.zeros(sizeX)
		self.avg_f = 0.0
		self.var_f = 0.0
		self.avg_vecG = np.zeros(sizeX)
		self.var_vecG = np.zeros(sizeX)
def evalSGD(vecXSeed, vecPSeed, prm=evalSGD_prm()):
	vecX = vecXSeed.copy()
	vecP = vecPSeed.copy()
	numSteps = 0
	avg_f = 0.0
	if (prm.calcDatOut):
		datOut = evalSGD_datOut()
	for stepIndex in range(prm.numSteps):
		if (prm.calcDatOut):
			datOut.avg_vecX[:] += vecX[:]
			datOut.var_vecX[:] += vecX[:]**2
		f, vecG = evalFG(vecX)
		vecP[:] = (prm.momentumCoeff * vecP[:]) - (prm.learningRate * vecG[:])
		vecX[:] += vecP[:]
		numSteps += 1
		avg_f += f
		if (prm.calcDatOut):
			datOut.avg_f += f
			datOut.var_f += f**2
			datOut.avg_vecG[:] += vecG[:]
			datOut.var_vecG[:] += vecG[:]**2
	# End step loop.
	if (prm.calcDatOut):
		datOut.numSteps = numSteps
		datOut.avg_vecX[:] /= numSteps
		datOut.var_vecX[:] /= numSteps
		datOut.avg_f /= numSteps
		datOut.var_f /= numSteps
		datOut.avg_vecG[:] /= numSteps
		datOut.var_vecG[:] /= numSteps
		datOut.var_f = danutil.var( datOut.avg_f, datOut.var_f )
		datOut.var_vecX = danutil.var( datOut.avg_vecX, datOut.var_vecX )
		datOut.var_vecG = danutil.var( datOut.avg_vecG, datOut.var_vecG )
		return vecX, vecP, avg_f, datOut
	else:
		return vecX, vecP, avg_f, None
# End evalSGD().
