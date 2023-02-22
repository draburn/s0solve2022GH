import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm
import demoproblem0221 as prob
import qnj0222 as qnj

vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
#
numSteps = 100
record_matX = np.zeros((sizeX, numSteps))
record_vecF = np.zeros(numSteps)
record_matG = np.zeros((sizeX, numSteps))
vecXSeed = vecX0.copy()
vecPSeed = vecP0.copy()
for n in range(numSteps):
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed)
	record_matX[:,n] = sgdDat.statsDat.avg_vecX[:]
	#record_vecF[n] = sgdDat.statsDat.avg_f
	record_vecF[n] = sgdDat.statsDat.hessfit_f
	record_matG[:,n] = sgdDat.statsDat.avg_vecG[:]
	vecXSeed[:] = vecXHarvest[:]
	vecPSeed[:] = vecPHarvest[:]
# End steps loop.
msg('')
chmPrm = qnj.calcHessModel_prm()
#chmPrm.dropRelThresh = 1.0e-4
if (True):
	msg('Using qnj.calcHessModel_basic()...')
	n = 0
	hm = qnj.calcHessModel_basic(
	  record_matX[:,n],
	  record_vecF[n],
	  record_matG[:,n],
	  record_matX,
	  record_vecF,
	  record_matG,
	  chmPrm )
elif (False):
	msg('Using qnj.calcHessModel_basicOracle()...')
	n = 0
	hm = qnj.calcHessModel_basicOracle(record_matX[:,n], record_matX, prob.evalFG, chmPrm)
else:
	msg('Using qnj.calcHessModel_fullspace()...')
	n = 0
	hm = qnj.calcHessModel_fullspace(record_matX[:,n], prob.evalFG, chmPrm)
msg(f'hm = {hm}...')
hm.dump()
vecY = np.array(range(hm.matV.shape[1]))
vecX = hm.vecXA + (hm.matV @ vecY)
f, vecG = prob.evalFG(vecX)
msg(f'Using prob.evalFG(vecX)...')
msg(f'  f = {f}')
msg(f'  vecG = {vecG}')
msg(f'  vecGamma = {hm.matV.T @ vecG}')
msg(f'  matH = ...\n{hm.matV.T @ prob.matH @ hm.matV}')
f, vecG, vecGamma = hm.evalFGGamma(vecX)
msg(f'Using hm.evalFGGamma(vecX)...')
msg(f'  f = {f}')
msg(f'  vecG = {vecG}')
msg(f'  vecGamma = {vecGamma}')
msg(f'  matH = ...\n{hm.matH}')
#
msg('')
msg('Bye.')
