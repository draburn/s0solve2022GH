import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm
import demoproblem0221 as prob
import qnj0221 as qnj

msg('*** WARNING: THIS FILE IS DEPRECATED. ***')
exit()

vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
#
numSteps = 5
record_matX = np.zeros((sizeX, numSteps))
record_vecF = np.zeros(numSteps)
record_matG = np.zeros((sizeX, numSteps))
vecXSeed = vecX0.copy()
vecPSeed = vecP0.copy()
for n in range(numSteps):
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed)
	record_matX[:,n] = sgdDat.statsDat.avg_vecX[:]
	record_vecF[n] = sgdDat.statsDat.avg_f
	record_matG[:,n] = sgdDat.statsDat.avg_vecG[:]
	vecXSeed[:] = vecXHarvest[:]
	vecPSeed[:] = vecPHarvest[:]
#
#record_matX = np.eye(sizeX)
#record_matX = np.diag(1+np.array(range(sizeX)))
#record_matX = np.random.randn(sizeX, sizeX)
#vecXAnchor = np.zeros(sizeX)
vecXAnchor = record_matX[:,0]
msg('')
#msg('Using qnj.calcHessModel_basic()...')
#hm = qnj.calcHessModel_basic(record_matX[:,0], record_vecF[0], record_matG[:,0], record_matX, record_matG)
msg('Using qnj.calcHessModel_basicOracle()...')
hm = qnj.calcHessModel_basicOracle(vecXAnchor, record_matX, prob.evalFG)
msg(f'hm = {hm}...')
hm.dump()
vecY = np.array(range(hm.matV.shape[1]))
#vecY = np.zeros(hm.matV.shape[1])
vecX = hm.vecXA + (hm.matV @ vecY)
#vecX = record_matX[:,0]
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
exit()



#
msg('')
msg('Using qnj.calcHessModel_basicOracle()...')
hm = qnj.calcHessModel_basicOracle(record_matX[:,0], record_matX, prob.evalFG)
msg(f'hm.matH = {hm.matH}')
msg(f'hm.matV.T @ prob.matH @ hm.matV = {hm.matV.T @ prob.matH @ hm.matV}')
f, vecG, vecGamma = hm.evalFGGamma(record_matX[:,-1])
msg(f'hm.evalFGGamma(record_matX[:,-1])[0] = {f}')
msg(f'hm.evalFGGamma(record_matX[:,-1])[1] = {vecG}')
msg(f'hm.evalFGGamma(record_matX[:,-1])[2] = {vecGamma}')
#
msg('')
msg('Bye.')
exit()
vecXA, matV, fA, vecGSSA, matHSS = qnj.genHessModelFD(record_matX[:,0], record_matX, prob.evalFG)
#msg(f'vecXA = {vecXA}')
#msg(f'matV = {matV}')
#msg(f'fA = {fA}')
#msg(f'vecGSSA = {vecGSSA}')
msg(f'matHSS = {matHSS}')
#msg(f'matV.T @ prob.matH @ matV = {matV.T @ prob.matH @ matV}')
#
vecXA, matV, fA, vecGSSA, matHSS = qnj.genHessModelBasic(
  record_matX[:,0],
  record_vecF[0],
  record_matG[:,0],
  record_matX,
  record_matG)
#msg(f'vecXA = {vecXA}')
#msg(f'matV = {matV}')
#msg(f'fA = {fA}')
#msg(f'vecGSSA = {vecGSSA}')
msg(f'matHSS = {matHSS}')
#msg(f'matV.T @ prob.matH @ matV = {matV.T @ prob.matH @ matV}')


msg('Bye.')
exit()

_, f0, vecG0, matH = qnj.genHessModel_fullSpace(vecX0, prob.evalFG)
msg(f'matH = {matH}')
exit()
