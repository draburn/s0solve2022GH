import danutil
from danutil import msg, msgtime
import numpy as np
from numpy.linalg import norm
#import demoproblem0221 as prob
import pytorchCIFAR10demo as prob
import qnj0222 as qnj

vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
#
numSuperPts = 5
record_matX = np.zeros((sizeX, numSuperPts))
record_vecF = np.zeros(numSuperPts)
record_matG = np.zeros((sizeX, numSuperPts))
vecXSeed = vecX0.copy()
vecPSeed = vecP0.copy()
sgdPrm = prob.evalSGD_prm()
sgdPrm.learning_rate = 1.0e-1
sgdPrm.momentum_coefficient = 0.9
for n in range(numSuperPts):
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed, sgdPrm)
	assert ( f >= 0.0 )
	record_matX[:,n] = sgdDat.statsDat.avg_vecX[:]
	if (sgdDat.statsDat.hessfit_f > 0.0):
		record_vecF[n] = sgdDat.statsDat.hessfit_f
	else:
		record_vecF[n] = sgdDat.statsDat.avg_f
	record_matG[:,n] = sgdDat.statsDat.avg_vecG[:]
	vecXSeed[:] = vecXHarvest[:]
	vecPSeed[:] = vecPHarvest[:]
	msg(f'{n}/{numSuperPts} {f}')
# End steps loop.
msg('')
chmPrm = qnj.calcHessModel_prm()

chmPrm.dropRelThresh = 1.0e-4

nAnchor = numSuperPts-1
msg('Using qnj.calcHessModel_basic()...')
hm_basic = qnj.calcHessModel_basic(
  record_matX[:,nAnchor],
  record_vecF[nAnchor],
  record_matG[:,nAnchor],
  record_matX,
  record_vecF,
  record_matG,
  chmPrm )
msg('Using qnj.calcHessModel_basicOracle()...')
hm_oracle = qnj.calcHessModel_basicOracle(record_matX[:,nAnchor], record_matX, prob.evalFG, chmPrm)
#
msg(f'hm_basic = {hm_basic}...')
hm_basic.dump()
msg(f'hm_oracle = {hm_oracle}...')
hm_oracle.dump()
#
curveDat = qnj.calcCurves(hm_basic)
coarse_vecXVals = curveDat.coarse_vecXVals
coarse_dVals = curveDat.coarse_dVals
coarse_numPts = coarse_vecXVals.shape[1]
coarse_fActualVals = np.zeros(coarse_numPts)
coarse_vecGActualVals = np.zeros((sizeX, coarse_numPts))
for n in range(coarse_numPts):
	msg(f'{n}')
	vecX = curveDat.coarse_vecXVals[:,n]
	f, vecG = prob.evalFG(vecX)
	coarse_fActualVals[n] = f
	coarse_vecGActualVals[:,n] = vecG
curveDat_oracle = qnj.calcCurves(hm_oracle)
coracle_vecXVals = curveDat_oracle.coarse_vecXVals
coracle_dVals = curveDat_oracle.coarse_dVals
coracle_numPts = coracle_vecXVals.shape[1]
coracle_fActualVals = np.zeros(coracle_numPts)
coracle_vecGActualVals = np.zeros((sizeX, coracle_numPts))
for n in range(coracle_numPts):
	msg(f'{n}')
	vecX = curveDat_oracle.coarse_vecXVals[:,n]
	f, vecG = prob.evalFG(vecX)
	coracle_fActualVals[n] = f
	coracle_vecGActualVals[:,n] = vecG
#
import matplotlib.pyplot as plt
plt.plot( 
  curveDat.dVals, curveDat.fVals, '+-',
  curveDat.dVals, curveDat.fWBVals, 'x-',
  coarse_dVals, coarse_fActualVals, '*-',
  curveDat_oracle.dVals, curveDat_oracle.fVals, '^-',
  curveDat_oracle.dVals, curveDat_oracle.fWBVals, 'v-',
  coracle_dVals, coracle_fActualVals, 's-' )
plt.legend([
  'f model orig',
  'fWB model orig',
  'fAct along orig',
  'f model oracle',
  'fWB model oracle',
  'fAct along oracle' ])
plt.grid(True)
plt.show()
#
msg('')
msgtime()
msg('Bye.')
