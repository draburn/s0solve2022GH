import danutil
from danutil import msg, msgtime
import numpy as np
from numpy.linalg import norm
#import demoproblem0221 as prob
import pytorchCIFAR10demo as prob
import hessmodel0224 as hessmodel

msgtime()
vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
#
#
numSuperPts = 100
maxNumRecords = numSuperPts
#
record_matX = np.zeros((sizeX, maxNumRecords))
record_vecF = np.zeros(maxNumRecords)
record_matG = np.zeros((sizeX, maxNumRecords))
numRecords = 0
vecXSeed = vecX0.copy()
vecPSeed = vecP0.copy()
sgdPrm = prob.evalSGD_prm()
sgdPrm.learningRate = 1.0e-2
sgdPrm.momentumCoefficient = 0.9
msg(f'sgdPrm = {sgdPrm}...')
sgdPrm.dump()
for n in range(numSuperPts):
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed, sgdPrm)
	assert ( f >= 0.0 )
	record_matX[:,1:] = record_matX[:,:-1]
	record_matG[:,1:] = record_matG[:,:-1]
	record_vecF[1:] = record_vecF[:-1]
	if ( numRecords < maxNumRecords ):
		numRecords += 1
	record_matX[:,0] = sgdDat.statsDat.avg_vecX[:]
	if (sgdDat.statsDat.hessfit_f > 0.0):
		record_vecF[0] = sgdDat.statsDat.hessfit_f
	else:
		record_vecF[0] = f
	record_matG[:,0] = sgdDat.statsDat.avg_vecG[:]
	
	vecXSeed[:] = vecXHarvest[:]
	vecPSeed[:] = vecPHarvest[:]
	msg(f'{n}/{numSuperPts} {f} {sgdDat.statsDat.hessfit_f}')
# End steps loop.

chmPrm = hessmodel.calcHessModel_prm()
nAnchor = 0

msg('')
msg('Calling hessmodel.calcHessModel()...')
msg(f'chmPrm = {chmPrm}...')
chmPrm.dump()
hm = hessmodel.calcHessModel(
  record_matX[:,nAnchor],
  record_vecF[nAnchor],
  record_matG[:,nAnchor],
  record_matX[:,0:numRecords],
  record_vecF[0:numRecords],
  record_matG[:,0:numRecords],
  chmPrm )
sizeK = hm.matV.shape[1]
msg(f'sizeK = {sizeK}')
msg('')
calcOracle = False
if (calcOracle):
	msg('Calling hessmodel.calcHessModel_basicOracle()...')
	chmPrm.dropRelThresh /= 10.0
	#chmPrm.epsRelFD = 1.0e-1 # DRaburn 2023-03-13-2100
	msg(f'chmPrm = {chmPrm}...')
	chmPrm.dump()
	oracle_hm = hessmodel.oracle_calcHessModel(prob.evalFG, record_matX[:,nAnchor], record_matX[:,0:numRecords], chmPrm)
	oracle_sizeK = oracle_hm.matV.shape[1]
	msg(f'oracle_sizeK = {oracle_sizeK}')
	#
	msg('')
	msg(f'hm = {hm}...')
	hm.dump()
	msg('')
	msg(f'oracle_hm = {oracle_hm}...')
	oracle_hm.dump()
if (False):
	msg('')
	msg(f'hm.matV.T @ prob.matH @ hm.matV = ...\n{hm.matV.T @ prob.matH @ hm.matV}')

msg('')
msg('Calling hessmodel.calcHessCurves(hm)...')
hc, hmPSD, hmWB, hmLS = hessmodel.calcHessCurves(hm)
# Note that hmLS should match hm with a shifted anchor.
msg('')
msg(f'hc = {hc}...')
hc.dump()
if (calcOracle):
	oracle_vecYLaunch = np.zeros(oracle_sizeK)
	msg('')
	msg('Calling hessmodel.calcHessCurves(oracle_hm)...')
	oracle_hc, oracle_hmPSD, oracle_hmWB, oracle_hmLS = hessmodel.calcHessCurves(oracle_hm)
	msg('')
	msg(f'oracle_hc = {oracle_hc}...')
	oracle_hc.dump()

if (False):
	vecXOptim = hessmodel.searchHessCurve(prob.evalFG, hc)
	fOptim, vecGOptim = prob.evalFG(vecXOptim)
	oracle_vecXOptim = hessmodel.searchHessCurve(prob.evalFG, oracle_hc)
	oracle_fOptim, oracle_vecGOptim = prob.evalFG(oracle_vecXOptim)
	msg(f'min f = {min(record_vecF)}, fOptim = {fOptim}, oracle_fOptim = {oracle_fOptim}')
	danutil.bye()

numVals = 100
#tVals = np.array(np.linspace(0.0,1.0,numVals))
#tVals = 1.0 - (1.0-(np.array(np.linspace(0.0,1.0,numVals))**3))**3
tExpLo = 2
tExpHi = 1
tVals = 1.0 - (1.0-(np.array(np.linspace(0.0,1.0,numVals))**tExpLo))**tExpHi
#tVals /= 100
tMax_cauchy = 1.0
#
muScl = min(hc.vecLambdaWB)
muVals = np.zeros(numVals)
levLS_dVals = np.zeros(numVals)
levLS_fLSVals = np.zeros(numVals)
levLS_fVals = np.zeros(numVals)
floor_dVals = np.zeros(numVals)
floor_fVals = np.zeros(numVals)
floor_fWBVals = np.zeros(numVals)
ccywb_dVals = np.zeros(numVals)
ccywb_fVals = np.zeros(numVals)
for n in range(numVals):
	msg(f' {n} / {numVals}...')
	t = tVals[n]
	if ( t == 0.0 ):
		muVals[n] = 0.0
		vecY_levLS  = np.zeros(sizeK)
		vecY_floor  = np.zeros(sizeK)
	else:
		mu = muScl*((1.0/t)-1.0)
		muVals[n] = mu
		vecY_levLS  = hc.calcYLevLSOfMu(mu)
		vecY_floor  = hc.calcYFloorOfMu(mu)
	vecY_ccywb = hc.calcYCauchyWBOfT( tMax_cauchy*(n/(numVals-1.0))**tExpLo )
	vecX_levLS = hm.vecXA + hm.matV @ vecY_levLS
	vecX_floor = hm.vecXA + hm.matV @ vecY_floor
	vecX_ccywb = hm.vecXA + hm.matV @ vecY_ccywb
	levLS_dVals[n] = norm(vecX_levLS - hm.vecXA)
	floor_dVals[n] = norm(vecX_floor - hm.vecXA)
	ccywb_dVals[n] = norm(vecX_ccywb - hm.vecXA)
	levLS_fVals[n],  _ = prob.evalFG(vecX_levLS)
	floor_fVals[n],  _ = prob.evalFG(vecX_floor)
	ccywb_fVals[n],  _ = prob.evalFG(vecX_ccywb)
	levLS_fLSVals[n], _, _ = hm.evalFGammaGPerpOfY( vecY_levLS )
	floor_fWBVals[n], _, _ = hmWB.evalFGammaGPerpOfY( vecY_floor )
msgtime()

import matplotlib.pyplot as plt
dVizMax = 2.0*floor_dVals[-1]
plt.plot( levLS_dVals[levLS_dVals<dVizMax], levLS_fLSVals[levLS_dVals<dVizMax], 'o-', markersize=28 )
plt.plot( levLS_dVals[levLS_dVals<dVizMax], levLS_fVals[levLS_dVals<dVizMax], 'o-', markersize=24 )
plt.plot( floor_dVals, floor_fWBVals, 's-', markersize=22 )
plt.plot( floor_dVals, floor_fVals, 's-', markersize=18 )
plt.plot( ccywb_dVals, ccywb_fVals, '^-', markersize=16 )
plt.grid(True)
plt.legend([
  'LevLS, LS model',
  'LevLS, actual',
  'Floor, WB model',
  'Floor, actual',
  'Cauchy, actual' ])
plt.xlabel('||delta||')
plt.ylabel('F')
plt.show()
