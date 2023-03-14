import danutil
from danutil import msg, msgtime
import numpy as np
from numpy.linalg import norm
import demoproblem0221 as prob
#import pytorchCIFAR10demo as prob
import hessmodel0224 as hessmodel

msgtime()
vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
#
#
numSuperPts = 5
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
msg('Calling hessmodel.calcHessModel_basicOracle()...')
chmPrm.dropRelThresh /= 10.0
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
oracle_vecYLaunch = np.zeros(oracle_sizeK)
msg('')
msg('Calling hessmodel.calcHessCurves(oracle_hm)...')
oracle_hc, oracle_hmPSD, oracle_hmWB, oracle_hmLS = hessmodel.calcHessCurves(oracle_hm)
#
msg('')
msg(f'hc = {hc}...')
hc.dump()
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

numVals = 200
#tVals = np.array(np.linspace(0.0,1.0,numVals))
tVals = 1.0 - (1.0-(np.array(np.linspace(0.0,1.0,numVals))**3))**3
#
muScl = min(hc.vecLambdaWB)
muVals = np.zeros(numVals)
olevLS_dVals  = np.zeros(numVals)
olevPSD_dVals = np.zeros(numVals)
olevWB_dVals  = np.zeros(numVals)
ofloor_dVals  = np.zeros(numVals)
olevLS_fVals  = np.zeros(numVals)
olevPSD_fVals = np.zeros(numVals)
olevWB_fVals  = np.zeros(numVals)
ofloor_fVals  = np.zeros(numVals)
for n in range(numVals):
	t = tVals[n]
	if ( t == 0.0 ):
		muVals[n] = 0.0
		vecY_levLS  = np.zeros(oracle_sizeK)
		vecY_levPSD = np.zeros(oracle_sizeK)
		vecY_levWB  = np.zeros(oracle_sizeK)
		vecY_floor  = np.zeros(oracle_sizeK)
	else:
		mu = muScl*((1.0/t)-1.0)
		muVals[n] = mu
		vecY_levLS  = oracle_hc.calcYLevLSOfMu(mu)
		vecY_levPSD = oracle_hc.calcYLevPSDOfMu(mu)
		vecY_levWB  = oracle_hc.calcYLevWBOfMu(mu)
		vecY_floor  = oracle_hc.calcYFloorOfMu(mu)
	olevLS_dVals[n]  = norm(vecY_levLS)
	olevPSD_dVals[n] = norm(vecY_levPSD)
	olevWB_dVals[n]  = norm(vecY_levWB)
	ofloor_dVals[n]  = norm(vecY_floor)
	olevLS_fVals[n], _  = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_levLS)
	olevPSD_fVals[n], _ = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_levPSD)
	olevWB_fVals[n], _  = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_levWB)
	ofloor_fVals[n], _  = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_floor)

msgtime()
import matplotlib.pyplot as plt
dWBNewt = olevWB_dVals[-1]
dVizMax = 2.0*dWBNewt
plt.plot(
  olevLS_dVals[olevLS_dVals<dVizMax],   olevLS_fVals[olevLS_dVals<dVizMax],   's-',
  olevPSD_dVals[olevPSD_dVals<dVizMax], olevPSD_fVals[olevPSD_dVals<dVizMax], '+-',
  olevWB_dVals,  olevWB_fVals, '.-',
  ofloor_dVals[ofloor_dVals<dVizMax],  ofloor_fVals[ofloor_dVals<dVizMax], '-' )
plt.grid(True)
plt.legend([
  'f (oLevLS)',
  'f (oLevPSD)',
  'f (oLevWB)',
  'f (oFloor)' ])
plt.show()
