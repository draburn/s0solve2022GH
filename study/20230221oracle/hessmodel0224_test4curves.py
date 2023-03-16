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

numVals = 20
#tVals = np.array(np.linspace(0.0,1.0,numVals))
#tVals = 1.0 - (1.0-(np.array(np.linspace(0.0,1.0,numVals))**3))**3
tExpLo = 1
tExpHi = 1
tVals = 1.0 - (1.0-(np.array(np.linspace(0.0,1.0,numVals))**tExpLo))**tExpHi
tVals /= 100
#
muScl = min(hc.vecLambdaWB)
muVals = np.zeros(numVals)
levLS_dVals  = np.zeros(numVals)
levPSD_dVals = np.zeros(numVals)
levWB_dVals  = np.zeros(numVals)
floor_dVals  = np.zeros(numVals)
levLS_fVals  = np.zeros(numVals)
levPSD_fVals = np.zeros(numVals)
levWB_fVals  = np.zeros(numVals)
floor_fVals  = np.zeros(numVals)
levWB_fLSVals = np.zeros(numVals)
olevLS_dVals  = np.zeros(numVals)
olevPSD_dVals = np.zeros(numVals)
olevWB_dVals  = np.zeros(numVals)
ofloor_dVals  = np.zeros(numVals)
olevLS_fVals  = np.zeros(numVals)
olevPSD_fVals = np.zeros(numVals)
olevWB_fVals  = np.zeros(numVals)
ofloor_fVals  = np.zeros(numVals)
olevWB_ofLSVals = np.zeros(numVals)
for n in range(numVals):
	msg(f' {n} / {numVals}...')
	t = tVals[n]
	if ( t == 0.0 ):
		muVals[n] = 0.0
		vecY_levLS  = np.zeros(sizeK)
		vecY_levPSD = np.zeros(sizeK)
		vecY_levWB  = np.zeros(sizeK)
		vecY_floor  = np.zeros(sizeK)
	else:
		mu = muScl*((1.0/t)-1.0)
		muVals[n] = mu
		vecY_levLS  = hc.calcYLevLSOfMu(mu)
		vecY_levPSD = hc.calcYLevPSDOfMu(mu)
		vecY_levWB  = hc.calcYLevWBOfMu(mu)
		vecY_floor  = hc.calcYFloorOfMu(mu)
	levLS_dVals[n]  = norm(vecY_levLS)
	levPSD_dVals[n] = norm(vecY_levPSD)
	levWB_dVals[n]  = norm(vecY_levWB)
	floor_dVals[n]  = norm(vecY_floor)
	levLS_fVals[n],  _ = prob.evalFG(hm.vecXA + hm.matV @ vecY_levLS)
	levPSD_fVals[n], _ = prob.evalFG(hm.vecXA + hm.matV @ vecY_levPSD)
	levWB_fVals[n],  _ = prob.evalFG(hm.vecXA + hm.matV @ vecY_levWB)
	floor_fVals[n],  _ = prob.evalFG(hm.vecXA + hm.matV @ vecY_floor)
	levWB_fLSVals[n], _, _ = hm.evalFGammaGPerpOfY( vecY_levWB )
	#
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
	olevLS_fVals[n],  _ = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_levLS)
	olevPSD_fVals[n], _ = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_levPSD)
	olevWB_fVals[n],  _ = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_levWB)
	ofloor_fVals[n],  _ = prob.evalFG(oracle_hm.vecXA + oracle_hm.matV @ vecY_floor)
	olevWB_ofLSVals[n], _, _ = oracle_hm.evalFGammaGPerpOfY( vecY_levWB )
msgtime()

msg('Finding min...')
shcPrm = hessmodel.searchHessCurve_prm()
#
levLSMin_vecX = hm.vecXA.copy()
levPSDMin_vecX = hm.vecXA.copy()
levWBMin_vecX = hm.vecXA.copy()
floorMin_vecX = hm.vecXA.copy()
if (True):
	shcPrm.curveSelector = "ls"
	msg(f'shcPrm = {shcPrm}...')
	shcPrm.dump()
	#levLSMin_vecX = hessmodel.searchHessCurve( prob.evalFG, hc, shcPrm )
	levLSMin_d = norm(levLSMin_vecX-hm.vecXA)
	levLSMin_f, _ = prob.evalFG(levLSMin_vecX)
	#
	shcPrm.curveSelector = "psd"
	msg(f'shcPrm = {shcPrm}...')
	shcPrm.dump()
	#levPSDMin_vecX = hessmodel.searchHessCurve( prob.evalFG, hc, shcPrm )
	levPSDMin_d = norm(levPSDMin_vecX-hm.vecXA)
	levPSDMin_f, _ = prob.evalFG(levPSDMin_vecX)
#
shcPrm.curveSelector = "wb"
msg(f'shcPrm = {shcPrm}...')
shcPrm.dump()
#levWBMin_vecX = hessmodel.searchHessCurve( prob.evalFG, hc, shcPrm )
levWBMin_d = norm(levWBMin_vecX-hm.vecXA)
levWBMin_f, _ = prob.evalFG(levWBMin_vecX)
#
shcPrm.curveSelector = "floor"
msg(f'shcPrm = {shcPrm}...')
shcPrm.dump()
#floorMin_vecX = hessmodel.searchHessCurve( prob.evalFG, hc, shcPrm )
floorMin_d = norm(floorMin_vecX-hm.vecXA)
floorMin_f, _ = prob.evalFG(floorMin_vecX)
#
oracle_vecX = hessmodel.multiSearchHessCurve( prob.evalFG, oracle_hc )
oracle_d = norm(oracle_vecX-oracle_hm.vecXA)
oracle_f, _ = prob.evalFG(oracle_vecX)
#
msgtime()

import matplotlib.pyplot as plt
dWBNewt = olevWB_dVals[-1]
dVizMax = 2.0*dWBNewt
plt.plot( levWB_dVals, levWB_fLSVals, 'p-', markersize=23 )
plt.plot( levLS_dVals[ levLS_dVals<dVizMax ], levLS_fVals[ levLS_dVals<dVizMax ], 's-', markersize=21 )
plt.plot( levPSD_dVals[levPSD_dVals<dVizMax], levPSD_fVals[levPSD_dVals<dVizMax], 'o-', markersize=19 )
plt.plot( levWB_dVals, levWB_fVals, '^-', markersize=17 )
plt.plot( floor_dVals[ floor_dVals<dVizMax], floor_fVals[ floor_dVals<dVizMax ], 'v-', markersize=15 )
plt.plot( olevWB_dVals, olevWB_ofLSVals, 'p-', markersize=13 )
plt.plot( olevLS_dVals[ olevLS_dVals<dVizMax ], olevLS_fVals[ olevLS_dVals<dVizMax ], 's-', markersize=11 )
plt.plot( olevPSD_dVals[olevPSD_dVals<dVizMax], olevPSD_fVals[olevPSD_dVals<dVizMax], 'o-', markersize=9 )
plt.plot( olevWB_dVals, olevWB_fVals, '^-', markersize=7 )
plt.plot( ofloor_dVals[ ofloor_dVals<dVizMax ], ofloor_fVals[ ofloor_dVals<dVizMax ], 'v-', markersize=5 )
plt.plot( 0, 0 )
plt.plot( levLSMin_d, levLSMin_f, '*', markersize=31 )
plt.plot( levPSDMin_d, levPSDMin_f, '*', markersize=29 )
plt.plot( levWBMin_d, levWBMin_f, '*', markersize=27 )
plt.plot( floorMin_d, floorMin_f, '*', markersize=25 )
plt.plot( oracle_d, oracle_f, '*', markersize=20 )
plt.grid(True)
plt.legend([
  'fLS (LevWB)',
  'f (LevLS)',
  'f (LevPSD)',
  'f (LevWB)',
  'f (Floor)',
  'ofLS (oLevWB)',
  'f (oLevLS)',
  'f (oLevPSD)',
  'f (oLevWB)',
  'f (oFloor)',
  '-',
  'f (LevLSMin)',
  'f (LevPSDMin)',
  'f (LevWBMin)',
  'f (FloorMin)',
  'f (oracleMin)' ])
plt.xlabel('||delta||')
plt.ylabel('F')
plt.show()
