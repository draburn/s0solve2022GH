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
#
msg('')
msg(f'hm = {hm}...')
hm.dump()

msg('')
msg('Calling hessmodel.calcHessCurves(hm)...')
hc, hmPSD, hmWB, hmLS = hessmodel.calcHessCurves(hm)
# Note that hmLS should match hm with a shifted anchor.
#
msg('')
msg(f'hc = {hc}...')
hc.dump()

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
levLS_fVals  = np.zeros(numVals)
levLS_fLSVals = np.zeros(numVals)
for n in range(numVals):
	msg(f' {n} / {numVals}...')
	t = tVals[n]
	if ( t == 0.0 ):
		muVals[n] = 0.0
		vecY_levLS  = hc.vecYLaunch.copy()
	else:
		mu = muScl*((1.0/t)-1.0)
		muVals[n] = mu
		vecY_levLS  = hc.calcYLevLSOfMu(mu)
	levLS_dVals[n] = norm(vecY_levLS)
	levLS_fVals[n],  _ = prob.evalFG(hm.vecXA + hm.matV @ vecY_levLS)
	levLS_fLSVals[n], _, _ = hm.evalFGammaGPerpOfY( vecY_levLS )
	msg(f'{tVals[n]}, {muVals[n]}, {levLS_dVals[n]}, {levLS_fVals[n]}, {levLS_fLSVals[n]}')
msgtime()

msg('Finding min...')
shcPrm = hessmodel.searchHessCurve_prm()
shcPrm.curveSelector = "ls"
msg(f'shcPrm = {shcPrm}...')
shcPrm.dump()
levLSMin_vecX = hessmodel.searchHessCurve( prob.evalFG, hc, shcPrm )
levLSMin_d = norm(levLSMin_vecX-hm.vecXA)
levLSMin_f, _ = prob.evalFG(levLSMin_vecX)
#
msgtime()

oracle_vecX = hessmodel.multiOracleMin(
  prob.evalFG,
  record_matX[:,nAnchor],
  record_vecF[nAnchor],
  record_matG[:,nAnchor],
  record_matX[:,0:numRecords],
  record_vecF[0:numRecords],
  record_matG[:,0:numRecords] )
oracle_d = norm( oracle_vecX - hm.vecXA )
oracle_f, _ = prob.evalFG(oracle_vecX)

import matplotlib.pyplot as plt
plt.plot( levLS_dVals, levLS_fLSVals, 'o-', markersize=15 )
plt.plot( levLS_dVals, levLS_fVals, 's-', markersize=10 )
plt.plot( levLSMin_d, levLSMin_f, 'x', markersize=20 )
plt.plot( oracle_d, oracle_f, '*', markersize=25 )
plt.grid(True)
plt.legend([
  'fLS (curve)',
  'f (curve)',
  'f (curveMin)',
  'f (oracle)' ])
plt.xlabel('||delta||')
plt.ylabel('F')
plt.show()
