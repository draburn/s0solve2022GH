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
#chm.dropRelThresh /= 10.0
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
oracle_hc, oracle_hmPSD, oracle_hmWB, oracle_debug_hmLS = hessmodel.calcHessCurves(oracle_hm)
#
msg('')
msg(f'hc = {hc}...')
hc.dump()
msg('')
msg(f'oracle_hc = {oracle_hc}...')
oracle_hc.dump()

numVals = 10
muScl = min(hc.vecLambdaWB)
tVals = np.array(np.linspace(0.0,1.0,numVals))
muVals = np.zeros(numVals)
levWB_vecYVals = np.zeros((sizeK, numVals))
levWB_dVals = np.zeros(numVals)
levWB_fVals = np.zeros(numVals)
levWB_fActVals = np.zeros(numVals)
levWB_fLSVals = np.zeros(numVals)
levWB_fPSDVals = np.zeros(numVals)
levWB_fWBVals = np.zeros(numVals)
levWB_vecGammaActVals = np.zeros((sizeK, numVals))
levWB_vecGammaVals = np.zeros((sizeK, numVals))
levWB_vecGammaLSVals = np.zeros((sizeK, numVals))
levWB_vecGammaPSDVals = np.zeros((sizeK, numVals))
levWB_vecGammaWBVals = np.zeros((sizeK, numVals))
levWB_vecGPerpActVals = np.zeros((sizeX, numVals))
levWB_vecGPerpVals = np.zeros((sizeX, numVals))
levWB_vecGPerpLSVals = np.zeros((sizeX, numVals))
levWB_vecGPerpPSDVals = np.zeros((sizeX, numVals))
levWB_vecGPerpWBVals = np.zeros((sizeX, numVals))
for n in range(numVals):
	t = tVals[n]
	if ( t < 1.0e-6 ):
		muVals[n] = 0.0
		vecY = np.zeros(sizeK)
	else:
		mu = muScl*((1.0/t)-1.0)
		muVals[n] = mu
		vecY = hc.calcYLevWBOfMu(mu)
	levWB_vecYVals[:,n] = vecY
	levWB_dVals[n] = norm(vecY)
	levWB_fVals[n], levWB_vecGammaVals[:,n], levWB_vecGPerpVals[:,n] = hm.evalFGammaGPerpOfY(vecY)
	levWB_fLSVals[n], levWB_vecGammaLSVals[:,n], levWB_vecGPerpLSVals[:,n] = hmLS.evalFGammaGPerpOfY(vecY)
	levWB_fPSDVals[n], levWB_vecGammaPSDVals[:,n], levWB_vecGPerpPSDVals[:,n] = hmPSD.evalFGammaGPerpOfY(vecY)
	levWB_fWBVals[n], levWB_vecGammaWBVals[:,n], levWB_vecGPerpWBVals[:,n] = hmWB.evalFGammaGPerpOfY(vecY)
	vecX = hm.vecXA + (hm.matV @ vecY)
	msg(f'{n} / {numVals}')
	f, vecG = prob.evalFG(vecX)
	vecGamma = hm.matV.T @ vecG
	vecGPerp = vecG - (hm.matV @ vecGamma)
	levWB_fActVals[n] = f
	levWB_vecGammaActVals[:,n] = vecGamma[:]
	levWB_vecGPerpActVals[:,n] = vecGPerp[:]
	assert (danutil.reldiff(levWB_fVals[n], levWB_fLSVals[n]) < 1.0e-6)
	assert (danutil.reldiff(levWB_vecGammaVals[:,n], levWB_vecGammaLSVals[:,n]) < 1.0e-6)
	assert (danutil.reldiff(levWB_vecGPerpVals[:,n], levWB_vecGPerpLSVals[:,n]) < 1.0e-6)

import matplotlib.pyplot as plt
plt.plot(
  levWB_dVals, levWB_fVals, 'x-',
  levWB_dVals, levWB_fPSDVals, 'o-',
  levWB_dVals, levWB_fWBVals, 's-',
  levWB_dVals, levWB_fActVals, '+-' )
plt.grid(True)
plt.legend([
  'f(levWB)',
  'fLS(levWB)',
  'fPSD(levWB)',
  'fWB(levWB)',
  'fAct(levWB)' ])
plt.show()
