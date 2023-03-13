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
numSuperPts = 50
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
	msg('')
	msg('')
	msg('Spot check (to compare with old, qnj0222_test2 code)...')
	msg('')
	msg('  Along levWB...')
	#mu  =  7.967038104510303
	#mu = 1.886649452148539685E-01
	mu = 1.905526115325545111E-01
	vecY = hc.calcYLevWBOfMu(mu)
	#vecX = hm.evalXOfY(vecY)
	vecX = hc.calcXLevWBOfMu(mu)
	f, vecG = prob.evalFG(vecX)
	msg(f'  At mu = {mu}...')
	msg(f'    vecY = {vecY}')
	msg(f'    vecX = {vecX}')
	msg(f'    f = {f}')
	msg(f'    vecG = {vecG}')
	# Req numSuperPts = 5?
	assert( abs(vecG[-1]-(-0.00519074)) < 1.0e-4 )
	msg('')
	msg(f' Along oracle_levWB...')
	#mu  = 7.967038851655461
	#mu = 1.876539040971718364E-01
	mu = 1.895312328559821957E-01
	vecY = oracle_hc.calcYLevWBOfMu(mu)
	#vecX = oracle_hm.evalXOfY(vecY)
	vecX = oracle_hc.calcXLevWBOfMu(mu)
	f, vecG = prob.evalFG(vecX)
	msg(f'  At mu = {mu}...')
	msg(f'    vecY = {vecY}')
	msg(f'    vecX = {vecX}')
	msg(f'    f = {f}')
	msg(f'    vecG = {vecG}')
	assert( abs(vecG[-1]-(-0.00549702)) < 1.0e-4 )
	msg('')
	danutil.bye()

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
	if ( t == 0.0 ):
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
	if ( (danutil.reldiff(levWB_fVals[n], levWB_fLSVals[n]) > 1.0e-6)
	  or (danutil.reldiff(levWB_vecGammaVals[:,n], levWB_vecGammaLSVals[:,n]) < 1.0e-6)
	  or (danutil.reldiff(levWB_vecGPerpVals[:,n], levWB_vecGPerpLSVals[:,n]) < 1.0e-6) ):
		msg('Warning: "LS" values did not (quite) match original.')
		msg('You may want to investigte this.')
#
oracle_muScl = min(oracle_hc.vecLambdaWB)
oracle_muVals = np.zeros(numVals)
oracle_levWB_vecYVals = np.zeros((oracle_sizeK, numVals))
oracle_levWB_dVals = np.zeros(numVals)
oracle_levWB_fVals = np.zeros(numVals)
oracle_levWB_fActVals = np.zeros(numVals)
oracle_levWB_fLSVals = np.zeros(numVals)
oracle_levWB_fPSDVals = np.zeros(numVals)
oracle_levWB_fWBVals = np.zeros(numVals)
oracle_levWB_vecGammaActVals = np.zeros((oracle_sizeK, numVals))
oracle_levWB_vecGammaVals = np.zeros((oracle_sizeK, numVals))
oracle_levWB_vecGammaLSVals = np.zeros((oracle_sizeK, numVals))
oracle_levWB_vecGammaPSDVals = np.zeros((oracle_sizeK, numVals))
oracle_levWB_vecGammaWBVals = np.zeros((oracle_sizeK, numVals))
oracle_levWB_vecGPerpActVals = np.zeros((sizeX, numVals))
oracle_levWB_vecGPerpVals = np.zeros((sizeX, numVals))
oracle_levWB_vecGPerpLSVals = np.zeros((sizeX, numVals))
oracle_levWB_vecGPerpPSDVals = np.zeros((sizeX, numVals))
oracle_levWB_vecGPerpWBVals = np.zeros((sizeX, numVals))
for n in range(numVals):
	t = tVals[n]
	if ( t == 0.0 ):
		oracle_muVals[n] = 0.0
		vecY = np.zeros(oracle_sizeK)
	else:
		mu = oracle_muScl*((1.0/t)-1.0)
		oracle_muVals[n] = mu
		vecY = oracle_hc.calcYLevWBOfMu(mu)
	oracle_levWB_vecYVals[:,n] = vecY
	oracle_levWB_dVals[n] = norm(vecY)
	oracle_levWB_fVals[n], oracle_levWB_vecGammaVals[:,n], oracle_levWB_vecGPerpVals[:,n] = oracle_hm.evalFGammaGPerpOfY(vecY)
	oracle_levWB_fLSVals[n], oracle_levWB_vecGammaLSVals[:,n], oracle_levWB_vecGPerpLSVals[:,n] = oracle_hmLS.evalFGammaGPerpOfY(vecY)
	oracle_levWB_fPSDVals[n], oracle_levWB_vecGammaPSDVals[:,n], oracle_levWB_vecGPerpPSDVals[:,n] = oracle_hmPSD.evalFGammaGPerpOfY(vecY)
	oracle_levWB_fWBVals[n], oracle_levWB_vecGammaWBVals[:,n], oracle_levWB_vecGPerpWBVals[:,n] = oracle_hmWB.evalFGammaGPerpOfY(vecY)
	vecX = oracle_hm.vecXA + (oracle_hm.matV @ vecY)
	msg(f'{n} / {numVals}')
	f, vecG = prob.evalFG(vecX)
	vecGamma = oracle_hm.matV.T @ vecG
	vecGPerp = vecG - (oracle_hm.matV @ vecGamma)
	oracle_levWB_fActVals[n] = f
	oracle_levWB_vecGammaActVals[:,n] = vecGamma[:]
	oracle_levWB_vecGPerpActVals[:,n] = vecGPerp[:]
	if ( (danutil.reldiff(oracle_levWB_fVals[n], oracle_levWB_fLSVals[n]) > 1.0e-6)
	  or (danutil.reldiff(oracle_levWB_vecGammaVals[:,n], oracle_levWB_vecGammaLSVals[:,n]) < 1.0e-6)
	  or (danutil.reldiff(oracle_levWB_vecGPerpVals[:,n], oracle_levWB_vecGPerpLSVals[:,n]) < 1.0e-6) ):
		msg('Warning: oracle "LS" values did not (quite) match original.')
		msg('You may want to investigte this.')

msgtime()
import matplotlib.pyplot as plt
plt.plot(
  levWB_dVals, levWB_fVals, 'x-',
  levWB_dVals, levWB_fPSDVals, 'o-',
  levWB_dVals, levWB_fWBVals, 's-',
  levWB_dVals, levWB_fActVals, '+-',
  oracle_levWB_dVals, oracle_levWB_fVals, 'x-',
  oracle_levWB_dVals, oracle_levWB_fPSDVals, 'o-',
  oracle_levWB_dVals, oracle_levWB_fWBVals, 's-',
  oracle_levWB_dVals, oracle_levWB_fActVals, '+-',
  levWB_dVals, levWB_fLSVals, '.-',
  oracle_levWB_dVals, oracle_levWB_fLSVals, '.-' )
plt.grid(True)
plt.legend([
  'fModel(levWB)',
  'fPSD(levWB)',
  'fWB(levWB)',
  'F(levWB)',
  'OfModel(OlevWB)',
  'OfPSD(OlevWB)',
  'OfWB(OlevWB)',
  'F(OlevWB)',
  'fLS(levWB)',
  'OfLS(OlevWB)' ])
plt.show()
