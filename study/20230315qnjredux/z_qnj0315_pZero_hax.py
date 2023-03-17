import time
import danutil
from danutil import msg, msgtime
import numpy as np
from numpy.linalg import norm
import hessmodel0316 as hessmodel

startTime = time.time()
justTest = False
msg(f'justTest = {justTest}')
if (justTest):
	import demoproblem0221 as prob
	maxNumSteps = 7
else:
	import pytorchCIFAR10demo as prob
	maxNumSteps = 500
vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
divergenceCoeff = 100.0
sgdPrm = prob.evalSGD_prm()
sgdPrm.learningRate = 1.0e-2
sgdPrm.momentumCoefficient = 0.9
msg(f'sgdPrm = {sgdPrm}...')
sgdPrm.dump()

maxNumRecords = 100
msg(f'maxNumRecords = {maxNumRecords}')
record_matX = np.zeros((sizeX, maxNumRecords))
record_vecF = np.zeros(maxNumRecords)
record_matG = np.zeros((sizeX, maxNumRecords))
numRecords = 0

useSMOP = True
msg(f'useSMOP = {useSMOP}')

usePReset = False
if (usePReset):
	sgdPrmFOnly = sgdPrm.copy()
	sgdPrmFOnly.doStats = False
	sgdPrmFOnly.doStore = False
	sgdPrm.storageSize = 0
	msg(f'sgdPrmFOnly = {sgdPrmFOnly}...')
	sgdPrmFOnly.dump()
	def funch_evalFG( x ):
		_, _, f, _ = prob.evalSGD(x, np.zeros(sizeX), sgdPrmFOnly)
		return f, np.ones(sizeX)
else:
	def funch_evalFG( x ):
		_, _, f, d = prob.evalSGD(x, np.zeros(sizeX), sgdPrm)
		return f, d.statsDat.avg_vecG

chmPrm = hessmodel.calcHessModel_prm()
chcPrm = hessmodel.calcHessCurves_prm()
shcPrm = hessmodel.searchHessCurve_prm()
#shcPrm.tMax = 0.1
msg(f'chmPrm = {chmPrm}...')
chmPrm.dump()
msg(f'chcPrm = {chcPrm}...')
chcPrm.dump()
msg(f'shcPrm = {shcPrm}...')
shcPrm.dump()

vecXSeed = vecX0.copy()
vecPSeed = vecP0.copy()
msg(f'Starting main loop...')
print('')
print('[')
for stepIndex in range(maxNumSteps):
	# Perform feval.
	if (usePReset):
		vecPSeed[:] = 0.0
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed, sgdPrm)
	print(f'[{time.time()-startTime:8.2f}, {stepIndex:5d}, {f:10.3e}, {norm(sgdDat.statsDat.avg_vecG[:]):10.3e}]')
	assert( f < 100.0 )
	assert ( f >= 0.0 )
	#
	# Add to records.
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
	#
	# Prepare for next iteration.
	#if (useSMOP and (0==((stepIndex+1)%10))):
	if (useSMOP):
		nAnchor = 0
		vecXSeed, vecPSeed, smopDat = hessmodel.searchMin_sgd_oracleP(
		  funch_evalFG,
		  record_matX[:,nAnchor],
		  record_vecF[nAnchor],
		  record_matG[:,nAnchor],
		  record_matX[:,0:numRecords],
		  record_vecF[0:numRecords],
		  record_matG[:,0:numRecords],
		  vecXHarvest,
		  vecPHarvest,
		  chmPrm,
		  chcPrm,
		  shcPrm )
	else:
		vecXSeed[:] = vecXHarvest[:]
		vecPSeed[:] = vecPHarvest[:]
# End main loop.
print('];')
msg('Finished main loop.')


import matplotlib.pyplot as plt
hc = smopDat.hc
hm = smopDat.hm
hmLS = smopDat.hmLS
hmWB = smopDat.hmWB
numVals = 200
tExpLo = 1
tExpHi = 1
tVals = 1.0 - (1.0-(np.array(np.linspace(0.0,1.0,numVals))**tExpLo))**tExpHi
tVals /= 10.0
muScl = min(hc.vecLambdaWB)
muVals = np.zeros(numVals)
levLS_dVals  = np.zeros(numVals)
levLS_fVals  = np.zeros(numVals)
levLS_fLSVals = np.zeros(numVals)
floor_dVals  = np.zeros(numVals)
floor_fVals  = np.zeros(numVals)
floor_fWBVals = np.zeros(numVals)
for n in range(numVals):
	msg(f' {n} / {numVals}...')
	t = tVals[n]
	if ( t == 0.0 ):
		muVals[n] = 0.0
		vecY_levLS = hc.vecYLaunch.copy()
		vecY_floor = hc.vecYLaunch.copy()
	else:
		mu = muScl*((1.0/t)-1.0)
		muVals[n] = mu
		vecY_levLS = hc.calcYLevLSOfMu(mu)
		vecY_floor = hc.calcYFloorOfMu(mu)
	vecX_levLS = hm.vecXA + hm.matV @ vecY_levLS
	vecX_floor = hm.vecXA + hm.matV @ vecY_floor
	levLS_dVals[n] = norm(vecX_levLS - vecXHarvest)
	levLS_fVals[n], _ = funch_evalFG(vecX_levLS)
	levLS_fLSVals[n], _, _ = hmLS.evalFGammaGPerpOfY( vecY_levLS )
	floor_dVals[n] = norm(vecX_floor - vecXHarvest)
	floor_fVals[n], _ = funch_evalFG(vecX_floor)
	floor_fWBVals[n], _, _ = hmWB.evalFGammaGPerpOfY( vecY_floor )
	msg(f'{tVals[n]}, {muVals[n]}, {levLS_dVals[n]}, {levLS_fVals[n]}, {levLS_fLSVals[n]}')
msgtime()
smop_vecX, smop_vecP, _ = hessmodel.searchMin_sgd_oracleP(
  funch_evalFG,
  record_matX[:,nAnchor],
  record_vecF[nAnchor],
  record_matG[:,nAnchor],
  record_matX[:,0:numRecords],
  record_vecF[0:numRecords],
  record_matG[:,0:numRecords],
  vecXHarvest,
  vecPHarvest )
smop_d = norm( smop_vecX - vecXHarvest ) #Hacky
smop_f, smop_vecG = funch_evalFG(smop_vecX)
msgtime()
plt.plot( levLS_dVals, levLS_fLSVals, 'o-', markersize=20 )
plt.plot( levLS_dVals, levLS_fVals, 's-', markersize=16 )
plt.plot( floor_dVals, floor_fWBVals, 'o-', markersize=12 )
plt.plot( floor_dVals, floor_fVals, 's-', markersize=8 )
plt.plot( 0*smop_d, f+0*smop_f, '*', markersize=16 )
plt.grid(True)
plt.legend([
  'fLS (lev curve)',
  'f (lev curve)',
  'fWB (ef curve)',
  'f (ef curve)',
  'f (smop)' ])
plt.xlabel('||delta||')
plt.ylabel('F')
plt.show()
