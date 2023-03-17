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

sgdPrmFOnly = sgdPrm.copy()
sgdPrmFOnly.doStats = False
sgdPrmFOnly.doStore = False
sgdPrm.storageSize = 0
chmPrm = hessmodel.calcHessModel_prm()
chcPrm = hessmodel.calcHessCurves_prm()
shcPrm = hessmodel.searchHessCurve_prm()
#shcPrm.tMax = 0.1
msg(f'sgdPrmFOnly = {sgdPrmFOnly}...')
sgdPrmFOnly.dump()
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
	vecPSeed[:] = 0.0
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed, sgdPrm)
	print(f'[{time.time()-startTime:8.2f}, {stepIndex:5d}, {f:10.3e}, {norm(sgdDat.statsDat.avg_vecG[:]):10.3e}]')
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
	if (useSMOP and (0==((stepIndex+1)%10))):
		def funch_evalFG( x ):
			_, _, f, _ = prob.evalSGD(x, np.zeros(sizeX), sgdPrmFOnly)
			return f, np.ones(sizeX)
		nAnchor = 0
		vecXSeed, vecPSeed, _ = hessmodel.searchMin_sgd_oracleP(
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
