import time
import danutil
from danutil import msg, msgtime
import numpy as np
from numpy.linalg import norm

startTime = time.time()
justTest = False
msg(f'justTest = {justTest}')
if (justTest):
	import demoproblem0221 as prob
	maxNumSteps = 10
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
sgdPrm.doStats = True
sgdPrm.doStore = False
msg(f'sgdPrm = {sgdPrm}...')
sgdPrm.dump()

vecXSeed = vecX0.copy()
vecPSeed = vecP0.copy()
msgtime()
msg(f'Starting main loop...')
print('')
print('[')
for stepIndex in range(maxNumSteps):
	# Perform feval.
	vecXHarvest, vecPHarvest, f, sgdDat = prob.evalSGD(vecXSeed, vecPSeed, sgdPrm)
	print(f'[', end='')
	print(f'  {time.time()-startTime:9.3f},', end='')
	print(f'  {stepIndex:4d},', end='')
	print(f'  {f:12.6e}', end='')
	if (sgdPrm.doStats):
		print(f'  {norm(sgdDat.statsDat.avg_vecG[:]):12.6e}', end='')
	print(f'  ]')
	#
	# Prepare for next iteration.
	vecXSeed[:] = vecXHarvest[:]
	vecPSeed[:] = vecPHarvest[:]
# End main loop.
print('];')
msgtime()
msg'Finished main loop.')
