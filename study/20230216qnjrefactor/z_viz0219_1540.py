# DRaburn 2023-02-17:
#  This is a driver for pytorchCIFAR10 in functional form.

import inspect
import time
start_time = time.time()
main_frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{main_frame.f_lineno:05d}]', *arguments, **keywords)
def msgtime():
	msg(f'It is now {time.asctime()}; time since start is {time.time()-start_time:0.3f}s.')
msgtime()
msg('')
import numpy as np
import pytorchCIFAR10demo
import matplotlib.pyplot as plt
from numpy.linalg import norm
msgtime()

def smooth( xIn, halfKernSize ):
	xOut = xIn.copy()
	for n in range(len(xIn)):
		k = halfKernSize
		if (n-k < 0):
			k = n
		if (n+k+1 > numPts2):
			k = numPts2-n-1
		nLo = max([ 0, n-k ])
		nHi = min([ numPts2, n+k+1 ])
		xOut[n] = np.sum(xIn[nLo:nHi])/(nHi-nLo)
	return xOut
# End smooth().

msg('')
vecX0 = pytorchCIFAR10demo.get_vecX()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
learning_rate = 0.001
momentum_coefficient = 0.9
fname_x0 = ''
fname_p0 = ''
dtype_x0 = np.float32
dtype_p0 = np.float32
#fname_x0 = 'in_vecX0.np'
#fname_p0 = 'in_vecP0.np'
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msg(f'  fname_x0 = "{fname_x0}"')
msg(f'  fname_p0 = "{fname_p0}"')
if (''!=fname_x0):
	msg(f'Reading x0 from disk using dtype = "{dtype_x0}".')
	vecX0[:] = np.fromfile(fname_x0,dtype=dtype_x0)[:]
if (''!=fname_p0):
	msg(f'Reading p0 from disk using dtype = "{dtype_p0}".')
	vecP0[:] = np.fromfile(fname_p0,dtype=dtype_p0)[:]
msg(f'  vecX0: self[0] = {vecX0[0]:0.18E}, ||self|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0: self[0] = {vecP0[0]:0.18E}, ||self|| = {np.linalg.norm(vecP0):0.18E}')
msgtime()


f0, vecG0, dat0 = pytorchCIFAR10demo.eval_epoch_sgd( vecX0, vecP0, learning_rate, momentum_coefficient, 0 )
print(f'[{time.time()-start_time:7.2f} {f0:12.6e} {np.linalg.norm(vecG0):12.6e}]')
vecPHarvest0 = dat0.vecPHarvest.copy()

dMin = 0.0
dMax = 1.0
g = norm(vecG0)

numPts = 10
if (numPts > 0):
	matDG = np.zeros(( sizeX, numPts ))
	rvcFG = np.zeros(( 1, numPts ))
	rvcFGLo = np.zeros(( 1, numPts ))
	rvcFGHi = np.zeros(( 1, numPts ))
	for n in range(numPts):
		temp_vecD = -vecG0.copy()*( dMin + (dMax-dMin)*n/(numPts-1.0) )/g
		temp_f, temp_fVar = pytorchCIFAR10demo.eval_epoch_loss(vecX0+temp_vecD, 0)
		print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e} {temp_fVar:12.6e}]')
		matDG[:, n] = temp_vecD[:]
		rvcFG[0, n] = temp_f
		rvcFGLo[0, n] = temp_f - temp_fVar
		rvcFGHi[0, n] = temp_f + temp_fVar

numPts2 = 50
matDG2 = np.zeros(( sizeX, numPts2 ))
rvcFG2 = np.zeros(( 1, numPts2 ))
for n in range(numPts2):
	temp_vecD = -vecG0.copy()*( dMin + (dMax-dMin)*n/(numPts2-1.0) )/g
	temp_f = pytorchCIFAR10demo.eval_batch_loss(vecX0+temp_vecD)
	print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e}]')
	matDG2[:, n] = temp_vecD[:]
	rvcFG2[0, n] = temp_f
msgtime()

rvcFG2K1 = np.zeros(( 1, numPts2 ))
kern = 1
for n in range(numPts2):
	k = kern
	if (n-k < 0):
		k = n
	if (n+k+1 > numPts2):
		k = numPts2-n-1
	nLo = max([ 0, n-k ])
	nHi = min([ numPts2, n+k+1 ])
	rvcFG2K1[0, n] = np.sum(rvcFG2[0,nLo:nHi])/(nHi-nLo)

vecFG2S = smooth( rvcFG2[0,:], 10 )

msg(f'Generating plot...')

plt.plot(np.sqrt(np.sum(matDG2*matDG2,0)),rvcFG2[0,:],'x')
plt.plot(np.sqrt(np.sum(matDG2*matDG2,0)),vecFG2S,'.-')
if (numPts > 0):
	plt.plot(np.sqrt(np.sum(matDG*matDG,0)),rvcFG[0,:],'s-')
	plt.plot(np.sqrt(np.sum(matDG*matDG,0)),rvcFGHi[0,:],'v-')
	plt.plot(np.sqrt(np.sum(matDG*matDG,0)),rvcFGLo[0,:],'^-')

plt.grid(True)
plt.show()
