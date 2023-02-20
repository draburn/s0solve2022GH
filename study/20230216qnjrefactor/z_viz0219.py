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

dMax = 1.0
numPts = 20

pMax = dMax/norm(vecPHarvest0)
matDP = np.zeros(( sizeX, numPts ))
rvcFP = np.zeros(( 1, numPts ))
for n in range(numPts):
	temp_vecD = (n*1.0/numPts)*(pMax)*vecPHarvest0.copy()
	#temp_f = pytorchCIFAR10demo.eval_batch_loss_turbo( vecX0+temp_vecD )
	temp_f = pytorchCIFAR10demo.eval_epoch_loss(vecX0+temp_vecD)
	#temp_f = pytorchCIFAR10demo.eval_batch_loss( vecX0+temp_vecD, 0 )
	#temp_f, temp_vecG, temp_dat = pytorchCIFAR10demo.eval_epoch_sgd( vecX0+temp_vecD, vecP0, learning_rate, momentum_coefficient, 0 )
	print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e}]')
	matDP[:, n] = temp_vecD[:]
	rvcFP[0, n] = temp_f

# "W" is for "walk".
vecW = dat0.vecXHarvest - dat0.vecXSeed
print(f'||vecXHarvest0 - vecXSeed0|| = {norm(vecW)}')
wMax = dMax/norm(vecW)
matDW = np.zeros(( sizeX, numPts ))
rvcFW = np.zeros(( 1, numPts ))
for n in range(numPts):
	temp_vecD = (n*1.0/numPts)*(wMax)*vecW.copy()
	#temp_f = pytorchCIFAR10demo.eval_batch_loss_turbo( vecX0+temp_vecD )
	temp_f = pytorchCIFAR10demo.eval_epoch_loss(vecX0+temp_vecD)
	#temp_f = pytorchCIFAR10demo.eval_batch_loss( vecX0+temp_vecD, 0 )
	#temp_f, temp_vecG, temp_dat = pytorchCIFAR10demo.eval_epoch_sgd( vecX0+temp_vecD, vecP0, learning_rate, momentum_coefficient, 0 )
	print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e}]')
	matDW[:, n] = temp_vecD[:]
	rvcFW[0, n] = temp_f

gMax = dMax/norm(vecG0)
matDG = np.zeros(( sizeX, numPts ))
rvcFG = np.zeros(( 1, numPts ))
for n in range(numPts):
	temp_vecD = (n*1.0/numPts)*(-gMax)*vecG0.copy()
	#temp_f = pytorchCIFAR10demo.eval_batch_loss_turbo( vecX0+temp_vecD )
	temp_f = pytorchCIFAR10demo.eval_epoch_loss(vecX0+temp_vecD)
	#temp_f = pytorchCIFAR10demo.eval_batch_loss( vecX0+temp_vecD, 0 )
	#temp_f, temp_vecG, temp_dat = pytorchCIFAR10demo.eval_epoch_sgd( vecX0+temp_vecD, vecP0, learning_rate, momentum_coefficient, 0 )
	print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e}]')
	matDG[:, n] = temp_vecD[:]
	rvcFG[0, n] = temp_f

plt.plot(np.sqrt(np.sum(matDP*matDP,0)),rvcFP[0,:],'x-')
plt.plot(np.sqrt(np.sum(matDW*matDW,0)),rvcFW[0,:],'^-')
plt.plot(np.sqrt(np.sum(matDG*matDG,0)),rvcFG[0,:],'o-')
plt.legend(['vecPHarvest0','vecXHarvest0-vecXSeed0','vecG0'])
plt.grid(True)
plt.show()
