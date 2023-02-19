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

numPts = 1000
matD = np.zeros(( sizeX, numPts ))
rvcF = np.zeros(( 1, numPts ))
matG = np.zeros(( sizeX, numPts ))
for n in range(numPts):
	temp_vecD = 20.0*(n*1.0/numPts)*(-vecG0.copy())
	matD[:, n] = temp_vecD[:]
	temp_f = pytorchCIFAR10demo.eval_batch_loss_turbo( vecX0+temp_vecD )
	print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e}]')
	#print(f'[{time.time()-start_time:7.2f} {n:4d} {temp_f:12.6e} {np.linalg.norm(temp_vecG):12.6e}]')
	rvcF[0, n] = temp_f
	#matG[:, n] = temp_vecG[:]
import matplotlib.pyplot as plt

plt.plot(rvcF[0,:],'o-')

plt.show()
