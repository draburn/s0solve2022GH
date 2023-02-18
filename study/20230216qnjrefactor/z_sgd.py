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
learning_rate = 0.01
momentum_coefficient = 0.9
msg(f'Parameters...')
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msg(f'  vecX0: self[0] = {vecX0[0]:0.18E}, ||self|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0: self[0] = {vecP0[0]:0.18E}, ||self|| = {np.linalg.norm(vecP0):0.18E}')
msgtime()

vecX = vecX0.copy()
vecP = vecP0.copy()
msg('')
msg('Starting main loop...')
print('[')
for epochIndex in range(100):
	f, vecG, dat = pytorchCIFAR10demo.eval_epoch_sgd( vecX, vecP, learning_rate, momentum_coefficient )
	print(f'[', end='')
	print(f'  ', end='')
	print(f'  {time.time()-start_time:9.2E} {epochIndex:5d} {f:15.6E}', end='')
	print(f']')
	vecX[:] = dat.vecXHarvest[:]
	vecP[:] = dat.vecPHarvest[:]
print('];')
print('')
msg('Finished main loop.')
msgtime()

msg('')
msg('Goodbye.')
