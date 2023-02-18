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
msg(f'Input values...')
msg(f'  sizeX = {sizeX}')
msg(f'  vecX0[0] = {vecX0[0]:0.18E}')
msg(f'  ||vecX0|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0[0] = {vecP0[0]:0.18E}')
msg(f'  ||vecP0|| = {np.linalg.norm(vecP0):0.18E}')
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msgtime()

msg('')
msg(f' Input...')
msg(f'  vecX0[0] = {vecX0[0]:0.18E}')
msg(f'  ||vecX0|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0[0] = {vecP0[0]:0.18E}')
msg(f'  ||vecP0|| = {np.linalg.norm(vecP0):0.18E}')
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msg(f' Calling pytorchCIFAR10demo.eval_epoch_sgd()...')
f0, vecG0, dat0 = pytorchCIFAR10demo.eval_epoch_sgd( vecX0, vecP0, learning_rate, momentum_coefficient )
msg(f' Output...')
msg(f'  f0 = {f0:0.18E}')
msg(f'  vecG0[0] = {vecG0[0]:0.18E}')
msg(f'  ||vecG0|| = {np.linalg.norm(vecG0):0.18E}')
msg(f'  dat0 = {dat0}')
dat0.dump()
msgtime()

msg('')
msg(f' Input...')
msg(f'  vecX0[0] = {vecX0[0]:0.18E}')
msg(f'  ||vecX0|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0[0] = {vecP0[0]:0.18E}')
msg(f'  ||vecP0|| = {np.linalg.norm(vecP0):0.18E}')
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msg(f' Calling pytorchCIFAR10demo.eval_epoch_sgd()...')
f0, vecG0, dat0 = pytorchCIFAR10demo.eval_epoch_sgd( vecX0, vecP0, learning_rate, momentum_coefficient )
msg(f' Output...')
msg(f'  f0 = {f0:0.18E}')
msg(f'  vecG0[0] = {vecG0[0]:0.18E}')
msg(f'  ||vecG0|| = {np.linalg.norm(vecG0):0.18E}')
msg(f'  dat0 = {dat0}')
dat0.dump()
msgtime()

msg('')
fname_x0 = 'in_vecX0.np'
dtype_x0 = np.float32
msg(f'Reading x0 from "{fname_x0}" using dtype {dtype_x0}.')
vecX0 = np.fromfile(fname_x0,dtype=dtype_x0)
vecP0 = np.zeros(sizeX)
learning_rate = 0.0
momentum_coefficient = 0.0
msg(f'  vecX0[0] = {vecX0[0]:0.18E}')
msg(f'  ||vecX0|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0[0] = {vecP0[0]:0.18E}')
msg(f'  ||vecP0|| = {np.linalg.norm(vecP0):0.18E}')
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msgtime()

msg('')
msg(f'Calling pytorchCIFAR10demo.eval_epoch_sgd()...')
f0, vecG0, dat0 = pytorchCIFAR10demo.eval_epoch_sgd( vecX0, vecP0, learning_rate, momentum_coefficient )
msg(f'  f0 = {f0:0.18E}')
msg(f'  vecG0[0] = {vecG0[0]:0.18E}')
msg(f'  ||vecG0|| = {np.linalg.norm(vecG0):0.18E}')
msg(f'  dat0 = {dat0}')
dat0.dump()
msgtime()

msg('')
msg(f'Calling pytorchCIFAR10demo.eval_epoch_sgd()...')
f0, vecG0, dat0 = pytorchCIFAR10demo.eval_epoch_sgd( vecX0, vecP0, learning_rate, momentum_coefficient )
msg(f'  f0 = {f0:0.18E}')
msg(f'  vecG0[0] = {vecG0[0]:0.18E}')
msg(f'  ||vecG0|| = {np.linalg.norm(vecG0):0.18E}')
msg(f'  dat0 = {dat0}')
dat0.dump()
msgtime()

msg('')
msg('Goodbye!')
