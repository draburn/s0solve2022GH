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
msg(f'Starting solver...')
vecX0 = pytorchCIFAR10demo.get_vecX()
sizeX = vecX0.shape[0]
vecP0 = np.zeros(sizeX)
learning_rate = 0.001
momentum_coefficient = 0.9
max_num_epochs = 10
fname_x0 = ''
fname_p0 = ''
dtype_x0 = np.float32
dtype_p0 = np.float32
#fname_x0 = 'in_vecX0.np'
#fname_p0 = 'in_vecP0.np'
dump_interval = 0
msg(f'  learning_rate = {learning_rate:0.18E}')
msg(f'  momentum_coefficient = {momentum_coefficient:0.18E}')
msg(f'  max_num_epochs = {max_num_epochs}')
msg(f'  fname_x0 = "{fname_x0}"')
msg(f'  fname_p0 = "{fname_p0}"')
msg(f'  dump_interval = {dump_interval}')
if (''!=fname_x0):
	msg(f'Reading x0 from disk using dtype = "{dtype_x0}".')
	vecX0[:] = np.fromfile(fname_x0,dtype=dtype_x0)[:]
if (''!=fname_p0):
	msg(f'Reading p0 from disk using dtype = "{dtype_p0}".')
	vecP0[:] = np.fromfile(fname_p0,dtype=dtype_p0)[:]
msg(f'  vecX0: self[0] = {vecX0[0]:0.18E}, ||self|| = {np.linalg.norm(vecX0):0.18E}')
msg(f'  vecP0: self[0] = {vecP0[0]:0.18E}, ||self|| = {np.linalg.norm(vecP0):0.18E}')
msgtime()

vecX = vecX0.copy()
vecP = vecP0.copy()
msg('')
msg('Starting main loop...')
print('[')
for epoch_index in range(max_num_epochs):
	f, vecG, dat = pytorchCIFAR10demo.eval_epoch_sgd( vecX, vecP, learning_rate, momentum_coefficient, 0 )
	#dat.dump()
	print(f'[{time.time()-start_time:7.2f} {epoch_index:3d} {f:12.6e} {np.linalg.norm(vecG):12.6e}', end='')
	print(f']')
	vecX[:] = dat.vecXHarvest[:]
	vecP[:] = dat.vecPHarvest[:]
	if ( (dump_interval>0) and (epoch_index>0) and (epoch_index%dump_interval == 0) ):
		vecX.tofile(f'out_vecX_{epoch_index:06d}.np')
		vecP.tofile(f'out_vecP_{epoch_index:06d}.np')
	# End saving progress.
# End epoch loop.
print('];')
print('')
msg('Finished main loop.')
msgtime()

msg('')
msg('Goodbye.')
