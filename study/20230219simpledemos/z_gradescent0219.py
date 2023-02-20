# DRaburn 2023-02-19:
#  This is a simple gradient-descent solver.

import inspect
import time
import numpy as np
from numpy.linalg import norm

start_time = time.time()
main_frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{main_frame.f_lineno:05d}]', *arguments, **keywords)
def msgtime():
	msg(f'It is {time.asctime()}; time since start is {time.time()-start_time:0.3f}s.')
msgtime()
import demolossfunc0219 as prob
#import pytorchCIFAR10demo
msgtime()

msg('')
learning_rate = 1.0e-4
epoch_limit = 100000
divergence_coeff = 100.0
report_interval = 10000
msg(f'learning_rate = {learning_rate:0.18e}')
msg(f'epoch_limit = {epoch_limit}')
msg(f'divergence_coeff = {divergence_coeff:0.18e}')
msg(f'report_interval = {report_interval}')
vecX0 = prob.get_vecX0()
f0, vecG0 = prob.fgeval(vecX0)
msg(f'vecX0[0] = {vecX0[0]:0.18e}, ||vecX0|| = {norm(vecX0):0.18e}')
msg(f'f0 = {f0:0.18E}')
msg(f'vecG0[0] = {vecG0[0]:0.18E}, ||vecG0|| = {norm(vecG0):0.18E}')
msgtime()

msg('')
msg('Starting main loop...')
vecX = vecX0
f = f0
vecG = vecG0
print('[')
print(f'[{time.time()-start_time:8.2f} {0:6d} {f0:12.6e} {norm(vecG0):12.6e}]')
for epoch_index in range(epoch_limit):
	vecX[:] += (-learning_rate) * vecG
	f, vecG = prob.fgeval(vecX)
	if (f > divergence_coeff*f0):
		print(f'[{time.time()-start_time:5.2f} {epoch_index:5d} {f:12.6e} {norm(vecG):12.6e}]')
		print('];')
		msg('We appear to be diverging...')
		msg(f'  {f:0.18e} > {divergence_coeff:0.18e} * {f0:0.18e}.')
		msg('The learning_rate ({learning_rate:0.18e}) should probably be reduced.')
		exit()
	if ((report_interval > 0) and ((epoch_index+1)%report_interval == 0)):
		print(f'[{time.time()-start_time:8.2f} {(epoch_index+1):6d} {f:12.6e} {norm(vecG):12.6e}]')
# End epoch loop.
print('];')
print('')
msgtime()
msg('Finished main loop.')

msg('')
vecXF = vecX.copy()
fF, vecGF = prob.fgeval(vecXF)
msg(f'vecXF[0] = {vecXF[0]:0.18e}, ||vecXF|| = {norm(vecXF):0.18e}')
msg(f'fF = {fF:0.18E}')
msg(f'vecGF[0] = {vecGF[0]:0.18E}, ||vecGF|| = {norm(vecGF):0.18E}')

msg('')
msg('Goodbye.')
msgtime()
