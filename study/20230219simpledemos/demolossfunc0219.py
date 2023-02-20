# DRaburn 2023-02-19:
#  This is a super simple demo loss function.

import inspect
import numpy as np
main_frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{main_frame.f_lineno:05d}]', *arguments, **keywords)
# End def msg().

# Set base params; modify these as needed.
report_initialization = True
np_seed = 0
prob_size = 10

# End if (report_initialization).
np.random.seed(np_seed)
vecXC = np.random.randn(prob_size)
fC = float(abs(np.random.randn(1)))
matA = np.random.randn(prob_size,prob_size)
matH = matA.T @ matA

# Do initialization.
if (report_initialization):
	msg(f'np_seed = {np_seed}')
	msg(f'prob_size = {prob_size}')
	msg(f'fC = {fC:0.18E}')

def fgeval( vecX ):
	vecD = vecX - vecXC
	vecG = matH @ vecD
	f = fC + ((vecG @ vecD )/2.0)
	return f, vecG
# End def fgeval().

def get_vecX0():
	return np.zeros(prob_size)
# End def getVecX0().
