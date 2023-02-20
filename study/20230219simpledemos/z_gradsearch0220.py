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
def ang(vecA, vecB):
	normA = norm(vecA)
	normB = norm(vecB)
	if((0.0==normA) or (0.0==normB)):
		return 0.0
	return (vecA @ vecB)/(normA*normB)
msgtime()
import demolossfunc0219 as prob
#import pytorchCIFAR10demo
import matplotlib.pyplot as plt
msgtime()

msg('')
epoch_limit = 10000
search_limit = 10
divergence_coeff = 100.0
report_interval = 1000
msg(f'epoch_limit = {epoch_limit}')
msg(f'search_limit = {search_limit}')
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
	gNorm = norm(vecG)
	assert(f > 0.0)
	assert(gNorm > 0.0)
	search_vecDelta = (-f/gNorm)*vecG
	search_sLo = 1.0e-14
	search_sHi = 2.0
	search_fLo, search_vecGLo = prob.fgeval(vecX+(search_sLo*search_vecDelta))
	search_fHi, search_vecGHi = prob.fgeval(vecX+(search_sHi*search_vecDelta))
	search_coeff = 0.1
	for search_index in range(search_limit):
		# Simple bisection.
		#temp_s = (search_sLo+search_sHi)/2.0
		temp_s = search_sLo + (search_coeff*(search_sHi-search_sLo))
		temp_f, temp_vecG = prob.fgeval(vecX+(temp_s*search_vecDelta))
		search_gNormLo = norm(search_vecGLo)
		search_gNormHi = norm(search_vecGHi)
		temp_gNorm = norm(temp_vecG)
		search_gangLo = ang(search_vecGLo, vecG)
		search_gangHi = ang(search_vecGHi, vecG)
		temp_gang = ang(temp_vecG, vecG)
		#msg(f'search: {search_sLo}, {temp_s}, {search_sHi}; {search_gangLo}, {temp_gang}, {search_gangHi}')
		#
		search_type = 2
		if (0 == search_type):
			if (temp_f >= search_fLo):
				temp_replaceHi = True
			elif (temp_f >= search_fHi):
				temp_replaceHi = False
			elif (search_fLo <= search_fHi):
				temp_replaceHi = True
			else:
				temp_replaceHi = False
		elif (1 == search_type):
			if (temp_gNorm >= search_gNormLo):
				temp_replaceHi = True
			elif (temp_gNorm >= search_gNormHi):
				temp_replaceHi = False
			elif (search_gNormLo <= search_gNormHi):
				temp_replaceHi = True
			else:
				temp_replaceHi = False
		elif (2 == search_type):
			if ( temp_gang * search_gangLo <= 0.0 ):
				temp_replaceHi = True
			elif ( temp_gang * search_gangHi <= 0.0 ):
				temp_repalceHi = False
			elif ( abs(search_gangLo) <= abs(search_gangHi) ):
				temp_replaceHi = True
			else:
				temp_replaceHi = False
		if (temp_replaceHi):
			#msg(f'Replacing Hi... search_sHi = {search_sHi}')
			search_sHi = temp_s
			search_fHi = temp_f
			search_vecGHi = temp_vecG
			#msg(f'                search_sHi = {search_sHi}')
		else:
			#msg(f'Replacing Lo... search_sLo = {search_sLo}')
			search_sLo = temp_s
			search_fLo = temp_f
			search_vecGLo = temp_vecG
			search_coeff = 0.5
			#msg(f'                search_sHi = {search_sLo}')
	# End search loop.
	doPlot137 = False
	if (doPlot137 and (137 == epoch_index)):
		msg(f'search_sLo = {search_sLo:0.18e}')
		msg(f'search_sHi = {search_sHi:0.18e}')
		numPts = 10000
		rvcS = np.zeros(numPts)
		rvcF = np.zeros(numPts)
		rvcGNorm = np.zeros(numPts)
		rvcGAng = np.zeros(numPts)
		for n in range (numPts):
			temp_s = 2.0*n/(numPts-1.0)
			temp_f, temp_vecG = prob.fgeval(vecX+(temp_s*search_vecDelta))
			rvcS[n] = temp_s
			rvcF[n] = temp_f
			rvcGNorm[n] = norm(temp_vecG)
			rvcGAng[n] = ang(temp_vecG, vecG)
		#plt.plot(rvcS,rvcF,'o-')
		plt.plot(rvcS,rvcGAng,'o-')
		plt.grid(True)
		plt.show()
		msg('HALT!')
		exit()
	#msg(f'search_sHi = {search_sHi:0.18e}')
	vecX[:] += search_sHi * search_vecDelta
	#vecX[:] += (-0.05)*vecG
	f, vecG = prob.fgeval(vecX)
	if (f > divergence_coeff*f0):
		print(f'[{time.time()-start_time:5.2f} {epoch_index:5d} {f:12.6e} {norm(vecG):12.6e}]')
		print('];')
		msg('We appear to be diverging...')
		msg(f'  {f:0.18e} > {divergence_coeff:0.18e} * {f0:0.18e}.')
		msg('The learning_rate ({learning_rate:0.18e}) should probably be reduced.')
		exit()
	if ((report_interval > 0) and ((epoch_index+1)%report_interval == 0)):
		print(f'[{time.time()-start_time:8.2f} {epoch_index+1:6d} {f:12.6e} {norm(vecG):12.6e}]')
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
