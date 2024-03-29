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
import qnj
#import demolossfunc0219 as prob
import pytorchCIFAR10demo as prob
import matplotlib.pyplot as plt
msgtime()

msg('')
learning_rate = 1.0e-5
epoch_limit = 50
max_num_records = 10
divergence_coeff = 100.0
report_interval = 1
msg(f'learning_rate = {learning_rate:0.18e}')
msg('*** UGH! HACKY! SOLVER SHOULD BE SUCH-THAT learning_rate IS UNNECESSARY! ***')
msg(f'epoch_limit = {epoch_limit}')
msg(f'max_num_records = {max_num_records}')
msg(f'divergence_coeff = {divergence_coeff:0.18e}')
msg(f'report_interval = {report_interval}')
vecX0 = prob.get_vecX()
sizeX = vecX0.shape[0]
f0, vecG0 = prob.fgeval(vecX0)
msg(f'vecX0[0] = {vecX0[0]:0.18e}, ||vecX0|| = {norm(vecX0):0.18e}')
msg(f'f0 = {f0:0.18E}')
msg(f'vecG0[0] = {vecG0[0]:0.18E}, ||vecG0|| = {norm(vecG0):0.18E}')
msgtime()

msg('')
msg('Starting main loop...')
record_matX = np.zeros((sizeX, max_num_records))
record_matG = np.zeros((sizeX, max_num_records))
record_rvcF = np.zeros((1, max_num_records))
num_records = 0
vecX = vecX0
f = f0
vecG = vecG0
ignore_vecPHarvest = np.zeros(sizeX)
tr_size = 1.0e-1
qnj_prm = qnj.prm()
print('[')
print(f'[{time.time()-start_time:8.2f} {0:6d} {f0:12.6e} {norm(vecG0):12.6e}]')
for epoch_index in range(epoch_limit):
	
	# Add record.
	# DRaburn 2023-02-16: Always rolling may be wasteful, but POITROME.
	record_matX = np.roll( record_matX, 1 )
	record_matG = np.roll( record_matG, 1 )
	record_rvcF = np.roll( record_rvcF, 1 )
	if ( num_records < max_num_records ):
		num_records += 1
	record_matX[:,0] = vecX[:]
	record_matG[:,0] = vecG[:]
	record_rvcF[0,0] = f
	
	# Calculate jump.
	qnj_prm.doViz = False
	qnjCode, vecXNext, ignore_vecPNext, sizeK, gammaRat, fPred = qnj.calcJump(
	  vecX,
	  ignore_vecPHarvest,
	  record_matX[:,0:num_records],
	  record_matG[:,0:num_records],
	  record_rvcF[:,0:num_records],
	  tr_size,
	  qnj_prm )
	msg(f'{qnjCode} {sizeK} {gammaRat}')
	if ((500==qnjCode) and (sizeK>5) and (epoch_index>5)):
		qnj_prm.doViz = True
	if (qnj_prm.doViz):
		qnj_prm.viz_numPts = 100
		qnjCode, vecXNext, ignore_vecPNext, sizeK, gammaRat, fPred, lev_matDelta, lev_vecD, lev_vecF, lev_vecFC = qnj.calcJump(
		  vecX,
		  ignore_vecPHarvest,
		  record_matX[:,0:num_records],
		  record_matG[:,0:num_records],
		  record_rvcF[:,0:num_records],
		  tr_size,
		  qnj_prm )
		lev_vecFActual = np.zeros(qnj_prm.viz_numPts)
		for n in range (qnj_prm.viz_numPts):
			msg(f'Viz: {n} / {qnj_prm.viz_numPts}')
			temp_f, temp_vecG = prob.fgeval( vecX + lev_matDelta[:,n])
			lev_vecFActual[n] = temp_f
		plt.plot(lev_vecD, lev_vecF, 'o-')
		plt.plot(lev_vecD, lev_vecFC, 's-')
		plt.plot(lev_vecD, lev_vecFActual, 'x-')
		plt.grid(True)
		plt.legend(['model (orig)','model (mod)','actual'])
		plt.xlabel('||delta||')
		plt.ylabel('f')
		plt.show()
		msg('HALT!')
		exit()
	
	vecX[:] = vecXNext[:] + ((-learning_rate)*vecG[:])
	f, vecG = prob.fgeval(vecX)
	if (f > divergence_coeff*f0):
		print(f'[{time.time()-start_time:8.2f} {epoch_index:6d} {f:12.6e} {norm(vecG):12.6e}]')
		print('];')
		print('')
		msg('We appear to be diverging...')
		msg(f'  {f:0.18e} > {divergence_coeff:0.18e} * {f0:0.18e}.')
		msg('The learning_rate ({learning_rate:0.18e}) should probably be reduced.')
		break
	if ( f < 1.0e-12*f0 ):
		print('];')
		print('')
		msg(f'Hit prompt f-convergence: {f:0.18e} < 1.0e-12 * {f0:0.18e}')
		break
	elif (norm(vecG) < 1.0e-12*norm(vecG0)):
		print('];')
		print('')
		msg(f'Hit prompt gNorm-convergence: {norm(vecG):0.18e} < 1.0e-12 * {norm(vecG0):0.18e}')
		break
	if ((report_interval > 0) and ((epoch_index+1)%report_interval == 0)):
		print(f'[{time.time()-start_time:8.2f} {epoch_index+1:6d} {f:12.6e} {norm(vecG):12.6e}]')
# End epoch loop.
print('];')
print('')
print(f'[{time.time()-start_time:8.2f} {epoch_index+1:6d} {f:12.6e} {norm(vecG):12.6e}]')
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
