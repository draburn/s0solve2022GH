import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm
import demoproblem0221 as prob
#import qnj0221

vecX0 = prob.genVecX0()
sizeX = vecX0.shape[0]
vecX0 = prob.genVecX0()
vecP0 = np.zeros(sizeX)

sgdPrm = prob.evalSGD_prm()
sgdPrm.doStore = True

vecX, vecP, f, sgdDat = prob.evalSGD(vecX0, vecP0, sgdPrm)
msg(f'vecX = {vecX}')
msg(f'vecP = {vecP}')
msg(f'f = {f}')
msg(f'sgdDat = {sgdDat}...')
sgdDat.dump()

vecX, vecP, f, sgdDat = prob.evalSGD(vecX0)
msg(f'vecX = {vecX}')
msg(f'vecP = {vecP}')
msg(f'f = {f}')
msg(f'sgdDat = {sgdDat}...')
sgdDat.dump()

msg('Bye.')
exit()

sgdPrm.doStats = False
sgdPrm.doStore = True
vecX, vecP, f, sgdDat = prob.evalSGD(vecX0, vecP0, sgdPrm)
msg(f'sgdDat = {sgdDat}...')
sgdDat.dump()

sgdPrm.doStats = False
sgdPrm.doStore = True
vecX, vecP, f, sgdDat = prob.evalSGD(vecX0)
msg(f'sgdDat = {sgdDat}...')
sgdDat.dump()

msg('Bye.')
exit()

sgdPrm = prob.evalSGD_prm()
msg(f'sgdPrm = {sgdPrm}')
sgdPrm.pie = 3.14159
msg(f'sgdPrm.pie = {sgdPrm.pie}')
sgdPrm2 = prob.evalSGD_prm()
sgdPrm2.pie = 2.0000
msg(f'sgdPrm2.pie = {sgdPrm2.pie}')
msg(f'sgdPrm.pie = {sgdPrm.pie}')

vecX, vecP, f, sgdDat = prob.evalSGD(vecX0, vecP0)
msg(f'vecX = {vecX}')
msg(f'vecP = {vecP}')
msg(f'f = {f}')
sgdDat.dump()

sgdPrm = prob.evalSGD_prm()
sgdPrm.genDatOut = False
vecX, vecP, f, _ = prob.evalSGD(vecX0, vecP0, sgdPrm)
msg(f'vecX = {vecX}')
msg(f'vecP = {vecP}')
msg(f'f = {f}')
