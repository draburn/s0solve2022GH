# DRaburn 2022-10-24
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image
# "Training a Classifier"

# Please run tut10_dl.py first.

#
# 1. Load and normalize CIFAR10
#

test_and_quit = False
print("Initializing...")
import time
start_time = time.time()

import torch
import torchvision
import torchvision.transforms as transforms

import numpy

torch.manual_seed(0)

transform = transforms.Compose(
    [transforms.ToTensor(),
     transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])

batch_size = 4

trainset = torchvision.datasets.CIFAR10(root='./data', train=True,
                                        download=False, transform=transform)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size,
                                          shuffle=True, num_workers=2)

testset = torchvision.datasets.CIFAR10(root='./data', train=False,
                                       download=False, transform=transform)
testloader = torch.utils.data.DataLoader(testset, batch_size=batch_size,
                                         shuffle=False, num_workers=2)

classes = ('plane', 'car', 'bird', 'cat',
           'deer', 'dog', 'frog', 'horse', 'ship', 'truck')

#
# 2. Define a Convolutional Neural Network
#

# Copy the neural network from the Neural Networks section
# before and modify it to take 3-channel images (instead of
# 1-channel images as it was defined).
import torch.nn as nn
import torch.nn.functional as F


class Net(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv2d(3, 6, 5)
        self.pool = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Conv2d(6, 16, 5)
        self.fc1 = nn.Linear(16 * 5 * 5, 120)
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 10)

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = torch.flatten(x, 1) # flatten all dimensions except batch
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x


net = Net()


#
# 3. Define a Loss function and optimizer
#

import torch.optim as optim

criterion = nn.CrossEntropyLoss()
optimizer = optim.SGD(net.parameters(), lr=0.001, momentum=0.9)


#
# 4. Train the network
#

def get_size_list(this_net):
   size_list = []
   size_list.append(this_net.conv1.bias.data.size())
   size_list.append(this_net.conv1.weight.data.size())
   size_list.append(this_net.conv2.bias.data.size())
   size_list.append(this_net.conv2.weight.data.size())
   size_list.append(this_net.fc1.bias.data.size())
   size_list.append(this_net.fc1.weight.data.size())
   size_list.append(this_net.fc2.bias.data.size())
   size_list.append(this_net.fc2.weight.data.size())
   size_list.append(this_net.fc3.bias.data.size())
   size_list.append(this_net.fc3.weight.data.size())
   return size_list

def get_numel_of_data(this_data):
    sz = this_data.size()
    num_dim = len(sz)
    if (1==num_dim):
        return sz[0]
    elif (2==num_dim):
        return sz[0]*sz[1]
    elif (3==num_dim):
        return sz[0]*sz[1]*sz[2]
    elif (4==num_dim):
        return sz[0]*sz[1]*sz[2]*sz[3]
    else:
        print("EXCEPTION: Unsupported number of dimensions.")
        raise BaseException

def get_numel_list(this_net):
    numel_list = [];
    numel_list.append(get_numel_of_data(this_net.conv1.bias.data))
    numel_list.append(get_numel_of_data(this_net.conv1.weight.data))
    numel_list.append(get_numel_of_data(this_net.conv2.bias.data))
    numel_list.append(get_numel_of_data(this_net.conv2.weight.data))
    numel_list.append(get_numel_of_data(this_net.fc1.bias.data))
    numel_list.append(get_numel_of_data(this_net.fc1.weight.data))
    numel_list.append(get_numel_of_data(this_net.fc2.bias.data))
    numel_list.append(get_numel_of_data(this_net.fc2.weight.data))
    numel_list.append(get_numel_of_data(this_net.fc3.bias.data))
    numel_list.append(get_numel_of_data(this_net.fc3.weight.data))
    return numel_list

def get_index_list(this_net):
    numel_list = get_numel_list(this_net)
    index_list = [0];
    n = 0
    for n in range(len(numel_list)):
        index_list.append(index_list[n] + numel_list[n])
    return index_list

def init_x_from_net(this_net):
    index_list = get_index_list(this_net)
    this_x = numpy.zeros(index_list[-1], dtype=numpy.float32)
    this_x[index_list[0]:index_list[1]] = numpy.reshape(this_net.conv1.bias.data.numpy(), -1)
    this_x[index_list[1]:index_list[2]] = numpy.reshape(this_net.conv1.weight.data.numpy(), -1)
    this_x[index_list[2]:index_list[3]] = numpy.reshape(this_net.conv2.bias.data.numpy(), -1)
    this_x[index_list[3]:index_list[4]] = numpy.reshape(this_net.conv2.weight.data.numpy(), -1)
    this_x[index_list[4]:index_list[5]] = numpy.reshape(this_net.fc1.bias.data.numpy(), -1)
    this_x[index_list[5]:index_list[6]] = numpy.reshape(this_net.fc1.weight.data.numpy(), -1)
    this_x[index_list[6]:index_list[7]] = numpy.reshape(this_net.fc2.bias.data.numpy(), -1)
    this_x[index_list[7]:index_list[8]] = numpy.reshape(this_net.fc2.weight.data.numpy(), -1)
    this_x[index_list[8]:index_list[9]] = numpy.reshape(this_net.fc3.bias.data.numpy(), -1)
    this_x[index_list[9]:index_list[10]] = numpy.reshape(this_net.fc3.weight.data.numpy(), -1)
    return this_x

def init_grad_from_net(this_net):
    index_list = get_index_list(this_net)
    this_grad = numpy.zeros(index_list[-1], dtype=numpy.float32)
    this_grad[index_list[0]:index_list[1]] = numpy.reshape(this_net.conv1.bias.data.numpy(), -1)
    this_grad[index_list[1]:index_list[2]] = numpy.reshape(this_net.conv1.weight.data.numpy(), -1)
    this_grad[index_list[2]:index_list[3]] = numpy.reshape(this_net.conv2.bias.data.numpy(), -1)
    this_grad[index_list[3]:index_list[4]] = numpy.reshape(this_net.conv2.weight.data.numpy(), -1)
    this_grad[index_list[4]:index_list[5]] = numpy.reshape(this_net.fc1.bias.data.numpy(), -1)
    this_grad[index_list[5]:index_list[6]] = numpy.reshape(this_net.fc1.weight.data.numpy(), -1)
    this_grad[index_list[6]:index_list[7]] = numpy.reshape(this_net.fc2.bias.data.numpy(), -1)
    this_grad[index_list[7]:index_list[8]] = numpy.reshape(this_net.fc2.weight.data.numpy(), -1)
    this_grad[index_list[8]:index_list[9]] = numpy.reshape(this_net.fc3.bias.data.numpy(), -1)
    this_grad[index_list[9]:index_list[10]] = numpy.reshape(this_net.fc3.weight.data.numpy(), -1)
    return this_grad

size_list = get_size_list(net)
index_list = get_index_list(net)
sxsolve_x = init_x_from_net(net)

# The following appears to make data point to the same memory as sxsolve_x.
# Documentation seems to suggest this is guaranteed for torch.from_numpy(),
#  though I do not see a guarantee for numpy.reshape().
# Not sure what happens to the memory originally allocated for data, and don't care.
net.conv1.bias.data   = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[0]:index_list[1]],size_list[0]))
net.conv1.weight.data = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[1]:index_list[2]],size_list[1]))
net.conv2.bias.data   = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[2]:index_list[3]],size_list[2]))
net.conv2.weight.data = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[3]:index_list[4]],size_list[3]))
net.fc1.bias.data   = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[4]:index_list[5]],size_list[4]))
net.fc1.weight.data = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[5]:index_list[6]],size_list[5]))
net.fc2.bias.data   = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[6]:index_list[7]],size_list[6]))
net.fc2.weight.data = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[7]:index_list[8]],size_list[7]))
net.fc3.bias.data   = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[8]:index_list[9]],size_list[8]))
net.fc3.weight.data = torch.from_numpy(numpy.reshape(sxsolve_x[index_list[9]:index_list[10]],size_list[9]))

#optimizer = optim.SGD(net.parameters(), lr=0.001, momentum=0.9)
sxsolve_lr = 0.001
sxsolve_momentum = 0.9
print(f"sxsolve_lr = {sxsolve_lr}")
print(f"sxsolve_momentum = {sxsolve_momentum}")
do_grad_init = True
sxsolve_step = numpy.zeros(index_list[-1],dtype=numpy.float32)

print("Initialization complete.")
print(f"  Elapsed time = {time.time()-start_time}s")
print(f"test_and_quit = {test_and_quit}")


# BEGIN QNJ INITIALIZATION.

import inspect
import numpy as np
from numpy.random import default_rng
from scipy import linalg
from scipy import optimize

# Init logging.
frame = inspect.currentframe()
def msg( *arguments, **keywords ):
    #print( f"[", __file__, ".", frame.f_lineno, "] ", *arguments, **keywords )
    print( f'[{__file__}.{frame.f_lineno:05d}]', *arguments, **keywords )
def reldiff( a, b ):
    sa = np.sum(np.abs(a))
    sb = np.sum(np.abs(b))
    if ( 0.0 == sa and 0.0 == sb ):
        return 0.0
    return ( np.sum(np.abs(a-b)) / ( sa + sb ) )
def utorthdrop( matA, dropRelThresh, dropAbsThresh ):
    #msg( 'Hey hey hey!' )
    matV = matA.copy() # Superfluous?
    sizeK = matA.shape[1]
    #msg( 'sizeK = ', sizeK )
    rvcDrop = np.zeros( (sizeK), dtype='bool' )
    rvcDrop[:] = False # Superfluous?
    #msg( 'matV = \n', matV )
    #msg( 'rvcDrop =', rvcDrop )
    for k in range ( 0, sizeK ):
        vNorm = linalg.norm( matV[:,k] )
        if ( vNorm <= dropAbsThresh ):
            rvcDrop[k] = True
            matV[:,k] = 0.0
        else:
            #msg( f'divding {k} by {vNorm}.' )
            matV[:,k] /= vNorm
    #msg( 'matV = \n', matV )
    #msg( 'rvcDrop =', rvcDrop )
    for k in range ( 1, sizeK ):
        ###if ( rvcDrop[k] ):
        ###    continue
        #msg( 'k = ', k )
        #msg( 'matV[:,k] = ', matV[:,k] )
        matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
        #msg( 'matV[:,k] = ', matV[:,k] )
        vNorm = linalg.norm( matV[:,k] )
        if ( vNorm <= dropRelThresh ):
            rvcDrop[k] = True
            matV[:,k] = 0.0
        else:
            matV[:,k] /= vNorm
            matV[:,k] -= matV[:,0:k] @ ( matV[:,0:k].T @ matV[:,k] )
            vNorm = linalg.norm( matV[:,k] )
            if ( vNorm <= dropRelThresh ):
                rvcDrop[k] = True
                matV[:,k] = 0.0
            else:
                matV[:,k] /= vNorm
                # Note: if dropThresh is too small, may end up keeping more than sizeX vectors.
    #msg( 'matV = \n', matV )
    #msg( 'rvcDrop =', rvcDrop )
    rvcKeep = ~rvcDrop
    matV = matV[:,rvcKeep]
    #msg( 'matV = \n', matV )
    return ( matV, rvcKeep )
def getLambdaFloor( f0, vecPhi, vecLambda, fFloor ):
    def fRes( lambdaF ):
        fRes = f0 - fFloor
        for k in range( 0, vecLambda.shape[0] ):
            fRes -= ((vecPhi[k])**2) / (2.0*max( vecLambda[k], lambdaF ))
        return fRes
    assert f0 > fFloor
    if ( min(vecLambda) > 0.0 ):
        if ( fRes( 0.0 ) >= 0.0 ):
            return 0.0
    lambdaHi = (vecPhi @ vecPhi) / ( 2.0 * ( f0 - fFloor )  )
    #msg( 'lambda: ', lambdaHi )
    if ( lambdaHi > max(vecLambda) ):
        return lambdaHi
    lambdaLo = MYEPS * max(vecLambda)
    #msg( 'lambda: ', lambdaHi, lambdaLo )
    if ( fRes( lambdaLo ) >= 0.0 ):
        return lambdaLo
    lambdaF = optimize.bisect( fRes, lambdaLo, lambdaHi )
    #msg( 'lambda: ', lambdaHi, lambdaLo, lambdaF )
    #msg( 'res = ', fRes(lambdaF) )
    #exit()
    return lambdaF
# Note:
#  "zeta" here was "phi" in 20230106stage\levsol0111.m;
#  "phi" here was "gamma" in 20230106stage\levsol0111.m.
def levsol( f0, vecPhi, matPsi, vecLambdaCurve, sMax, vecS, dMax, vecLambdaObjf, fMin ):
    def zetaOfP( p ):
        if ( p < MYEPS ):
            vecZeta = p * vecPhi / np.min( vecLambdaCurve )
        else:
            mu = np.min( vecLambdaCurve ) * ( (1.0/p) - 1.0 )
            vecZeta = vecPhi / ( vecLambdaCurve + mu )
        return vecZeta
    def deltaYOfP( p ):
        vecZeta = zetaOfP( p )
        vecDeltaY = ( matPsi @ vecZeta ) / vecS
        return vecDeltaY
    def fOfP( p ):
        vecZeta = zetaOfP( p )
        f = f0 - ( vecZeta @ vecPhi ) + (( vecZeta @ ( vecLambdaObjf * vecZeta ))/2.0)
        return f
    def sPastMaxOfP( p ):
        return linalg.norm( zetaOfP(p) ) - sMax
    def dPastMaxOfP( p ):
        return linalg.norm( deltaYOfP(p) ) - dMax
    def fTillMinOfP( p ):
        return fOfP( p ) - fMin
    assert np.min(vecS) > 0.0
    assert np.min(vecLambdaCurve) > 0.0
    assert linalg.norm(vecPhi) > 0.0
    p1 = 1.0
    if ( sMax > 0.0 ):
        assert sPastMaxOfP(0.0) < 0.0
        if ( sPastMaxOfP(p1) > 0.0 ):
            p1New = optimize.bisect( sPastMaxOfP, 0.0, p1 )
            p1 = p1New
    elif ( 0.0 == sMax ):
        p1 = 0.0
    if ( dMax > 0.0 ):
        assert dPastMaxOfP(0.0) < 0.0
        if ( dPastMaxOfP(p1) > 0.0 ):
            p1New = optimize.bisect( dPastMaxOfP, 0.0, p1 )
            p1 = p1New
    elif ( 0.0 == dMax ):
        p1 = 0.0
    # Ugh. Just require fMin.
    assert fTillMinOfP(0.0) > 0.0
    if ( fTillMinOfP(p1) < 0.0 ):
        p1New = optimize.bisect( fTillMinOfP, 0.0, p1 )
        p1New = p1
    return deltaYOfP( p1 )
# Init problem.
MYEPS = 1.0E-8
sizeX = index_list[-1]
# DRaburn 2023-01-27, pytorchDanmo: Not ANY of you.
#rngSeed = 0
#sizeX = 50
#sizeF = sizeX
#msg( 'rngSeed = ', rngSeed )
#msg( 'sizeX = ', sizeX )
#msg( 'sizeF = ', sizeF )
#rng = default_rng(rngSeed)
#matA = rng.standard_normal(( sizeF, sizeX ))
#matA = np.diag(np.linspace(1.0,sizeX,sizeX))
#matHCrit = np.zeros(( sizeX, sizeX ))
#matHCrit[:,:] = matA.T @ matA
#vecXCrit = rng.standard_normal(( sizeX ))
#vecXCrit = np.ones(( sizeX ))
#fCrit = 10.0
#noiseX = 1.0E-5
#noiseX = 0.0
#def funcFG( x ):
#    #d = x - vecXCrit
#    d = x - vecXCrit + noiseX*rng.standard_normal(( sizeX ))
#    g = matHCrit @ d
#    f = fCrit + (( d @ g )/2.0)
#    return ( f, g )
vecX0 = sxsolve_x.copy()
#f0, vecG0 = funcFG( vecX0 )
#msg( 'f0 = ', f0 )
#msg( '||vecG0|| = ', linalg.norm(vecG0) )

# Init SGD solver.
# DRaburn 2023-01-27, pytorchDanmo: I failed to write the tolerances integrably.
#fBail = f0 * 1E8
fBail = 1.0E8
fevalLimit = 1000000
learningRate = sxsolve_lr
momentumFactor = sxsolve_momentum
msg( 'fevalLimit = ', fevalLimit )
msg( 'learningRate = ', learningRate )
msg( 'momentumFactor = ', momentumFactor )
fevalCount = 0
# DRaburn 2023-01-27, pytorchDanmo: Hook-up to pre-existing code.
vecX = sxsolve_x
vecP = sxsolve_step

# Init superPt.
numFevalPerSuperPt = 1000
superPtLimit = 1000
# DRaburn 2023-01-27, pytorchDanmo: I failed to write the tolerances integrably.
#fTol = f0*1.0E-12
#gTol = linalg.norm(vecG0)*1.0E-6
fTol = 1.0E-4
gTol = 1.0E-12
xTol = sizeX * 1.0E-12
#fTol = 1.0E-6
#gTol = 1.0E-6
msg( 'numFevalPerSuperPt = ', numFevalPerSuperPt)
msg( 'superPtLimit = ', superPtLimit )
msg( 'fTol = ', fTol)
msg( 'gTol = ', gTol )
running_fevalCount = 0
#running_fMin = -1.0
running_fTot = 0.0
running_fSqTot = 0.0
running_xtgTot = 0.0
running_vecGTot = np.zeros(( sizeX ))
running_vecXTot = np.zeros(( sizeX ))
running_fSqTot = 0.0
running_xtgSqTot = 0.0
running_vecGSqTot = np.zeros(( sizeX ))
running_vecXSqTot = np.zeros(( sizeX ))
superPtCount = 0
vecXSeed = vecX.copy()
vecPSeed = vecP.copy()
vecXHarvest = np.zeros(( sizeX ))
vecPHarvest = np.zeros(( sizeX ))
superPt_f = 0.0
superPt_fVar = 0.0
superPt_vecG = np.zeros(( sizeX ))
superPt_vecX = np.zeros(( sizeX ))

# Init minf and best...
#  "best" is superPt with min ||vecG|| sbjt f not too much larger than minf_f.
coeff_best_minf = 1.0
coeff_best_best = 1.0
coeff_best_curr = 1.0
forceNewAsBest = True
msg( 'coeff_best_minf = ', coeff_best_minf )
msg( 'coeff_best_best = ', coeff_best_best )
msg( 'coeff_best_curr = ', coeff_best_curr )
if (forceNewAsBest):
    msg( '*** WARNING: forceNewAsBest = ', forceNewAsBest )
else:
    msg( 'forceNewAsBest = ', forceNewAsBest )
minf_present = False
minf_f = 0.0
minf_fVar = 0.0
minf_vecG = np.zeros(( sizeX ))
minf_vecX = np.zeros(( sizeX ))
best_present = False
best_f = 0.0
best_fVar = 0.0
best_vecG = np.zeros(( sizeX ))
best_vecX = np.zeros(( sizeX ))
best_vecXHarvest = np.zeros(( sizeX ))
best_vecPHarvest = np.zeros(( sizeX ))
badCount = 0

# Init records.
maxNumRecords = 20
msg( 'maxNumRecords = ', maxNumRecords )
record_matX = np.zeros(( sizeX, maxNumRecords ))
record_matG = np.zeros(( sizeX, maxNumRecords ))
record_rvcF = np.zeros(( 1, maxNumRecords ))
numRecords = 0

# Init QNJ.
useQNJ = True # Unless...
useQNJ = False
maxSubspaceSize = maxNumRecords
qnj_dropThresh = 0.1
msg( 'useQNJ = ', useQNJ )
msg( 'maxSubspaceSize = ', maxSubspaceSize )
msg( 'qnj_dropThresh = ', qnj_dropThresh )
# Pre-alloc workspaces.
###matD = np.zeros(( sizeX, maxNumRecords ))
###matV = np.zeros(( sizeX, maxSubspaceSize ))
sizeK = 0
qnj_havePrev = False
qnj_sPrev = -1.0
qnj_sMax = 1.0
qnj_sMax_btCoeff = 0.1
qnj_sMax_ftCoeff = 1.2
qnj_dPrev = -1.0
qnj_dMax = -1.0
qnj_dMax_btCoeff = 0.1
qnj_dMax_ftCoeff = 1.2

doMainLoop = True

# END QNJ INITIALIZATION



print("Main loop...")
for epoch in range(2):  # loop over the dataset multiple times
    
    running_loss = 0.0
    running_feval_count = 0
    running_time0 = time.time()
    for i, data in enumerate(trainloader, 0):
        
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        
        #
        sxsolve_f = loss
        if (do_grad_init):
            sxsolve_gradloc = init_grad_from_net(net)
            net.conv1.bias.grad   = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[0]:index_list[1]],size_list[0]))
            net.conv1.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[1]:index_list[2]],size_list[1]))
            net.conv2.bias.grad   = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[2]:index_list[3]],size_list[2]))
            net.conv2.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[3]:index_list[4]],size_list[3]))
            net.fc1.bias.grad   = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[4]:index_list[5]],size_list[4]))
            net.fc1.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[5]:index_list[6]],size_list[5]))
            net.fc2.bias.grad   = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[6]:index_list[7]],size_list[6]))
            net.fc2.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[7]:index_list[8]],size_list[7]))
            net.fc3.bias.grad   = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[8]:index_list[9]],size_list[8]))
            net.fc3.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_gradloc[index_list[9]:index_list[10]],size_list[9]))            
            sxsolve_f_min = sxsolve_f
            do_grad_init = False
        
        # DRaburn 2023-01-27, pytorchDanmo: Interface
        fevalCount += 1
        #vecX = sxsolve_x.copy()
        vecG = sxsolve_gradloc.copy()
        f = float( sxsolve_f.detach().numpy() )
        #msg( type( vecX ) )
        #msg( type( vecG ) )
        #msg( type( f ) )
        #exit()
        
        # DRaburn 2023-01-27, pytorchDanmo: QNJ
        
        # Update superPt running totals before updating SGD.
        xtg = vecX @ vecG
        running_fevalCount += 1
        running_fTot += f
        running_fSqTot += f*f
        running_xtgTot += xtg
        running_vecGTot[:] += vecG[:]
        running_vecXTot[:] += vecX[:]
        running_xtgSqTot += xtg*xtg
        running_vecGSqTot[:] += vecG[:]**2
        running_vecXSqTot[:] += vecX[:]**2
        #if ( running_fMin < 0.0 or f < running_fMin ):
        #    running_fMin = f
        
        # Update SGD.
        vecP[:] = ( momentumFactor * vecP[:] ) - ( learningRate * vecG[:] )
        vecX[:] += vecP[:]
        #print( 'vecX =\n', vecX )
        
        # Check per-feval stop crit.
        if ( f > fBail ):
            msg( 'IMPOSED STOP: f > fBail. This strongly indicates divergence.' )
            doMainLoop = False
        elif ( fevalLimit > 0 and fevalCount >= fevalLimit ):
            msg( 'IMPOSED STOP: fevalCount >= fevalLimit.' )
            doMainLoop = False
        # Check elapsed time?.
        # Check for "stop signal on disk"?
        if ( not doMainLoop ):
            break
        
        # Have we finished the super point?
        if ( running_fevalCount < numFevalPerSuperPt ):
            continue
        vecXHarvest[:] = vecX[:]
        vecPHarvest[:] = vecP[:]
        
        # Do super-point analysis.
        superPtCount += 1
        #
        superPt_vecG[:] = running_vecGTot[:] / running_fevalCount
        superPt_vecX[:] = running_vecXTot[:] / running_fevalCount
        superPt_fAvg = running_fTot / running_fevalCount
        superPt_xtgAvg = running_xtgTot / running_fevalCount
        ###superPt_f = superPt_fAvg - (( superPt_xtgAvg - ( superPt_vecX @ superPt_vecG ) )/2.0)
        superPt_f = superPt_fAvg
        assert superPt_f >= -fTol
        superPt_fSqVar = (running_fSqTot/running_fevalCount) - (superPt_fAvg**2)
        if ( 0.0 < superPt_fSqVar ):
            superPt_fVar = np.sqrt( superPt_fSqVar )
        else:
            superPt_fVar = 0.0
        # Note that part of fVar is due to f actually varying along the path;
        #  this is not desirable, but is probably acceptable.
        superPt_vecXSqVar = (running_vecXSqTot[:] / running_fevalCount) - superPt_vecX[:]**2
        superPt_vecGSqVar = (running_vecGSqTot[:] / running_fevalCount) - superPt_vecG[:]**2
        superPt_gVar = np.sqrt( np.sum( np.abs(superPt_vecGSqVar) ) )
        
        #
        running_fevalCount = 0
        running_fTot = 0.0
        running_xtgTot = 0.0
        running_vecGTot[:] = 0.0
        running_vecXTot[:] = 0.0
        running_xtgSqTot = 0.0
        running_vecGSqTot[:] = 0.0
        running_vecXSqTot[:] = 0.0
        
        # Do minf & best analysis.
        newIsMinf = False # Unless...
        newIsBest = False # Unless...
        if ( forceNewAsBest ):
            newIsMinf = True
            newIsBest = True
        elif ( not minf_present ):
            newIsMinf = True
            newIsBest = True
        elif ( superPt_f < minf_f ):
            newIsMinf = True
            newIsBest = True
        else:
            fBestThresh = (  minf_f
              + ( coeff_best_minf * minf_fVar )
              + ( coeff_best_best * best_fVar )
              + ( coeff_best_curr * superPt_fVar )  )
            if ( superPt_f <= fBestThresh and linalg.norm(superPt_vecG) < linalg.norm(best_vecG) ):
                newIsBest = True
        if ( newIsMinf ):
            minf_f = superPt_f
            minf_fVar = superPt_fVar
            minf_vecX[:] = superPt_vecX[:]
            minf_vecG[:] = superPt_vecG[:]
            minf_present = True
        if ( newIsBest ):
            best_f = superPt_f
            best_fVar = superPt_fVar
            best_vecX[:] = superPt_vecX[:]
            best_vecG[:] = superPt_vecG[:]
            best_vecXHarvest[:] = vecXHarvest[:]
            best_vecPHarvest[:] = vecPHarvest[:]
            best_present = True
        else:
            badCount += 1
        
        # Print progress log.
        #if ( 0 == superPtCount % 10 ):
        if ( True ):
            if ( newIsMinf ):
                progLogSymbol = '*'
            elif ( newIsBest ):
                progLogSymbol = '.'
            else:
                progLogSymbol = 'X'
            msg(
              f' {superPtCount:4d} ({badCount:4d}X), {fevalCount:7d}:',
              f' {sizeK:3d} / {numRecords:3d}:'
              f'  {linalg.norm( best_vecX - vecX0 ):8.2E};',
              f'  {linalg.norm( vecXHarvest - vecXSeed ):8.2E}, {qnj_dPrev:9.2E} / {qnj_dMax:9.2E};',
              f'  {best_f:8.2E};',
              f'  {linalg.norm(best_vecG):8.2E} +/- {superPt_gVar:8.2E}',
              progLogSymbol )
        
        # Check superPt stop crit.
        if ( linalg.norm(superPt_vecG) <= gTol ):
            msg( 'SUCCESS: linalg.norm(superPt_vecG) <= gTol.' )
            doMainLoop = False
        elif ( superPt_f <= fTol ):
            msg( 'SUCCESS: superPt_f <= fTol.' )
            doMainLoop = False
        elif ( superPtCount > 0 and superPtCount >= superPtLimit ):
            msg( 'IMPOSED STOP: superPtCount >= superPtLimit.' )
            doMainLoop = False
        elif ( qnj_havePrev and ( qnj_dPrev < xTol ) and (not newIsBest) ):
            msg( 'IMPOSED STOP: Failed to improve with a QNJ step smaller than xTol.' )
            doMainLoop = False
        # Check superPt_vecX vs prev?
        # Check vecXHarvest vs vecXSeed?
        # Check superPt_f vs prev?
        if ( not doMainLoop ):
            break
        
        # Prepare for next iteration.
        # Record seed for posterity.
        # This will almost certainly be modified if use a quasi-newton jump.
        sizeK = 0
        
        forceBasisGen = True # For comparison to Octave code.
        if ( (not useQNJ) and (not forceBasisGen) ):
            vecXSeed[:] = vecX[:]
            vecPSeed[:] = vecP[:]
            sxsolve_x[:] = vecX[:]
            continue
        
        # Add information to records.
        # Always rolling is wasteful. POITROME.
        record_matX = np.roll( record_matX, 1 )
        record_matG = np.roll( record_matG, 1 )
        record_rvcF = np.roll( record_rvcF, 1 )
        # Does this not require unnecessary mem alloc and copy?
        if ( numRecords < maxNumRecords ):
            numRecords += 1
        record_matX[:,0] = superPt_vecX[:]
        record_matG[:,0] = superPt_vecG[:]
        record_rvcF[0,0] = superPt_f
        
        #msg( 'record_matX =\n', record_matX )
        
        # Finally, QNJ!... or not.
        if ( numRecords < 2 ):
            vecXSeed[:] = vecX[:]
            vecPSeed[:] = vecP[:]
            sxsolve_x[:] = vecX[:]
            continue
        elif ( not newIsBest ):
            vecX[:] = best_vecXHarvest[:]
            vecP[:] = best_vecPHarvest[:]
            vecXSeed[:] = vecX[:]
            vecPSeed[:] = vecP[:]
            sxsolve_x[:] = vecX[:]
            continue
            # This is 'grad-if-bad'.
        
        # Generate basis.
        # DRaburn 2023-01-24: This is a crude two-pass QR method,
        #  involving an unfortunate cap to the number of records before the QR calculation.
        # POITROME
        vecXAnchor = best_vecX # Shallow copy / reference only / DO NOT MODIFY!
        vecGAnchor = best_vecG # Shallow copy / reference only / DO NOT MODIFY!
        fAnchor = best_f
        ###matD[:,0:maxNumRecords] = record_matX[:,0:maxNumRecords] - np.reshape( vecXAnchor, (sizeX,1) ) # Autobroadcast.
        #msg( 'numRecords = ', numRecords )
        matD = record_matX[:,0:numRecords].copy() - np.reshape( vecXAnchor, (sizeX,1) ) # Autobroadcast.
        matG = record_matG[:,0:numRecords].copy()
        #msg( 'vecXAnchor = ', vecXAnchor )
        #msg( 'vecGAnchor = ', vecGAnchor )
        #msg( 'fAnchor = ', fAnchor )
        #msg( 'matD =\n', matD )
        #msg( 'matG =\n', matG )
        # We want an equivalent of my Octave "utorthdrop":
        #  construct a basis upper-triangularly, dropping any vectors that are below some threshold in orthogonality.
        matQ, rvcKeep = utorthdrop( matD, qnj_dropThresh, 1.0E-16 )
        matD = matD[:,rvcKeep]
        matG = matG[:,rvcKeep]
        matR = np.triu( matQ.T @ matD )
        sizeK = matQ.shape[1]
        #msg( 'matQ.T @ matQ =\n', matQ.T @ matQ )
        #msg( 'matQ =\n', matQ )
        #msg( 'matR =\n', matR )
        #msg( 'sizeK = ', sizeK )
        #msg( 'D =\n', matD )
        #msg( 'Q*R =\n', matQ @ matR )
        if ( 0 == sizeK ):
            vecXSeed[:] = vecX[:]
            vecPSeed[:] = vecP[:]
            sxsolve_x[:] = vecX[:]
            continue
        matGamma = matQ.T @ matG
        vecGammaAnchor = matQ.T @ vecGAnchor
        if ( not useQNJ ):
            vecXSeed[:] = vecX[:]
            vecPSeed[:] = vecP[:]
            sxsolve_x[:] = vecX[:]
            continue
        
        # Generate fit.
        # 2023-02-24: This is simplisic but reasonable.
        #  However, see "hessfit.m".
        vecGammaFit = vecGammaAnchor
        fFit = fAnchor
        matA = linalg.solve( matR.T, (matGamma - np.reshape( vecGammaAnchor, (sizeK,1) ) ).T )
        matHFit = ( matA.T + matA )/2.0
        #
        #vecGammaTrue = matQ.T @ matHCrit @ ( vecXAnchor - vecXCrit )
        #matHTrue = matQ.T @ matHCrit @ matQ
        #msg( 'vecGammaFit = ', vecGammaFit )
        #msg( 'vecGammaTrue = ', vecGammaTrue )
        #msg( 'matHFit =\n', matHFit )
        #msg( 'matHTrue = \n', matHTrue )
        #if ( np.sum(np.abs(matHFit-matHTrue)) >= 1.0e-8*( np.sum(np.abs(matHFit)) + np.sum(np.abs(matHTrue)) ) ):
        #    msg( 'ERROR: fit is not true.' )
        #    doMainLoop = False
        #    break
        
        #vecDelta = matQ @ linalg.solve( matHFit, -vecGammaFit )
        #vecX[:] = vecXAnchor + vecDelta
        #vecXSeed[:] = vecX
        #vecPSeed[:] = vecP
        #continue
        
        # Update trust region and scaling.

        if ( qnj_havePrev ):
            if ( newIsBest ):
                qnj_sMax = qnj_sPrev * qnj_sMax_ftCoeff
                qnj_dMax = qnj_dPrev * qnj_dMax_ftCoeff
            else:
                qnj_sMax = qnj_sPrev * qnj_sMax_btCoeff
                qnj_dMax = qnj_dPrev * qnj_dMax_btCoeff
        else:
            if ( qnj_dMax < linalg.norm( vecXHarvest - vecXSeed ) ):
                qnj_dMax = linalg.norm( vecXHarvest - vecXSeed )
        #msg( 'caps: ', qnj_sMax, qnj_dMax )
        #msg( 'matR =\n', matR )
        vecCap = np.max( np.abs(matR), 1 )
        vecCap[:] += np.sqrt(MYEPS)*np.max(vecCap)
        #msg( 'vecCap =', vecCap )
        vecS = 1.0 / vecCap.copy()
        #msg( 'vecS = ', vecS )
        matS = np.diag(vecS)
        matSInv = np.diag(1.0/vecS)
        vecGammaScl = matSInv @ vecGammaFit
        matHScl = matSInv @ matHFit @ matSInv
        
        # Apply scaling and do eigenfactorization
        vecLambdaC, matPsi = linalg.eig( matHScl )
        for n in range( 0, vecLambdaC.shape[0]):
            assert np.isreal(vecLambdaC[n])
        vecLambdaOrig = np.real( vecLambdaC )
        vecPhi = matPsi.T @ (-vecGammaScl)
        # So, now:
        #  matM = matLambda + mu * matI
        #  vecZ = matSInv @ matPsi @ ( matM \ vecPhi )
        #  s = np.norm( matS * vecZ ) = np.norm( matM \ vecPhi )
        #  vecDelta = matV @ vecZ
        #  d = np.norm( vecDelta ) = np.norm( vecZ )
        #msg( "But... we want to do all of this from LAUNCH not ANCHOR?" )
        
        # Calculate lambdaMod,
        #  lambda perturbed so that Hessian is pos-def and fModMin >= 0.0
        # Note: we might consider doing this from our launch rather than anchor. Oh well.
        #msg( 'vecLambdaOrig = ', vecLambdaOrig )
        doLambdaFloorTest = False
        if (doLambdaFloorTest):
            fFit = 1.0
            vecPhi = np.array([1.0,1.0,1.0])
            #vecLambdaOrig = np.array([10.0,0.01,-1.0])
            vecLambdaOrig = np.array([0.0,0.0,-1.0])
            vecLambdaOrig = np.array([10.0,10.0,10.0])
            vecLambdaOrig = np.array([10.0,10.0,0.0])
            vecLambdaOrig = np.array([10.0,10.0,0.56])
            vecLambdaOrig = np.array([10.0,10.0,0.55])
            lambdaFloor = getLambdaFloor( fFit, vecPhi, vecLambdaOrig, 0.0 )
            msg( 'lambdaFloor =', lambdaFloor )
            vecLambdaMod = vecLambdaOrig.copy()
            for k in range ( 0, vecLambdaMod.shape[0] ):
                if ( vecLambdaMod[k] < lambdaFloor ):
                    vecLambdaMod[k] = lambdaFloor
            fRes = fFit - (np.sum( vecPhi * vecPhi / vecLambdaMod )/2.0)
            msg( 'fRes = ', fRes )
            exit()
        lambdaFloor = getLambdaFloor( fFit, vecPhi, vecLambdaOrig, -0.01*fFit )
        vecLambdaMod = vecLambdaOrig.copy()
        for k in range ( 0, vecLambdaMod.shape[0] ):
            if ( vecLambdaMod[k] < lambdaFloor ):
                vecLambdaMod[k] = lambdaFloor
        matHMod = matS @ matPsi @ np.diag(vecLambdaMod) @ (matPsi.T) @ matS
        
        # Decompose "launch".
        vecXLaunch = best_vecXHarvest.copy()
        vecPLaunch = best_vecPHarvest.copy()
        ###vecXLaunch = vecXHarvest.copy()
        ###vecPLaunch = vecPHarvest.copy()
        #
        vecDLaunch = vecXLaunch - vecXAnchor
        vecYLaunch = matQ.T @ vecDLaunch
        vecXPerp = vecDLaunch - ( matQ @ vecYLaunch )
        vecGammaLaunch = vecGammaFit + ( matHMod @ vecYLaunch )
        fLaunch = fFit + ( vecYLaunch @ vecGammaFit ) + (( vecYLaunch @ vecGammaLaunch )/2.0)
        vecT = matQ.T @ vecPLaunch
        vecPPerp = vecPLaunch - ( matQ @ vecT )
        assert linalg.norm( vecGammaLaunch ) > 0.0
        coeffPG = ( vecGammaLaunch @ vecT ) / ( vecGammaLaunch @ vecGammaLaunch )
        vecGammaPerp = vecT - ( coeffPG * vecGammaLaunch )
        testDecomp = True
        if ( testDecomp ):
            assert reldiff( vecXLaunch, vecXAnchor + (matQ @ vecYLaunch) + vecXPerp ) < MYEPS
            assert reldiff( vecPLaunch, (matQ @ ( (coeffPG*vecGammaLaunch) + vecGammaPerp )) + vecPPerp ) <= MYEPS
        if ( coeffPG > 0.0 ):
            coeffPG = 0.0
        
        # Calculate step.
        # TODO.
        # Placeholder: take a Newton step without any trust region.
        #vecZ = linalg.solve( matHMod, -vecGammaLaunch )
        vecLambdaCurve = vecLambdaMod  # Shallow copy / reference only / DO NOT MODIFY!
        vecLambdaObjf = vecLambdaOrig  # Shallow copy / reference only / DO NOT MODIFY!
        qnj_fMin = -0.01*fFit
        #msg( 'caps = ', qnj_sMax, qnj_dMax )
        vecZ = levsol( fFit, vecPhi, matPsi, vecLambdaCurve, qnj_sMax, vecS, qnj_dMax, vecLambdaObjf, qnj_fMin )
        ###vecZ = np.zeros( sizeK )
        
        # Generate new seed.
        vecYNew = vecYLaunch + vecZ
        vecGammaNew = vecGammaLaunch + ( matHMod @ vecZ )
        fNew = fLaunch + ( vecZ @ vecGammaLaunch ) + (( vecZ @ vecGammaNew )/2.0)
        assert fLaunch > 0.0
        assert linalg.norm(vecGammaLaunch/vecS) > 0.0
        alphaF = fNew / fLaunch
        alphaG = linalg.norm( vecGammaNew / vecS ) / linalg.norm( vecGammaLaunch / vecS )
        vecDelta = matQ @ vecZ
        #
        vecX[:] = vecXAnchor + ( matQ @ vecYNew ) + ( alphaG * vecXPerp )
        vecP[:] = (matQ @ ( (coeffPG*vecGammaNew) + (alphaG*vecGammaPerp) )) + (alphaF*vecPPerp)
        
        vecXSeed[:] = vecX[:]
        vecPSeed[:] = vecP[:]
        qnj_havePrev = True
        qnj_dPrev = linalg.norm( vecDelta )
        qnj_sPrev = linalg.norm( matS @ vecZ )
        
        continue
        # DRaburn 2023-01-27, pytorchDanmo: Pre-QNJ
        
        # TAKE STEP
        #optimizer.step()
        sxsolve_step = (sxsolve_momentum * sxsolve_step) - (sxsolve_lr * sxsolve_gradloc)
        sxsolve_x += sxsolve_step
        if (sxsolve_f<sxsolve_f_min):
            sxsolve_f_min = sxsolve_f
        
        # print statistics
        running_loss += loss.item()
        running_feval_count += 1
        #if i % 2000 == 1999:    # print every 2000 mini-batches
        #print_interval = 2000
        print_interval = 200
        if i % print_interval == (print_interval-1):    # print every few mini-batches
        #if (time.time()-running_time0>3.0):
            print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss * 1.0 / running_feval_count:.3f}')
            print(f"  sxsolve_f_min = {sxsolve_f_min:17e}")
            print(f"  Avg time per feval = { (time.time()-running_time0)*1.0 / running_feval_count }")
            print(f"  Elapsed time = {time.time()-start_time}s")
            running_loss = 0.0
            running_time0 = time.time()
            running_feval_count = 0
    if ( not doMainLoop ):
        break

print('Finished Training Demo')
# Look at results.
progLogSymbol = 'F'
msg(
  f' {superPtCount:4d} ({badCount:4d}X), {fevalCount:7d}:',
  f' {sizeK:3d} / {numRecords:3d}:'
  f'  {linalg.norm( superPt_vecX - vecX0 ):8.2E};',
  f'  {linalg.norm( vecXHarvest - vecXSeed ):8.2E}, {qnj_dPrev:9.2E} / {qnj_dMax:9.2E};',
  f'  {superPt_f:8.2E};',
  f'  {linalg.norm(best_vecG):8.2E} +/- {superPt_gVar:8.2E}',
  progLogSymbol )
vecXF = best_vecX
vecGF = best_vecG
fF = best_f
#fF, vecGF = funcFG( vecXF )
msg( '||vecXF - vecX0|| = ', linalg.norm( vecXF - vecX0 ) )
msg( 'fF = ', fF )
msg( '||vecGF|| = ', linalg.norm(vecGF) )
