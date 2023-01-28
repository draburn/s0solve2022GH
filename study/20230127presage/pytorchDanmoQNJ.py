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
    if ( dMax > 0.0 ):
        assert dPastMaxOfP(0.0) < 0.0
        if ( dPastMaxOfP(p1) > 0.0 ):
            p1New = optimize.bisect( dPastMaxOfP, 0.0, p1 )
            p1 = p1New
    # Ugh. Just require fMin.
    assert fTillMinOfP(0.0) > 0.0
    if ( fTillMinOfP(p1) < 0.0 ):
        p1New = optimize.bisect( fTillMinOfP, 0.0, p1 )
        p1New = p1
    return deltaYOfP( p1 )
# Init problem.
MYEPS = 1.0E-8
rngSeed = 0
sizeX = 50
sizeF = sizeX
msg( 'rngSeed = ', rngSeed )
msg( 'sizeX = ', sizeX )
msg( 'sizeF = ', sizeF )
rng = default_rng(rngSeed)
matA = rng.standard_normal(( sizeF, sizeX ))
#matA = np.diag(np.linspace(1.0,sizeX,sizeX))
matHCrit = np.zeros(( sizeX, sizeX ))
matHCrit[:,:] = matA.T @ matA
vecXCrit = rng.standard_normal(( sizeX ))
#vecXCrit = np.ones(( sizeX ))
fCrit = 10.0
noiseX = 1.0E-5
#noiseX = 0.0
# DRaburn 2023-01-27, pytorchDanmo: Not you.
#def funcFG( x ):
#    #d = x - vecXCrit
#    d = x - vecXCrit + noiseX*rng.standard_normal(( sizeX ))
#    g = matHCrit @ d
#    f = fCrit + (( d @ g )/2.0)
#    return ( f, g )
#vecX0 = np.zeros(( sizeX ))
#f0, vecG0 = funcFG( vecX0 )
#msg( 'f0 = ', f0 )
#msg( '||vecG0|| = ', linalg.norm(vecG0) )

# Init SGD solver.
# DRaburn 2023-01-27, pytorchDanmo: I failed to write the tolerances integrably.
#fBail = f0 * 1E8
fBail = 1.0E8
fevalLimit = 100000
# DRaburn 2023-01-27, pytorchDanmo: Gradient behavior is controlled by pre-existing (pytorchDemod) code.
#learningRate = 0.01
#learningRate = 0.001
#momentumFactor = 0.9
msg( 'fevalLimit = ', fevalLimit )
#msg( 'learningRate = ', learningRate )
#msg( 'momentumFactor = ', momentumFactor )
fevalCount = 0
# DRaburn 2023-01-27, pytorchDanmo: Hook-up to pre-existing code.
vecX = sxsolve_x
vecP = sxsolve_step

# Init superPt.
numFevalPerSuperPt = 100
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
running_fTot = 0.0
running_fSqTot = 0.0
running_xtgTot = 0.0
running_vecGTot = np.zeros(( sizeX ))
running_vecXTot = np.zeros(( sizeX ))
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
msg( 'coeff_best_minf = ', coeff_best_minf )
msg( 'coeff_best_best = ', coeff_best_best )
msg( 'coeff_best_curr = ', coeff_best_curr )
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
#useQNJ = False
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
qnj_sMax = 3.0
qnj_sMax_btCoeff = 0.1
qnj_sMax_ftCoeff = 2.0
qnj_dPrev = -1.0
qnj_dMax = -1.0
qnj_dMax_btCoeff = 0.1
qnj_dMax_ftCoeff = 2.0

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

print('Finished Training Demo')
