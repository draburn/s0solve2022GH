# DRaburn 2023-02-16
# Based on "Training a Classifier":
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image

# Import.
import danutil
from danutil import msg
import torch
import torchvision
import numpy as np

################################################################################
#
# PERFORM INITIALIZATION ON IMPORT
#

# Set base params.
# Modify these as needed.
report_initialization = True
torch_seed = 0
CIFAR10_root = '../../dat/CIFAR10data'
trainset_size = 5000
batch_size = 500
placeholder_lr = 0.0
placeholder_momentum = 0.0
num_workers = 2 # num_workers impacts threads for reading from disk, or something.

# Our neural network.
class Net(torch.nn.Module):
	def __init__(self):
		super().__init__()
		self.conv1 = torch.nn.Conv2d(3, 6, 5)
		self.pool = torch.nn.MaxPool2d(2, 2)
		self.conv2 = torch.nn.Conv2d(6, 16, 5)
		self.fc1 = torch.nn.Linear(16 * 5 * 5, 120)
		self.fc2 = torch.nn.Linear(120, 84)
		self.fc3 = torch.nn.Linear(84, 10)
	def forward(self, x):
		x = self.pool(torch.nn.functional.relu(self.conv1(x)))
		x = self.pool(torch.nn.functional.relu(self.conv2(x)))
		x = torch.flatten(x, 1)
		x = torch.nn.functional.relu(self.fc1(x))
		x = torch.nn.functional.relu(self.fc2(x))
		x = self.fc3(x)
		return x

# Stuff for injection/extraction of Torch data.
# DRaburn 2023-02-16:
#  I'd imagine a proper Pytorch API exists, but I couldn't find information about it.
#  A "brute-force comprehensive Python" approach is likely also possible, but, mise;
#   this code is limited to the specific NN defined above.
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
def get_numel_list(this_net):
	numel_list = [this_net.conv1.bias.data.numel()]
	numel_list.append(this_net.conv1.weight.data.numel())
	numel_list.append(this_net.conv2.bias.data.numel())
	numel_list.append(this_net.conv2.weight.data.numel())
	numel_list.append(this_net.fc1.bias.data.numel())
	numel_list.append(this_net.fc1.weight.data.numel())
	numel_list.append(this_net.fc2.bias.data.numel())
	numel_list.append(this_net.fc2.weight.data.numel())
	numel_list.append(this_net.fc3.bias.data.numel())
	numel_list.append(this_net.fc3.weight.data.numel())
	return numel_list
def get_cumel_list(this_net):
	 numel_list = get_numel_list(this_net)
	 cumel_list = [0];
	 for n in range(len(numel_list)):
		  cumel_list.append(cumel_list[n] + numel_list[n])
	 return cumel_list
def init_x_from_net(this_net):
	cumel_list = get_cumel_list(this_net)
	this_x = np.zeros(cumel_list[-1], dtype=np.float32)
	this_x[cumel_list[0]:cumel_list[1]] = np.reshape(this_net.conv1.bias.data.numpy(), -1)
	this_x[cumel_list[1]:cumel_list[2]] = np.reshape(this_net.conv1.weight.data.numpy(), -1)
	this_x[cumel_list[2]:cumel_list[3]] = np.reshape(this_net.conv2.bias.data.numpy(), -1)
	this_x[cumel_list[3]:cumel_list[4]] = np.reshape(this_net.conv2.weight.data.numpy(), -1)
	this_x[cumel_list[4]:cumel_list[5]] = np.reshape(this_net.fc1.bias.data.numpy(), -1)
	this_x[cumel_list[5]:cumel_list[6]] = np.reshape(this_net.fc1.weight.data.numpy(), -1)
	this_x[cumel_list[6]:cumel_list[7]] = np.reshape(this_net.fc2.bias.data.numpy(), -1)
	this_x[cumel_list[7]:cumel_list[8]] = np.reshape(this_net.fc2.weight.data.numpy(), -1)
	this_x[cumel_list[8]:cumel_list[9]] = np.reshape(this_net.fc3.bias.data.numpy(), -1)
	this_x[cumel_list[9]:cumel_list[10]] = np.reshape(this_net.fc3.weight.data.numpy(), -1)
	return this_x
def init_g_from_net(this_net):
	cumel_list = get_cumel_list(this_net)
	this_grad = np.zeros(cumel_list[-1], dtype=np.float32)
	this_grad[cumel_list[0]:cumel_list[1]] = np.reshape(this_net.conv1.bias.grad.numpy(), -1)
	this_grad[cumel_list[1]:cumel_list[2]] = np.reshape(this_net.conv1.weight.grad.numpy(), -1)
	this_grad[cumel_list[2]:cumel_list[3]] = np.reshape(this_net.conv2.bias.grad.numpy(), -1)
	this_grad[cumel_list[3]:cumel_list[4]] = np.reshape(this_net.conv2.weight.grad.numpy(), -1)
	this_grad[cumel_list[4]:cumel_list[5]] = np.reshape(this_net.fc1.bias.grad.numpy(), -1)
	this_grad[cumel_list[5]:cumel_list[6]] = np.reshape(this_net.fc1.weight.grad.numpy(), -1)
	this_grad[cumel_list[6]:cumel_list[7]] = np.reshape(this_net.fc2.bias.grad.numpy(), -1)
	this_grad[cumel_list[7]:cumel_list[8]] = np.reshape(this_net.fc2.weight.grad.numpy(), -1)
	this_grad[cumel_list[8]:cumel_list[9]] = np.reshape(this_net.fc3.bias.grad.numpy(), -1)
	this_grad[cumel_list[9]:cumel_list[10]] = np.reshape(this_net.fc3.weight.grad.numpy(), -1)
	return this_grad
# DRaburn 2023-02-19:
#  From stack overflow, a start at a more general interface:
#total_params = 0
#for name, parameter in Net().named_parameters():
#	if not parameter.requires_grad: continue
#	params = parameter.numel()
#	print(name, params)
#	total_params+=params
#print(f"Total Trainable Params: {total_params}")
if (report_initialization):
	msg(f'Initializing pytorchCIFAR10demo...')
	msg(f'  torch_seed = {torch_seed}')
	msg(f'  trainset_size = {trainset_size}')
	msg(f'  CIFAR10_root = "{CIFAR10_root}"')

# Init Torch, NN, etc.
torch.manual_seed(torch_seed)
transform = torchvision.transforms.Compose([
  torchvision.transforms.ToTensor(),
  torchvision.transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])
if (trainset_size > 0):
	full_dataset = torchvision.datasets.CIFAR10(root=CIFAR10_root, train=True, download=False, transform=transform)
	full_num_samples = len(full_dataset)
	msg(f'*** WARNING: Only using {trainset_size} / {full_num_samples} samples.')
	trainset, dummyset = torch.utils.data.random_split(full_dataset, [trainset_size, full_num_samples-trainset_size])
else:
	trainset = torchvision.datasets.CIFAR10(root=CIFAR10_root, train=True, download=False, transform=transform)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size, shuffle=True, num_workers=num_workers)
classes = ('plane', 'car', 'bird', 'cat', 'deer', 'dog', 'frog', 'horse', 'ship', 'truck')
net = Net()
torch_optim_SGD = torch.optim.SGD(net.parameters(), lr=placeholder_lr, momentum=placeholder_momentum)
loss_criterion = torch.nn.CrossEntropyLoss()
size_list = get_size_list(net)
cumel_list = get_cumel_list(net)
num_unknowns = cumel_list[-1]
num_samples = len(trainset)
num_batches_in_epoch = round( num_samples / batch_size )

# Initialize our "x".
# DRaburn 2023-02-16:
#  There may be a better way to do this, but, I don't know what it is.
#  First, we'll create a "shared_vecX" and populate it.
#  Then, we make the weights (and biases) in the net point to our "shared_vecX".
shared_vecX = init_x_from_net(net)
net.conv1.bias.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[0]:cumel_list[1]],size_list[0]))
net.conv1.weight.data = torch.from_numpy(np.reshape(shared_vecX[cumel_list[1]:cumel_list[2]],size_list[1]))
net.conv2.bias.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[2]:cumel_list[3]],size_list[2]))
net.conv2.weight.data = torch.from_numpy(np.reshape(shared_vecX[cumel_list[3]:cumel_list[4]],size_list[3]))
net.fc1.bias.data     = torch.from_numpy(np.reshape(shared_vecX[cumel_list[4]:cumel_list[5]],size_list[4]))
net.fc1.weight.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[5]:cumel_list[6]],size_list[5]))
net.fc2.bias.data     = torch.from_numpy(np.reshape(shared_vecX[cumel_list[6]:cumel_list[7]],size_list[6]))
net.fc2.weight.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[7]:cumel_list[8]],size_list[7]))
net.fc3.bias.data     = torch.from_numpy(np.reshape(shared_vecX[cumel_list[8]:cumel_list[9]],size_list[8]))
net.fc3.weight.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[9]:cumel_list[10]],size_list[9]))
# I'm not sure what happens to original memory in net, and I don't particularly care.

# Initialize our gradient.
# DRaburn 2023-02-16:
#  There may be a better way to do this, but, I don't know what it is.
#  This approach is basically the same as for shared_vecX, but,
#   it seems we need to do at least one gradient calculation first.
torch_optim_SGD.zero_grad()
for batch_index, batch_data in enumerate(trainloader, 0):
	batch_inputs, batch_labels = batch_data
	batch_outputs = net(batch_inputs)
	batch_loss = loss_criterion(batch_outputs, batch_labels)
	batch_loss.backward()
	break
shared_vecG = init_g_from_net(net)
net.conv1.bias.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[0]:cumel_list[1]],size_list[0]))
net.conv1.weight.grad = torch.from_numpy(np.reshape(shared_vecG[cumel_list[1]:cumel_list[2]],size_list[1]))
net.conv2.bias.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[2]:cumel_list[3]],size_list[2]))
net.conv2.weight.grad = torch.from_numpy(np.reshape(shared_vecG[cumel_list[3]:cumel_list[4]],size_list[3]))
net.fc1.bias.grad     = torch.from_numpy(np.reshape(shared_vecG[cumel_list[4]:cumel_list[5]],size_list[4]))
net.fc1.weight.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[5]:cumel_list[6]],size_list[5]))
net.fc2.bias.grad     = torch.from_numpy(np.reshape(shared_vecG[cumel_list[6]:cumel_list[7]],size_list[6]))
net.fc2.weight.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[7]:cumel_list[8]],size_list[7]))
net.fc3.bias.grad     = torch.from_numpy(np.reshape(shared_vecG[cumel_list[8]:cumel_list[9]],size_list[8]))
net.fc3.weight.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[9]:cumel_list[10]],size_list[9]))

if (report_initialization):
	msg(f'  num_unknowns = {num_unknowns}')
	msg(f'  num_samples = {num_samples}')
	msg(f'  batch_size = {batch_size}')
	if ( num_batches_in_epoch * batch_size != num_samples ):
		msg(f'*** WARNING: The number of samples is not divisible by the batch size. ***')
	msg(f'  num_workers = {num_workers}')
	msg(f'  placeholder_lr = {placeholder_lr:0.18E}')
	msg(f'  placeholder_momentum = {placeholder_momentum:0.18E}')
	msg(f'  num_batches_in_epoch = {num_batches_in_epoch}')
	msg('Finished pytorchCIFAR10demo initialization.')

################################################################################
#
# INTERFACE ROUTINES FOR EXTERNAL OPTIMIZER.
#
def genVecX0():
	return shared_vecX.copy()
# End get_vecX().

class evalSGD_prm():
	def __init__(self):
		self.learningRate = 1.0e-3
		self.momentumCoeff = 0.9
		self.doStats = True
		self.doStore = False
		self.storageSize = num_batches_in_epoch + 1 # +1 to be safe.
	def dump(self):
		msg(f'Begin evalSGD_prm.dump()...')
		msg(f'self = {self}')
		msg(f'  learningRate = {self.learningRate}')
		msg(f'  momentumCoeff = {self.momentumCoeff}')
		msg(f'  doStats = {self.doStats}')
		msg(f'  doStore = {self.doStore}')
		msg(f'  storageSize = {self.storageSize}')
		msg(f'End evalSGD_prm.dump().')
	def copy(self):
		cp = evalSGD_prm()
		cp.learningRate = self.learningRate
		cp.momentumCoeff = self.momentumCoeff
		cp.doStats = self.doStats
		cp.doStore = self.doStore
		cp.storageSize = self.storageSize
		return cp
class evalSGD_statsDat():
	def __init__(self):
		self.numSteps = 0
		self.avg_vecX = np.zeros(num_unknowns)
		self.var_vecX = np.zeros(num_unknowns)
		self.avg_f = 0.0
		self.var_f = 0.0
		self.avg_xtg = 0.0
		self.var_xtg = 0.0
		self.avg_vecG = np.zeros(num_unknowns)
		self.var_vecG = np.zeros(num_unknowns)
		self.hessfit_f = 0.0
	def dump(self):
		msg(f'Begin evalSGD_statsDat.dump()...')
		msg(f'self = {self}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  avg_vecX = {self.avg_vecX}')
		msg(f'  var_vecX = {self.var_vecX}')
		msg(f'  avg_f = {self.avg_f}')
		msg(f'  var_f = {self.var_f}')
		msg(f'  avg_xtg = {self.avg_xtg}')
		msg(f'  var_xtg = {self.var_xtg}')
		msg(f'  avg_vecG = {self.avg_vecG}')
		msg(f'  var_vecG = {self.var_vecG}')
		msg(f'  hessfit_f = {self.hessfit_f}')
		msg(f'End evalSGD_statsDat.dump().')
	def absorb(self, vecX, f, vecG):
		self.numSteps +=1
		self.avg_vecX[:] += vecX[:]
		self.var_vecX[:] += vecX[:]**2
		self.avg_f += f
		self.var_f += f**2
		xtg = vecX @ vecG
		self.avg_xtg += xtg
		self.var_xtg += xtg**2
		self.avg_vecG[:] += vecG[:]
		self.var_vecG[:] += vecG[:]**2
	def finalize(self):
		numSteps = self.numSteps
		self.avg_vecX[:] /= numSteps
		self.var_vecX[:] /= numSteps
		self.avg_f /= numSteps
		self.var_f /= numSteps
		self.avg_xtg /= numSteps
		self.var_xtg /= numSteps
		self.avg_vecG[:] /= numSteps
		self.var_vecG[:] /= numSteps
		self.var_vecX = danutil.var(self.avg_vecX, self.var_vecX)
		self.var_f = danutil.var(self.avg_f, self.var_f)
		self.var_xtg = danutil.var(self.avg_xtg, self.var_xtg)
		self.var_vecG = danutil.var(self.avg_vecG, self.var_vecG)
		self.hessfit_f = self.avg_f - (( self.avg_xtg - (self.avg_vecX @ self.avg_vecG))/2.0)
class evalSGD_storeDat():
	def __init__(self, storageSize):
		self.numSteps = 0
		self.matX = np.zeros((num_unknowns, storageSize))
		self.vecF = np.zeros(storageSize)
		self.matG = np.zeros((num_unknowns, storageSize))
	def dump(self):
		msg(f'Begin evalSGD_storeDat.dump()...')
		msg(f'self = {self}')
		msg(f'  numSteps = {self.numSteps}')
		msg(f'  matX = ...\n{self.matX}')
		msg(f'  vecF = {self.vecF}')
		msg(f'  matG = ...\n{self.matG}')
		msg(f'End evalSGD_storeDat.dump().')
	def absorb(self, vecX, f, vecG):
		self.matX[:,self.numSteps] = vecX[:]
		self.vecF[self.numSteps] = f
		self.matG[:,self.numSteps] = vecG[:]
		self.numSteps += 1
	def finalize(self):
		pass
class evalSGD_datOut():
	def __init__(self, doStats, doStore, storageSize):
		self.statsDat = None # Unless...
		self.storeDat = None # Unless...
		if (doStats):
			self.statsDat = evalSGD_statsDat()
		if (doStore):
			self.storeDat = evalSGD_storeDat(storageSize)
	def dump(self):
		msg(f'Begin evalSGD_datOut.dump()...')
		msg(f'self = {self}')
		msg(f'  statsDat = {self.statsDat}')
		if (not(None == self.statsDat)):
			self.statsDat.dump()
		msg(f'  storeDat = {self.storeDat}')
		if (not(None == self.storeDat)):
			self.storeDat.dump()
		msg(f'End evalSGD_datOut.dump().')
	def absorb(self, vecX, f, vecG):
		if (not(None == self.statsDat)):
			self.statsDat.absorb(vecX, f, vecG)
		if (not(None == self.storeDat)):
			self.storeDat.absorb(vecX, f, vecG)
	def finalize(self):
		if (not(None == self.statsDat)):
			self.statsDat.finalize()
		if (not(None == self.storeDat)):
			self.storeDat.finalize()
# DRaburn 2023-02-23:
# The concepts for evalSGD(), evalFG(), and evalF() here are a bit inconsistent with
#  what I originally had in demoproblem0221.py regarding "numSteps":
#  consistency would be for evalSGD to go for numSteps epochs.
# Oh well. This is just for proof-of-principle anyway.
def evalSGD( vecXSeed, vecPSeed=np.zeros(num_unknowns), prm=evalSGD_prm() ):
	vecX = vecXSeed.copy()
	vecP = vecPSeed.copy()
	stepCount = 0
	avg_f = 0.0
	datOut = evalSGD_datOut(prm.doStats, prm.doStore, prm.storageSize)
	for batch_index, batch_data in enumerate(trainloader, 0):
		shared_vecX[:] = vecX[:]
		torch_optim_SGD.zero_grad()
		batch_inputs, batch_labels = batch_data
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		batch_f = batch_loss.item()
		datOut.absorb( vecX, batch_f, shared_vecG )
		stepCount += 1
		avg_f += batch_f
		vecP[:] = (prm.momentumCoeff * vecP[:]) - (prm.learningRate * shared_vecG[:])
		vecX[:] += vecP[:]
	datOut.finalize()
	avg_f /= stepCount
	return vecX, vecP, avg_f, datOut
# End evalSGD().

def evalFG( vecX ):
	stepCount = 0
	avg_f = 0.0
	avg_vecG = np.zeros(num_unknowns)
	for batch_index, batch_data in enumerate(trainloader, 0):
		shared_vecX[:] = vecX[:]
		torch_optim_SGD.zero_grad()
		batch_inputs, batch_labels = batch_data
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		batch_f = batch_loss.item()
		stepCount += 1
		avg_f += batch_f
		avg_vecG[:] += shared_vecG[:]
	avg_f /= stepCount
	avg_vecG /= stepCount
	return avg_f, avg_vecG
# End def evalFG().

def evalF( vecX ):
	avg_f, _ = evalFG(vecX)
	return avg_f
# End def evalFG().
