# DRaburn 2023-02-16
# Based on "Training a Classifier":
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image

# Import.
import inspect
import torch
import torchvision
import numpy as np
main_frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{main_frame.f_lineno:05d}]', *arguments, **keywords)

################################################################################
#
# PERFORM INITIALIZATION ON IMPORT
#

# Set base params.
# Modify these as needed.
report_initialization = True
torch_seed = 0
CIFAR10_root = '../../dat/CIFAR10data'
batch_size = 500
placeholder_lr = 0.0
placeholder_momentum = 0.0
num_workers = 2

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

# Report(?)
if (report_initialization):
	msg(f'Initializing pytorchCIFAR10demo...')
	msg(f'  torch_seed = {torch_seed}')
	msg(f'  CIFAR10_root = "{CIFAR10_root}"')

# Init Torch, NN, etc.
torch.manual_seed(torch_seed)
transform = torchvision.transforms.Compose([
  torchvision.transforms.ToTensor(),
  torchvision.transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])
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
#
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

def get_vecX():
	return shared_vecX.copy()
# End get_vecX().

class eval_epoch_sgd_datOut():
	def __init__(self):
		self.vecXSeed = np.zeros(num_unknowns)
		self.vecPSeed = np.zeros(num_unknowns)
		self.batch_count = 0
		self.avg_f = 0.0
		self.avg_vecX = np.zeros(num_unknowns)
		self.avg_vecG = np.zeros(num_unknowns)
		self.avg_vecP = np.zeros(num_unknowns)
		self.avg_fSq = 0.0
		self.avg_vecXSq = np.zeros(num_unknowns)
		self.avg_vecGSq = np.zeros(num_unknowns)
		self.avg_vecPSq = np.zeros(num_unknowns)
		self.vecXHarvest = np.zeros(num_unknowns)
		self.vecPHarvest = np.zeros(num_unknowns)
	def dump(self):
		msg('Begin eval_epoch_sgd_datOut.dump()...')
		print(f'  batch_count = {self.batch_count}')
		print(f'  f: avg = {self.avg_f:25.18E}, sqAvg = {self.avg_fSq:25.18E}')
		print(f'  vecG: avg[0] = {self.avg_vecG[0]:25.18E}, avgSq[0] = {self.avg_vecGSq[0]:25.18E}, ||avg|| = {np.linalg.norm(self.avg_vecG):25.18E}, ||avgSq|| = {np.linalg.norm(self.avg_vecGSq):25.18E}')
		print(f'  vecX: avg[0] = {self.avg_vecX[0]:25.18E}, avgSq[0] = {self.avg_vecXSq[0]:25.18E}, ||avg|| = {np.linalg.norm(self.avg_vecX):25.18E}, ||avgSq|| = {np.linalg.norm(self.avg_vecXSq):25.18E}, hms[0] = {self.vecXHarvest[0]-self.vecXSeed[0]:25.18E}, ||hms|| = {np.linalg.norm(self.vecXHarvest-self.vecXSeed):25.18E}')
		print(f'  vecP: avg[0] = {self.avg_vecP[0]:25.18E}, avgSq[0] = {self.avg_vecPSq[0]:25.18E}, ||avg|| = {np.linalg.norm(self.avg_vecP):25.18E}, ||avgSq|| = {np.linalg.norm(self.avg_vecPSq):25.18E}, hms[0] = {self.vecPHarvest[0]-self.vecPSeed[0]:25.18E}, ||hms|| = {np.linalg.norm(self.vecPHarvest-self.vecPSeed):25.18E}')
		print(f'  vecPSeed:    self[0] = {self.vecPSeed[0]:25.18E}, ||self|| = {np.linalg.norm(self.vecPSeed):25.18E}')
		print(f'  vecPHarvest: self[0] = {self.vecPHarvest[0]:25.18E}, ||self|| = {np.linalg.norm(self.vecPHarvest):25.18E}')
		msg('End eval_epoch_sgd_datOut.dump().')
# End eval_epoch_sgd_datOut().

def eval_epoch_sgd( vecXSeed, vecPSeed, learning_rate, momentum_coefficient, num_batches ):
	vecX = vecXSeed.copy()
	vecP = vecPSeed.copy()
	batch_count = 0
	dat = eval_epoch_sgd_datOut()
	dat.vecXSeed = vecXSeed.copy()
	dat.vecPSeed = vecPSeed.copy()
	for batch_index, batch_data in enumerate(trainloader, 0):
		# Accumulate pre-feval data.
		batch_count += 1
		dat.avg_vecX[:] += vecX[:]
		dat.avg_vecP[:] += vecP[:]
		dat.avg_vecXSq[:] += vecX[:]**2
		dat.avg_vecPSq[:] += vecP[:]**2
		# Do feval.
		shared_vecX[:] = vecX[:]
		torch_optim_SGD.zero_grad()
		batch_inputs, batch_labels = batch_data
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		batch_f = batch_loss.item()
		# Accumulate post-feval data.
		dat.avg_f += batch_f
		dat.avg_vecG[:] += shared_vecG[:]
		dat.avg_fSq += batch_f**2
		dat.avg_vecGSq[:] += shared_vecG[:]**2
		# Take step.
		vecP[:] = (momentum_coefficient * vecP[:]) - (learning_rate * shared_vecG[:])
		vecX[:] += vecP[:]
		#
		if ((num_batches > 0) and (batch_index+1 >= num_batches)):
			break
	if ((num_batches > 0) and (batch_index+1 < num_batches)):
		msg(f'*** WARNING: {num_batches} batches were requested but only {batch_index+1} were performed. ***')
	# End batch loop.
	dat.batch_count = batch_count
	dat.avg_f /= batch_count
	dat.avg_vecX[:] /= batch_count
	dat.avg_vecG[:] /= batch_count
	dat.avg_vecP[:] /= batch_count
	dat.avg_fSq /= batch_count
	dat.avg_vecXSq[:] /= batch_count
	dat.avg_vecGSq[:] /= batch_count
	dat.avg_vecPSq[:] /= batch_count
	dat.vecXHarvest[:] = vecX[:]
	dat.vecPHarvest[:] = vecP[:]
	return ( dat.avg_f, dat.avg_vecG, dat )
# End eval_epoch_sgd__x_forced().
