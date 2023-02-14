# DRaburn 2023-02-13:
# This code is based on the pytorch "Training a Classifier" tutorial at:
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image
# This code includes considerable "hacks" to interface pytorch with my algorithms.


#
# INIT
#
# Import stuff.
import inspect
import numpy
import time
import torch
import torchvision as tv
import torchvision.transforms as tv_transforms
import torch.nn as nn
import torch.optim
import torch.nn.functional as F

# Init basics.
start_time = time.time()
frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{frame.f_lineno:05d}]', *arguments, **keywords)
msg(f'time.asctime() = {time.asctime()}')

# Set universal parameters.
torch_seed = 0
batch_size = 500
CIFAR10_root = './data'
fname_x0 = ''
dtype_x0 = numpy.float32
fname_p0 = ''
dtype_p0 = numpy.float32
msg(f'torch_seed = {torch_seed}')
msg(f'batch_size = {batch_size}')
msg(f'CIFAR10_root = "{CIFAR10_root}"')

torch.manual_seed(torch_seed)
transform = tv_transforms.Compose([tv_transforms.ToTensor(), tv_transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])
trainset = tv.datasets.CIFAR10(root=CIFAR10_root, train=True, download=False, transform=transform)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size, shuffle=True, num_workers=2)
classes = ('plane', 'car', 'bird', 'cat', 'deer', 'dog', 'frog', 'horse', 'ship', 'truck')
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

# DRaburn 2023-02-13:
# The following is hack-ish code to freely interface (inject/extract) with pytorch;
#  a "proper Pytorch" approach likely exists, but I couldn't find it readliy;
#  a "brute-force comprehensive Python" approach is likely also possible, but, mise.
# This is copied from older code.
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
	this_x = numpy.zeros(cumel_list[-1], dtype=numpy.float32)
	this_x[cumel_list[0]:cumel_list[1]] = numpy.reshape(this_net.conv1.bias.data.numpy(), -1)
	this_x[cumel_list[1]:cumel_list[2]] = numpy.reshape(this_net.conv1.weight.data.numpy(), -1)
	this_x[cumel_list[2]:cumel_list[3]] = numpy.reshape(this_net.conv2.bias.data.numpy(), -1)
	this_x[cumel_list[3]:cumel_list[4]] = numpy.reshape(this_net.conv2.weight.data.numpy(), -1)
	this_x[cumel_list[4]:cumel_list[5]] = numpy.reshape(this_net.fc1.bias.data.numpy(), -1)
	this_x[cumel_list[5]:cumel_list[6]] = numpy.reshape(this_net.fc1.weight.data.numpy(), -1)
	this_x[cumel_list[6]:cumel_list[7]] = numpy.reshape(this_net.fc2.bias.data.numpy(), -1)
	this_x[cumel_list[7]:cumel_list[8]] = numpy.reshape(this_net.fc2.weight.data.numpy(), -1)
	this_x[cumel_list[8]:cumel_list[9]] = numpy.reshape(this_net.fc3.bias.data.numpy(), -1)
	this_x[cumel_list[9]:cumel_list[10]] = numpy.reshape(this_net.fc3.weight.data.numpy(), -1)
	return this_x
def init_grad_from_net(this_net):
	cumel_list = get_cumel_list(this_net)
	this_grad = numpy.zeros(cumel_list[-1], dtype=numpy.float32)
	this_grad[cumel_list[0]:cumel_list[1]] = numpy.reshape(this_net.conv1.bias.grad.numpy(), -1)
	this_grad[cumel_list[1]:cumel_list[2]] = numpy.reshape(this_net.conv1.weight.grad.numpy(), -1)
	this_grad[cumel_list[2]:cumel_list[3]] = numpy.reshape(this_net.conv2.bias.grad.numpy(), -1)
	this_grad[cumel_list[3]:cumel_list[4]] = numpy.reshape(this_net.conv2.weight.grad.numpy(), -1)
	this_grad[cumel_list[4]:index_list[5]] = numpy.reshape(this_net.fc1.bias.grad.numpy(), -1)
	this_grad[cumel_list[5]:cumel_list[6]] = numpy.reshape(this_net.fc1.weight.grad.numpy(), -1)
	this_grad[cumel_list[6]:cumel_list[7]] = numpy.reshape(this_net.fc2.bias.grad.numpy(), -1)
	this_grad[cumel_list[7]:cumel_list[8]] = numpy.reshape(this_net.fc2.weight.grad.numpy(), -1)
	this_grad[cumel_list[8]:cumel_list[9]] = numpy.reshape(this_net.fc3.bias.grad.numpy(), -1)
	this_grad[cumel_list[9]:cumel_list[10]] = numpy.reshape(this_net.fc3.weight.grad.numpy(), -1)
	return this_grad

net = Net()
numel_list = get_numel_list(net)
cumel_list = get_cumel_list(net)
sxsolve_vecX = init_x_from_net(net)
sizeX = cumel_list[-1]
msg(f'numel_list = {numel_list}')
msg(f'cumel_list = {cumel_list}')
msg(f'sizeX = {sizeX}')
sxsolve_vecP = numpy.zeros(sizeX, dtype=numpy.float32)

# DRaburn 2022:
# The following appears to make data point to the same memory as sxsolve_vecX.
# Documentation seems to suggest this is guaranteed for torch.from_numpy(),
#  though I do not see a guarantee for numpy.reshape().
# Not sure what happens to the memory originally allocated for data, and don't care.
net.conv1.bias.data	= torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[0]:cumel_list[1]],numel_list[0]))
net.conv1.weight.data = torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[1]:cumel_list[2]],numel_list[1]))
net.conv2.bias.data	= torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[2]:cumel_list[3]],numel_list[2]))
net.conv2.weight.data = torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[3]:cumel_list[4]],numel_list[3]))
net.fc1.bias.data	= torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[4]:cumel_list[5]],numel_list[4]))
net.fc1.weight.data = torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[5]:cumel_list[6]],numel_list[5]))
net.fc2.bias.data	= torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[6]:cumel_list[7]],numel_list[6]))
net.fc2.weight.data = torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[7]:cumel_list[8]],numel_list[7]))
net.fc3.bias.data	= torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[8]:cumel_list[9]],numel_list[8]))
net.fc3.weight.data = torch.from_numpy(numpy.reshape(sxsolve_vecX[cumel_list[9]:cumel_list[10]],numel_list[9]))
# We'll wait to init vecG until it's calculated, to get the pointers right.
do_grad_init = True

# Read from disk, if necessary.
msg(f'Elapsed time = {time.time()-start_time}s.')
if (''!=fname_x0):
	msg(f'Reading x0 from disk using dtype = "{dtype_x0}".')
	foo = numpy.fromfile(fname_x0,dtype=dtype_x0)
	sxsolve_vecX[:] = foo[:]
	foo = []
if (''!=fname_p0):
	msg(f'Reading p0 from disk using dypte = "{dtype_p0}".')
	foo = numpy.fromfile(fname_p0,dtype=dtype_p0)
	sxsolve_vecP[:] = foo[:]
	foo = []
msg(f'Elapsed time = {time.time()-start_time}s.')

# Initialize placeholder optimizer, SGD, and loss function.
placeholder_epoch_limit = 5
placeholder_lr = 1.0
placeholder_momentum = 0.0
msg(f'placeholder_epoch_limit = {placeholder_epoch_limit}')
msg(f'placeholder_lr = {placeholder_lr:0.17E}')
msg(f'placeholder_momentum = {placeholder_momentum:0.17E}')
placeholder_optimizer = torch.optim.SGD(net.parameters(), lr=placeholder_lr, momentum=placeholder_momentum)
loss_criterion = torch.nn.CrossEntropyLoss()


#
# MAIN LOOP
#
msg(f'Elapsed time = {time.time()-start_time}s.')
msg(f'Starting main loop...')
print('')
print('[')
feval_count = 0
running_feval_count = 0
running_loss = 0.0
for epoch_index in range(placeholder_epoch_limit):
	for batch_index, batch_data in enumerate(trainloader, 0):
		batch_inputs, batch_labels = batch_data
		placeholder_optimizer.zero_grad()
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		feval_count += 1
		
		if (do_grad_init):
			sxsolve_vecG = init_grad_from_net(net)
			net.conv1.bias.grad	= torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[0]:cumel_list[1]],numel_list[0]))
			net.conv1.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[1]:cumel_list[2]],numel_list[1]))
			net.conv2.bias.grad	= torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[2]:cumel_list[3]],numel_list[2]))
			net.conv2.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[3]:cumel_list[4]],numel_list[3]))
			net.fc1.bias.grad	= torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[4]:cumel_list[5]],numel_list[4]))
			net.fc1.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[5]:cumel_list[6]],numel_list[5]))
			net.fc2.bias.grad	= torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[6]:cumel_list[7]],numel_list[6]))
			net.fc2.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[7]:cumel_list[8]],numel_list[7]))
			net.fc3.bias.grad	= torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[8]:cumel_list[9]],numel_list[8]))
			net.fc3.weight.grad = torch.from_numpy(numpy.reshape(sxsolve_vecG[cumel_list[9]:cumel_list[10]],numel_list[9]))			
			do_grad_init = False
		
		#optimizer.step()
		sxsolve_vecP[:] = (sxsolve_momentum * sxsolve_vecP[:]) - (sxsolve_lr * sxsolve_vecG[:])
		sxsolve_vecX += sxsolve_vecP
		running_loss += loss.item()
		running_feval_count += 1
		
	avg_loss = running_loss / running_feval_count
	print(f'[ {time.time()-start_time:10.2F}, {feval_count:10d}, {epoch_index+1:6d}, {avg_loss:0.17E}]')
	if (epoch_index % 100 == 0):
		sxsolve_vecX.tofile(f'out_vecX{epoch:06d}.np')
		sxsolve_vecP.tofile(f'out_vecP{epoch:06d}.np')
	running_feval_count = 0
	running_loss = 0.0
print('];')
print('')
print('Finished main loop.')
sxsolve_x.tofile(f'out_vecX_fin.np')
sxsolve_step.tofile(f'out_vecP_fin.np')
