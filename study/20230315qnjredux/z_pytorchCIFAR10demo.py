# DRaburn 2023-03-21.
# Based on "Training a Classifier":
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image

# Import.
import time
import danutil
from danutil import msg
from danutil import msgtime
import torch
import torchvision
#import numpy as np
startTime = time.time()
msgtime()

################################################################################
#
# PERFORM INITIALIZATION ON IMPORT
#

# Set base params.
# Modify these as needed.
report_initialization = True
torch_seed = 0
CIFAR10_root = '../../dat/CIFAR10data'
#trainset_size = 5000
trainset_size = -1
batch_size = 500
optimizer_lr = 0.01
optimizer_momentum = 0.9
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
	if ( trainset_size != full_num_samples ):
		msg(f'*** WARNING: Only using {trainset_size} / {full_num_samples} samples.')
	trainset, dummyset = torch.utils.data.random_split(full_dataset, [trainset_size, full_num_samples-trainset_size])
else:
	trainset = torchvision.datasets.CIFAR10(root=CIFAR10_root, train=True, download=False, transform=transform)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size, shuffle=True, num_workers=num_workers)
classes = ('plane', 'car', 'bird', 'cat', 'deer', 'dog', 'frog', 'horse', 'ship', 'truck')
net = Net()
torch_optim_SGD = torch.optim.SGD(net.parameters(), lr=optimizer_lr, momentum=optimizer_momentum)
loss_criterion = torch.nn.CrossEntropyLoss()
num_samples = len(trainset)
num_batches_in_epoch = round( num_samples / batch_size )
num_epochs = 500

if (report_initialization):
	msg(f'  num_samples = {num_samples}')
	msg(f'  batch_size = {batch_size}')
	if ( num_batches_in_epoch * batch_size != num_samples ):
		msg(f'*** WARNING: The number of samples is not divisible by the batch size. ***')
	msg(f'  num_workers = {num_workers}')
	msg(f'  optimizer_lr = {optimizer_lr:0.18E}')
	msg(f'  optimizer_momentum = {optimizer_momentum:0.18E}')
	msg(f'  num_batches_in_epoch = {num_batches_in_epoch}')
	msg('Finished initialization.')

msgtime()
# Solve.
print('[')
for epoch_index in range(num_epochs):
	avg_f = 0.0
	batch_count = 0
	for batch_index, batch_data in enumerate(trainloader, 0):
		torch_optim_SGD.zero_grad()
		batch_inputs, batch_labels = batch_data
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		avg_f += batch_loss.item()
		batch_count += 1
		torch_optim_SGD.step()
	avg_f /= batch_count
	print(f'[', end='')
	print(f'  {time.time()-startTime:9.3f},', end='')
	print(f'  {epoch_index:4d},', end='')
	print(f'  {avg_f:12.6e}', end='')
	print(f'  ]')
print('];')
msgtime()
