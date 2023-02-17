# DRaburn 2023-02-16
# Based on "Training a Classifier":
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image

# Import.
import inspect
import time
import torch
import torchvision
import numpy as np

# Startup.
start_time = time.time()
main_frame = inspect.currentframe()
def msg(*arguments, **keywords):
	print(f'[{__file__}.{main_frame.f_lineno:05d}]', *arguments, **keywords)
def msgtime():
	print(f'[{__file__}.{main_frame.f_lineno:05d}] It is currently {time.asctime()}, {time.time()-start_time:0.3f}s elapsed.')
msgtime()
msg(f'Initializing...')

# Set base params.
torch_seed = 0
CIFAR10_root = '../../dat/CIFAR10data'
batch_size = 500
learning_rate = 0.1
momentum_coefficient = 0.9
max_num_epochs = 10
msg(f'torch_seed = {torch_seed}')
msg(f'batch_size = {batch_size}')
msg(f'CIFAR10_root = "{CIFAR10_root}"')
msg(f'learning_rate = {learning_rate:0.9E}')
msg(f'momentum_coefficient = {momentum_coefficient:0.9E}')
msg(f'max_num_epochs = {max_num_epochs}')

# Init Torch, NN, etc.
torch.manual_seed(torch_seed)
transform = torchvision.transforms.Compose([
  torchvision.transforms.ToTensor(),
  torchvision.transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])
trainset = torchvision.datasets.CIFAR10(root=CIFAR10_root, train=True, download=False, transform=transform)
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
net = Net()
torch_optim_SGD = torch.optim.SGD(net.parameters(), lr=learning_rate, momentum=momentum_coefficient)
loss_criterion = torch.nn.CrossEntropyLoss()

# Main loop.
msg('Finished initialization.')
msgtime()
msg('Starting main loop...')
print('')
print('[')
for epoch in range(max_num_epochs):
	running_loss = 0.0
	running_batch_count = 0
	for batch_index, batch_data in enumerate(trainloader, 0):
		# Prep.
		torch_optim_SGD.zero_grad()
		
		# Calculate f (loss) and gradient.
		batch_inputs, batch_labels = batch_data
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		
		# Grab some data.
		running_loss += batch_loss.item()
		running_batch_count += 1
		
		# Take step and cleanup.
		torch_optim_SGD.step()
	# End batch loop.
	print(f'[', end='')
	print(f' {time.time()-start_time:10.3f} {epoch:5d}', end='')
	print(f'  ', end='')
	print(f'  {running_loss / running_batch_count:15.9E}', end='')
	print(f' ]')
print(']')
print('')
msg('Finished main loop.')

# Exit.
msgtime()
msg( "Goodbye!" )
