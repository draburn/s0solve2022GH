# DRaburn 2022-10-24
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image
# "Training a Classifier"

# Please run tut10_dl.py first.

#
# 1. Load and normalize CIFAR10
#

import torch
import torchvision
import torchvision.transforms as transforms

import time
import inspect
import numpy
start_time = time.time()

# Init logging.
frame = inspect.currentframe()
def msg( *arguments, **keywords ):
    print( f'[{__file__}.{frame.f_lineno:05d}]', *arguments, **keywords )
msg( f'time.asctime() = {time.asctime()}' )

torch.manual_seed(0)

transform = transforms.Compose([transforms.ToTensor(), transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])

num_eptots = 1
batch_size = 50000
lr = 1.0
max_iter = 5000
max_eval = round(max_iter * 1.25)
tolerance_grad = 1.0E-5
tolerance_change = 1.0E-9
history_size = 100
line_search_fn = None
line_search_fn = 'strong_wolfe'
fname_x_initial = ''
#fname_x_initial = 'in_vecX.np'

msg( f'num_eptots = {num_eptots}' )
msg( f'batch_size = {batch_size}' )
msg( f'lr = {lr:0.14E}' )
msg( f'max_iter = {max_iter}' )
msg( f'max_eval = {max_eval}' )
msg( f'tolerance_grad = {tolerance_grad:0.14E}' )
msg( f'tolerance_change = {tolerance_change:0.14E}' )
msg( f'line_search_fn = \'{line_search_fn}\'' )
msg( f'history_size = {history_size}' )
msg( f'fname_x_initial = \'{fname_x_initial}\'' )

trainset = torchvision.datasets.CIFAR10(root='./data', train=True, download=False, transform=transform)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size, shuffle=True, num_workers=2)
classes = ('plane', 'car', 'bird', 'cat', 'deer', 'dog', 'frog', 'horse', 'ship', 'truck')

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
#optimizer = optim.LBFGS(net.parameters())
optimizer = optim.LBFGS(
  net.parameters(),
  lr=lr,
  max_iter=max_iter,
  max_eval=max_eval,
  tolerance_grad=tolerance_grad,
  tolerance_change=tolerance_change,
  history_size=history_size,
  line_search_fn=line_search_fn )

msg( 'Begin my hacks for injection/extraction...' )
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

size_list = get_size_list(net)
index_list = get_index_list(net)
sxsolve_x = init_x_from_net(net)
sizeX = index_list[-1]
msg( f'sizeX = {sizeX}' )

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

if (''!=fname_x_initial):
    msg( 'Reading from x disk...' )
    foo = numpy.fromfile(fname_x_initial,dtype=numpy.float32)
    msg( f'foo.size = {foo.size}' )
    sxsolve_x[:] = foo[:]

msg( 'Begin training...' )
print( '' )
print( '[' )
for eptot in range(num_eptots):
    for i, data in enumerate(trainloader, 0):
        inputs, labels = data
        def closure():
            optimizer.zero_grad()
            outputs = net(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            print( f'[ {time.time()-start_time:8.2f}, {eptot:4d}, {i:5d}, {loss:0.14e}]' )
            return loss
        optimizer.step(closure)
        #outputs = net(inputs)
        #loss = criterion(outputs, labels)
        #print( f'[ {time.time()-start_time:8.2f}, {eptot:4d}, {i:5d}, {loss:0.14e}]' )
print( '];' )
print( '' )
msg( 'Training complete.' )

# Let's quickly save...
#PATH = './temp.pth'
#torch.save(net.state_dict(), PATH)
#import numpy
#import matplotlib.pyplot as plt
