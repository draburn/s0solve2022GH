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
