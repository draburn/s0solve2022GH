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

def print_net_data(d):
    sz = d.size()
    num_dim = len(sz)
    if (1==num_dim):
        for n0 in range(sz[0]):
            print(f"{d[n0]:.17e}")
    elif (2==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                print(f"{d[n0][n1]:.17e}")
    elif (3==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    print(f"{d[n0][n1][n2]:.17e}")
    elif (4==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        print(f"{d[n0][n1][n2][n3]:.17e}")
    else:
        print("EXCEPTION: Unsupported number of dimensions.")
        raise BaseException

def print_x(dummy_net):
    print_net_data(dummy_net.conv1.bias.data)
    print_net_data(dummy_net.conv1.weight.data)
    print_net_data(dummy_net.conv2.bias.data)
    print_net_data(dummy_net.conv2.weight.data)
    print_net_data(dummy_net.fc1.bias.data)
    print_net_data(dummy_net.fc1.weight.data)
    print_net_data(dummy_net.fc2.bias.data)
    print_net_data(dummy_net.fc2.weight.data)
    print_net_data(dummy_net.fc2.bias.data)
    print_net_data(dummy_net.fc3.weight.data)

def print_grad(dummy_net):
    print_net_data(dummy_net.conv1.bias.grad)
    print_net_data(dummy_net.conv1.weight.grad)
    print_net_data(dummy_net.conv2.bias.grad)
    print_net_data(dummy_net.conv2.weight.grad)
    print_net_data(dummy_net.fc1.bias.grad)
    print_net_data(dummy_net.fc1.weight.grad)
    print_net_data(dummy_net.fc2.bias.grad)
    print_net_data(dummy_net.fc2.weight.grad)
    print_net_data(dummy_net.fc3.bias.grad)
    print_net_data(dummy_net.fc3.weight.grad)


def get_numel_of_data(d):
    sz = d.size()
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

def get_numel_list(dummy_net):
    numel_list = [];
    numel_list.append(get_numel_of_data(dummy_net.conv1.bias.data))
    numel_list.append(get_numel_of_data(dummy_net.conv1.weight.data))
    numel_list.append(get_numel_of_data(dummy_net.conv2.bias.data))
    numel_list.append(get_numel_of_data(dummy_net.conv2.weight.data))
    numel_list.append(get_numel_of_data(dummy_net.fc1.bias.data))
    numel_list.append(get_numel_of_data(dummy_net.fc1.weight.data))
    numel_list.append(get_numel_of_data(dummy_net.fc2.bias.data))
    numel_list.append(get_numel_of_data(dummy_net.fc2.weight.data))
    numel_list.append(get_numel_of_data(dummy_net.fc3.bias.data))
    numel_list.append(get_numel_of_data(dummy_net.fc3.weight.data))
    return numel_list

def get_index_list(dummy_net):
    numel_list = get_numel_list(dummy_net)
    index_list = [0];
    n = 0
    for n in range(len(numel_list)):
        index_list.append(index_list[n] + numel_list[n])
    return index_list

def get_assigned_net_to_x(this_net, this_x):
    index_list = get_index_list(this_net)
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

def get_assigned_net_to_grad(this_net, this_grad):
    index_list = get_index_list(this_net)
    this_grad[index_list[0]:index_list[1]] = numpy.reshape(this_net.conv1.bias.grad.numpy(), -1)
    this_grad[index_list[1]:index_list[2]] = numpy.reshape(this_net.conv1.weight.grad.numpy(), -1)
    this_grad[index_list[2]:index_list[3]] = numpy.reshape(this_net.conv2.bias.grad.numpy(), -1)
    this_grad[index_list[3]:index_list[4]] = numpy.reshape(this_net.conv2.weight.grad.numpy(), -1)
    this_grad[index_list[4]:index_list[5]] = numpy.reshape(this_net.fc1.bias.grad.numpy(), -1)
    this_grad[index_list[5]:index_list[6]] = numpy.reshape(this_net.fc1.weight.grad.numpy(), -1)
    this_grad[index_list[6]:index_list[7]] = numpy.reshape(this_net.fc2.bias.grad.numpy(), -1)
    this_grad[index_list[7]:index_list[8]] = numpy.reshape(this_net.fc2.weight.grad.numpy(), -1)
    this_grad[index_list[8]:index_list[9]] = numpy.reshape(this_net.fc3.bias.grad.numpy(), -1)
    this_grad[index_list[9]:index_list[10]] = numpy.reshape(this_net.fc3.weight.grad.numpy(), -1)
    return this_grad

def get_dat_from_x(this_x, this_index0, this_index1, this_dat ):
    #return torch.from_numpy( numpy.reshape( numpy.asarray(this_x[this_index0:this_index1]), this_dat.size() ) )
    ###return torch.from_numpy( numpy.reshape( numpy.asarray(float(this_x[this_index0:this_index1])), this_dat.size() ) )
    #return torch.from_numpy( numpy.float(numpy.reshape( numpy.asarray(this_x[this_index0:this_index1]), this_dat.size() )) ).
    return torch.Tensor( numpy.reshape( numpy.asarray(this_x[this_index0:this_index1]), this_dat.size() ) )


index_list = get_index_list(net)
sxsolve_x = []
sxsolve_grad = []
sxsolve_gradwmom = []
for n in range(index_list[-1]):
    sxsolve_x.append(1.0)
    sxsolve_grad.append(1.0)
    sxsolve_gradwmom.append(1.0)

if (test_and_quit):
    sxsolve_x = get_assigned_net_to_x(net, sxsolve_x)
    pass
else:
    sxsolve_x = get_assigned_net_to_x(net, sxsolve_x)
    for n in range(len(sxsolve_x)):
        sxsolve_x[n] = sxsolve_x[n] + 0.0
    pass

#optimizer = optim.SGD(net.parameters(), lr=0.001, momentum=0.9)
sxsolve_lr = 0.01
sxsolve_momentum = 0.9
sxsolve_have_initial_grad = 0

tottime_assign_x_to_net = 0.0
tottime_get_assigned_net_to_grad = 0.0

print("Initialization complete.")
print(f"  Elapsed time = {time.time()-start_time}s")
print(f"test_and_quit = {test_and_quit}")
print("Main loop...")
for epoch in range(2):  # loop over the dataset multiple times
    
    running_loss = 0.0
    running_feval_count = 0
    running_time0 = time.time()
    for i, data in enumerate(trainloader, 0):
        
        # ASSIGN X
        time_tic = time.time()
        if (test_and_quit):
            print(f"len(sxsolve_x) = {len(sxsolve_x)}")
            print(f"sxsolve_x[0] = {sxsolve_x[0]:.17e}")
            print(f"sxsolve_x[-1] = {sxsolve_x[-1]:.17e}")
            #
            print("Before 'assign_x_to_net'...")
            print(f"net.conv1.bias.data[0] = {net.conv1.bias.data[0]:.17e}")
            print(f"net.fc1.weight.data[0][0] = {net.fc3.weight.data[0][0]:.17e}")
            print(f"net.fc3.weight.data[-1][-1] = {net.fc3.weight.data[-1][-1]:.17e}")
            print(f"net.conv1.bias.data.size()    = {net.conv1.bias.data.size()}")
            print(f"net.conv1.weight.data.size() = {net.conv1.weight.data.size()}")
            print(f"net.conv2.bias.data.size()    = {net.conv2.bias.data.size()}")
            print(f"net.conv2.weight.data.size() = {net.conv2.weight.data.size()}")
            print(f"net.fc1.bias.data.size()    = {net.fc1.bias.data.size()}")
            print(f"net.fc1.weight.data.size() = {net.fc1.weight.data.size()}")
            print(f"net.fc2.bias.data.size()    = {net.fc2.bias.data.size()}")
            print(f"net.fc2.weight.data.size() = {net.fc2.weight.data.size()}")
            print(f"net.fc3.bias.data.size()    = {net.fc3.bias.data.size()}")
            print(f"net.fc3.weight.data.size() = {net.fc3.weight.data.size()}")
            net.conv1.bias.data   = get_dat_from_x(sxsolve_x, index_list[0], index_list[1], net.conv1.bias.data )
            net.conv1.weight.data = get_dat_from_x(sxsolve_x, index_list[1], index_list[2], net.conv1.weight.data )
            net.conv2.bias.data   = get_dat_from_x(sxsolve_x, index_list[2], index_list[3], net.conv2.bias.data )
            net.conv2.weight.data = get_dat_from_x(sxsolve_x, index_list[3], index_list[4], net.conv2.weight.data )
            net.fc1.bias.data   = get_dat_from_x(sxsolve_x, index_list[4], index_list[5], net.fc1.bias.data )
            net.fc1.weight.data = get_dat_from_x(sxsolve_x, index_list[5], index_list[6], net.fc1.weight.data )
            net.fc2.bias.data   = get_dat_from_x(sxsolve_x, index_list[6], index_list[7], net.fc2.bias.data )
            net.fc2.weight.data = get_dat_from_x(sxsolve_x, index_list[7], index_list[8], net.fc2.weight.data )
            net.fc3.bias.data   = get_dat_from_x(sxsolve_x, index_list[8], index_list[9], net.fc3.bias.data )
            net.fc3.weight.data = get_dat_from_x(sxsolve_x, index_list[9], index_list[10], net.fc3.weight.data )
            print("After 'assign_x_to_net'...")
            print(f"net.conv1.bias.data[0] = {net.conv1.bias.data[0]:.17e}")
            print(f"net.fc1.weight.data[0][0] = {net.fc1.weight.data[0][0]:.17e}")
            print(f"net.fc3.weight.data[-1][-1] = {net.fc3.weight.data[-1][-1]:.17e}")
            print(f"net.fc3.weight.data[-1][-1] = {net.fc3.weight.data[-1][-1]:.17e}")
            print(f"net.conv1.bias.data.size()    = {net.conv1.bias.data.size()}")
            print(f"net.conv1.weight.data.size() = {net.conv1.weight.data.size()}")
            print(f"net.conv2.bias.data.size()    = {net.conv2.bias.data.size()}")
            print(f"net.conv2.weight.data.size() = {net.conv2.weight.data.size()}")
            print(f"net.fc1.bias.data.size()    = {net.fc1.bias.data.size()}")
            print(f"net.fc1.weight.data.size() = {net.fc1.weight.data.size()}")
            print(f"net.fc2.bias.data.size()    = {net.fc2.bias.data.size()}")
            print(f"net.fc2.weight.data.size() = {net.fc2.weight.data.size()}")
            print(f"net.fc3.bias.data.size()    = {net.fc3.bias.data.size()}")
            print(f"net.fc3.weight.data.size() = {net.fc3.weight.data.size()}")
        else:
            ###assign_x_to_net(sxsolve_x, net)
            net.conv1.bias.data   = get_dat_from_x(sxsolve_x, index_list[0], index_list[1], net.conv1.bias.data )
            net.conv1.weight.data = get_dat_from_x(sxsolve_x, index_list[1], index_list[2], net.conv1.weight.data )
            net.conv2.bias.data   = get_dat_from_x(sxsolve_x, index_list[2], index_list[3], net.conv2.bias.data )
            net.conv2.weight.data = get_dat_from_x(sxsolve_x, index_list[3], index_list[4], net.conv2.weight.data )
            net.fc1.bias.data   = get_dat_from_x(sxsolve_x, index_list[4], index_list[5], net.fc1.bias.data )
            net.fc1.weight.data = get_dat_from_x(sxsolve_x, index_list[5], index_list[6], net.fc1.weight.data )
            net.fc2.bias.data   = get_dat_from_x(sxsolve_x, index_list[6], index_list[7], net.fc2.bias.data )
            net.fc2.weight.data = get_dat_from_x(sxsolve_x, index_list[7], index_list[8], net.fc2.weight.data )
            net.fc3.bias.data   = get_dat_from_x(sxsolve_x, index_list[8], index_list[9], net.fc3.bias.data )
            net.fc3.weight.data = get_dat_from_x(sxsolve_x, index_list[9], index_list[10], net.fc3.weight.data )
            pass
        tottime_assign_x_to_net += time.time() - time_tic
        
        
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        
        #
        # GRAB F AND GRAD
        time_tic = time.time()
        if (test_and_quit):
            print(f"loss = {loss:.17e}")
            print(f"net.conv1.bias.grad[0] = {net.conv1.bias.grad[0]:.17e}")
            print(f"net.fc3.weight.grad[-1][-1] = {net.fc3.weight.grad[-1][-1]:.17e}")
            #
            sxsolve_f = loss
            print(f"sxsolve_f = {sxsolve_f:.17e}")
            print(f"len(sxsolve_grad) = {len(sxsolve_grad)}")
            #
            print("Before get_assigned_net_to_grad()...")
            print(f"sxsolve_grad[0] = {sxsolve_grad[0]:.17e}")
            print(f"sxsolve_grad[-1] = {sxsolve_grad[-1]:.17e}")
            sxsolve_grad = get_assigned_net_to_grad(net, sxsolve_grad)
            print("After get_assigned_net_to_grad()...")
            print(f"sxsolve_grad[0] = {sxsolve_grad[0]:.17e}")
            print(f"sxsolve_grad[-1] = {sxsolve_grad[-1]:.17e}")
            exit()
        else:
            sxsolve_f = loss
            sxsolve_grad = get_assigned_net_to_grad(net, sxsolve_grad)
            tottime_get_assigned_net_to_grad += time.time() - time_tic
        
        # TAKE STEP
        use_original_step = False
        if (use_original_step):
            optimizer.step()
        else:
            if (sxsolve_have_initial_grad):
                for n in range(len(sxsolve_grad)):
                    sxsolve_gradwmom[n] = sxsolve_momentum*sxsolve_gradwmom[n] + (1.0-sxsolve_momentum)*sxsolve_grad[n]
            else:
                sxsolve_have_initial_grad = True
                for n in range(len(sxsolve_grad)):
                    sxsolve_gradwmom[n] = sxsolve_grad[n]
            for n in range(len(sxsolve_x)):
                sxsolve_x[n] -= sxsolve_lr*sxsolve_gradwmom[n]
        
        # print statistics
        running_loss += loss.item()
        running_feval_count += 1
        #if i % 2000 == 1999:    # print every 2000 mini-batches
        #print_interval = 2000
        print_interval = 200
        if i % print_interval == (print_interval-1):    # print every few mini-batches
        #if (time.time()-running_time0>3.0):
            print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss * 1.0 / running_feval_count:.3f}')
            print(f"  Avg time per feval = { (time.time()-running_time0)*1.0 / running_feval_count }")
            print(f"  Elapsed time = {time.time()-start_time}s")
            print(f"  tottime_assign_x_to_net = {tottime_assign_x_to_net}s")
            print(f"  tottime_get_assigned_net_to_grad = {tottime_get_assigned_net_to_grad}s")
            running_loss = 0.0
            running_time0 = time.time()
            running_feval_count = 0

print('Finished Training Demo')
