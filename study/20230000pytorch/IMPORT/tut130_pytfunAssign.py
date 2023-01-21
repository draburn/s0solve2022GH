# DRaburn 2022-10-24
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image
# "Training a Classifier"

# Please run tut10_dl.py first.

#
# 1. Load and normalize CIFAR10
#

print("Initializing...")

import torch
import torchvision
import torchvision.transforms as transforms

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
print("Initialization complete.")
print("")


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
    print(numel_list)
    index_list = [0];
    n = 0
    for n in range(len(numel_list)):
        index_list.append(index_list[n] + numel_list[n])
    print(index_list)
    return index_list


def get_vec_from_data__append(d):
    sz = d.size()
    num_dim = len(sz)
    v = []
    if (1==num_dim):
        for n0 in range(sz[0]):
            v.append(float(d[n0]))
    elif (2==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                v.append(float(d[n0][n1]))
    elif (3==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    v.append(float(d[n0][n1][n2]))
    elif (4==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        v.append(float(d[n0][n1][n2][n3]))
    else:
        print("EXCEPTION: Unsupported number of dimensions.")
        raise BaseException
    return v


def get_vec_from_data__prealloc(d):
    numel = get_numel_of_data(d)
    #v = range(numel) # This doesn't work. Ugh.
    v = []
    for nv in range(numel):
        v.append(0.0)
    sz = d.size()
    num_dim = len(sz)
    nv = 0
    if (1==num_dim):
        for n0 in range(sz[0]):
            v[nv] = d[n0]
            nv += 1
    elif (2==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                v[nv] = d[n0][n1]
                nv += 1
    elif (3==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    v[nv] = d[n0][n1][n2]
                    nv += 1
    elif (4==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        v[nv] = d[n0][n1][n2][n3]
                        nv += 1
    else:
        print("EXCEPTION: Unsupported number of dimensions.")
        raise BaseException
    return v



def get_vec_from_data(d):
    v = get_vec_from_data__append(d)
    return v

def get_appended_vec_from_data(x,d):
    sz = d.size()
    num_dim = len(sz)
    if (1==num_dim):
        for n0 in range(sz[0]):
            x.append(float(d[n0]))
    elif (2==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                x.append(float(d[n0][n1]))
    elif (3==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    x.append(float(d[n0][n1][n2]))
    elif (4==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        x.append(float(d[n0][n1][n2][n3]))
    else:
        print("EXCEPTION: Unsupported number of dimensions.")
        raise BaseException
    return x

def get_x_from_net(dummy_net):
    x = []
    x = get_appended_vec_from_data(x, dummy_net.conv1.bias.data)
    x = get_appended_vec_from_data(x, dummy_net.conv1.weight.data)
    x = get_appended_vec_from_data(x, dummy_net.conv2.bias.data)
    x = get_appended_vec_from_data(x, dummy_net.conv2.weight.data)
    x = get_appended_vec_from_data(x, dummy_net.fc1.bias.data)
    x = get_appended_vec_from_data(x, dummy_net.fc1.weight.data)
    x = get_appended_vec_from_data(x, dummy_net.fc2.bias.data)
    x = get_appended_vec_from_data(x, dummy_net.fc2.weight.data)
    x = get_appended_vec_from_data(x, dummy_net.fc3.bias.data)
    x = get_appended_vec_from_data(x, dummy_net.fc3.weight.data)
    return x

def assign_x_to_data(x, d):
    sz = d.size()
    num_dim = len(sz)
    nx = 0
    if (1==num_dim):
        for n0 in range(sz[0]):
            d[n0] = x[nx]
            nx += 1
    elif (2==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                d[n0][n1] = x[nx]
                nx += 1
    elif (3==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    d[n0][n1][n2] = x[nx]
                    nx += 1
    elif (4==num_dim):
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        d[n0][n1][n2][n3] = x[nx]
                        nx += 1
    else:
        print("EXCEPTION: Unsupported number of dimensions.")
        raise BaseException


def assign_x_to_net(this_x, this_net):
    index_list = get_index_list(this_net)
    assign_x_to_data(this_x[index_list[0]:index_list[1]], this_net.conv1.bias.data)
    assign_x_to_data(this_x[index_list[1]:index_list[2]], this_net.conv1.weight.data)
    assign_x_to_data(this_x[index_list[2]:index_list[3]], this_net.conv2.bias.data)
    assign_x_to_data(this_x[index_list[3]:index_list[4]], this_net.conv2.weight.data)
    assign_x_to_data(this_x[index_list[4]:index_list[5]], this_net.fc1.bias.data)
    assign_x_to_data(this_x[index_list[5]:index_list[6]], this_net.fc1.weight.data)
    assign_x_to_data(this_x[index_list[6]:index_list[7]], this_net.fc2.bias.data)
    assign_x_to_data(this_x[index_list[7]:index_list[8]], this_net.fc2.weight.data)
    assign_x_to_data(this_x[index_list[8]:index_list[9]], this_net.fc3.bias.data)
    assign_x_to_data(this_x[index_list[9]:index_list[10]], this_net.fc3.weight.data)



vec_x = get_x_from_net(net)

for epoch in range(2):  # loop over the dataset multiple times
    
    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):
        
        for n in range(len(vec_x)):
            vec_x[n] = 0.0
        assign_x_to_net(vec_x,net)
        print_x(net)
        for n in range(len(vec_x)):
            vec_x[n] = 1.0
        assign_x_to_net(vec_x,net)
        print_x(net)
        exit()
        
        #
        #        
        # DO THE ASSIGNMENT THING
        #
        # To (help?) ensure assignment works, we best assign the values directly.
        nx = 0;
        #
        # conv1
        d = net.conv1.bias.data
        sz = d.size()
        num_dim = len(sz)
        if (1==num_dim):
            for n0 in range(sz[0]):
                d[n0] = float(my_x[nx])
                nx += 1
        elif (2==num_dim):
            for n0 in range(sz[0]):
                for n1 in range(sz[1]):
                    d[n0][n1] = my_x[nx]
                    nx += 1
        elif (3==num_dim):
            for n0 in range(sz[0]):
                for n1 in range(sz[1]):
                    for n2 in range(sz[2]):
                        d[n0][n1][n2] = my_x[nx]
                        nx += 1
        elif (4==num_dim):
            for n0 in range(sz[0]):
                for n1 in range(sz[1]):
                    for n2 in range(sz[2]):
                        for n3 in range(sz[3]):
                            d[n0][n1][n2][n3] = my_x[nx]
                            nx += 1
        else:
            print("EXCEPTION: Unsupported number of dimensions.")
            raise BaseException
        #
        exit()
        # conv2
        # ...
        
        #
        # FINISHED THE ASSIGNMENT THING
        #
        #
        
        
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        
        numel_grad = get_numel_grad(net)
        print(f"numel_grad = {numel_grad}")
        print("GOODBYE!")
        exit()
        
        print("Printing loss...")
        print(f"{loss:.17e}")
        print("Finished printing loss. (It's just one non-negative scalar.)")
        
        print("Calling print_grad()...")
        print_grad(net)
        print("Finished print_grad().")
        
        break
        
        optimizer.step()
        
        print("Calling print_x() for *UPDATED* net...")
        print_x(net)
        print("Finished print_x() for *UPDATED* net.")
        
        break
        # print statistics
        running_loss += loss.item()
        if i % 2000 == 1999:    # print every 2000 mini-batches
            print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss / 2000:.3f}')
            running_loss = 0.0
    break;

print('Finished Training Demo')
