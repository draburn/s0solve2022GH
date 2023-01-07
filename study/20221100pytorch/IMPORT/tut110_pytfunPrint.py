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


for epoch in range(2):  # loop over the dataset multiple times

    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):
        
        #print("Calling print_x()...")
        #print_x(net)
        #print("Finished print_x().")
        
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        
        print("Printing loss...")
        print(f"{loss:.17e}")
        print("Finished printing loss. (It's just one non-negative scalar.)")
        exit()
        
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
