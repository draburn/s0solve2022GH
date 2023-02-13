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
start_time = time.time()
torch.manual_seed(0)

transform = transforms.Compose([transforms.ToTensor(), transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])


num_eptots = 5
batch_size = 50000
lr = 1.0
max_iter = 20
max_eval = round(max_iter * 1.25)
tolerance_grad = 1.0E-5
tolerance_change = 1.0E-9
history_size = 100
line_search_fn = None
print( f'num_eptots = {num_eptots}' )
print( f'batch_size = {batch_size}' )
print( f'lr = {lr:0.14E}' )
print( f'max_iter = {max_iter}' )
print( f'max_eval = {max_eval}' )
print( f'tolerance_grad = {tolerance_grad:0.14E}' )
print( f'tolerance_change = {tolerance_change:0.14E}' )
print( f'line_search_fn = {line_search_fn}' )
print( f'history_size = {history_size}' )

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

print( 'Begin training...' )
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
print( ']' )
print( 'Training complete.' )

# Let's quickly save...
#PATH = './temp.pth'
#torch.save(net.state_dict(), PATH)
#import numpy
#import matplotlib.pyplot as plt
