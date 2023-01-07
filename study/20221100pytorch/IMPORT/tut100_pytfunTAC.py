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

for epoch in range(2):  # loop over the dataset multiple times

    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):
        print("Performing one feval...")
        
        
        sz = net.conv1.weight.data.size()
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        net.conv1.weight.data[n0][n1][n2][n3] = 0.0
                        
        
        print("")
        print("vvvvvvvvvvvvvvvvvvvv Weights vvvvvvvvvvvvvvvvvvvv")
        print("Weights: Conv1...")
        print(net.conv1.weight.data.size())
        d = net.conv1.weight.data
        sz = d.size()
        for n0 in range(sz[0]):
            for n1 in range(sz[1]):
                for n2 in range(sz[2]):
                    for n3 in range(sz[3]):
                        print(f"{d[n0][n1][n2][n3]:.16e}")
        exit()
        print(net.conv1.weight.data)
        print("")
        print("Aborting...")
        exit()
        #print("Weights: Pool...")
        #print(net.pool.weight.data.size())
        #print(net.pool.weight.data)
        #print(""
        print("Weights: Conv2...")
        print(net.conv2.weight.data.size())
        print(net.conv2.weight.data)
        print("")
        print("Weights: FC1...")
        print(net.fc1.weight.data.size())
        print(net.fc1.weight.data)
        print("")
        print("Weights: FC2...")
        print(net.fc2.weight.data.size())
        print(net.fc2.weight.data)
        print("")
        print("Weights: FC3...")
        print(net.fc3.weight.data.size())
        print(net.fc3.weight.data)
        print("")
        print("^^^^^^^^^^^^^^^^^^^^ Weights ^^^^^^^^^^^^^^^^^^^^")
        print("")
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        
        print("")
        print("vvvvvvvvvvvvvvvvvvvv Loss vvvvvvvvvvvvvvvvvvvv")
        print(loss)
        print("^^^^^^^^^^^^^^^^^^^^ Loss ^^^^^^^^^^^^^^^^^^^^")
        print("")
        
        print("")
        print("vvvvvvvvvvvvvvvvvvvv Gradients vvvvvvvvvvvvvvvvvvvv")
        print("Gradients: Conv1...")
        print(net.conv1.weight.grad.size())
        print(net.conv1.weight.grad)
        print("")
        #print("Gradients: Pool...")
        #print(net.pool.weight.grad)
        #print("")
        print("Gradients: Conv2...")
        print(net.conv2.weight.grad.size())
        print(net.conv2.weight.grad)
        print("")
        print("Gradients: FC1...")
        print(net.fc1.weight.grad.size())
        print(net.fc1.weight.grad)
        print("")
        print("Gradients: FC2...")
        print(net.fc2.weight.grad.size())
        print(net.fc2.weight.grad)
        print("")
        print("Gradients: FC3...")
        print(net.fc3.weight.grad.size())
        print(net.fc3.weight.grad)
        print("")
        print("^^^^^^^^^^^^^^^^^^^^ Gradients ^^^^^^^^^^^^^^^^^^^^")
        print("")
        
        
        print("")
        print("  Feval complete.");
        print("")
        
        # break before we step
        optimizer.step()
        
        print("")
        print("vvvvvvvvvvvvvvvvvvvv UPDATED Weights vvvvvvvvvvvvvvvvvvvv")
        print("UPDATED Weights: Conv1...")
        print(net.conv1.weight.data.size())
        print(net.conv1.weight.data)
        print("")
        #print("Weights: Pool...")
        #print(net.pool.weight.data)
        #print(""
        print("UPDATED Weights: Conv2...")
        print(net.conv2.weight.data.size())
        print(net.conv2.weight.data)
        print("")
        print("UPDATED Weights: FC1...")
        print(net.fc1.weight.data.size())
        print(net.fc1.weight.data)
        print("")
        print("UPDATED Weights: FC2...")
        print(net.fc2.weight.data.size())
        print(net.fc2.weight.data)
        print("")
        print("UPDATED Weights: FC3...")
        print(net.fc3.weight.data.size())
        print(net.fc3.weight.data)
        print("")
        print("^^^^^^^^^^^^^^^^^^^^ UPDATED Weights ^^^^^^^^^^^^^^^^^^^^")
        print("")
        
        
        break;

        # print statistics
        running_loss += loss.item()
        if i % 2000 == 1999:    # print every 2000 mini-batches
            print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss / 2000:.3f}')
            running_loss = 0.0
    break;

print('Finished Training Demo')
