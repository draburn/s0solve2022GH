# DRaburn 2023-02-16
# Based on "Training a Classifier":
# https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html?highlight=display%20image

# Import.
import inspect
import time
import torch
import torchvision
import numpy as np
import danutil
import qnj

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
learning_rate = 0.01
momentum_coefficient = 0.9
max_num_records = 200
#tr_accel_coeff = 2.0
max_num_epochs = 200
fname_x0 = ''
fname_p0 = ''
#fname_x0 = 'in_vecX0.np'
#fname_p0 = 'in_vecP0.np'
dtype_x0 = np.float32
dtype_p0 = np.float32
msg(f'torch_seed = {torch_seed}')
msg(f'batch_size = {batch_size}')
msg(f'CIFAR10_root = "{CIFAR10_root}"')
msg(f'learning_rate = {learning_rate:0.9E}')
msg(f'momentum_coefficient = {momentum_coefficient:0.9E}')
msg(f'max_num_records = {max_num_records}')
#msg(f'tr_accel_coeff = {tr_accel_coeff:0.9E}')
msg(f'max_num_epochs = {max_num_epochs}')
msg(f'fname_x0 = "{fname_x0}"')
msg(f'fname_p0 = "{fname_p0}"')

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

# Stuff for injection/extraction of Torch data.
# DRaburn 2023-02-16:
#  I'd imagine a proper Pytorch API exists, but I couldn't find information about it.
#  A "brute-force comprehensive Python" approach is likely also possible, but, mise;
#   this code is limited to the specific NN defined above.
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
	this_x = np.zeros(cumel_list[-1], dtype=np.float32)
	this_x[cumel_list[0]:cumel_list[1]] = np.reshape(this_net.conv1.bias.data.numpy(), -1)
	this_x[cumel_list[1]:cumel_list[2]] = np.reshape(this_net.conv1.weight.data.numpy(), -1)
	this_x[cumel_list[2]:cumel_list[3]] = np.reshape(this_net.conv2.bias.data.numpy(), -1)
	this_x[cumel_list[3]:cumel_list[4]] = np.reshape(this_net.conv2.weight.data.numpy(), -1)
	this_x[cumel_list[4]:cumel_list[5]] = np.reshape(this_net.fc1.bias.data.numpy(), -1)
	this_x[cumel_list[5]:cumel_list[6]] = np.reshape(this_net.fc1.weight.data.numpy(), -1)
	this_x[cumel_list[6]:cumel_list[7]] = np.reshape(this_net.fc2.bias.data.numpy(), -1)
	this_x[cumel_list[7]:cumel_list[8]] = np.reshape(this_net.fc2.weight.data.numpy(), -1)
	this_x[cumel_list[8]:cumel_list[9]] = np.reshape(this_net.fc3.bias.data.numpy(), -1)
	this_x[cumel_list[9]:cumel_list[10]] = np.reshape(this_net.fc3.weight.data.numpy(), -1)
	return this_x
def init_g_from_net(this_net):
	cumel_list = get_cumel_list(this_net)
	this_grad = np.zeros(cumel_list[-1], dtype=np.float32)
	this_grad[cumel_list[0]:cumel_list[1]] = np.reshape(this_net.conv1.bias.grad.numpy(), -1)
	this_grad[cumel_list[1]:cumel_list[2]] = np.reshape(this_net.conv1.weight.grad.numpy(), -1)
	this_grad[cumel_list[2]:cumel_list[3]] = np.reshape(this_net.conv2.bias.grad.numpy(), -1)
	this_grad[cumel_list[3]:cumel_list[4]] = np.reshape(this_net.conv2.weight.grad.numpy(), -1)
	this_grad[cumel_list[4]:cumel_list[5]] = np.reshape(this_net.fc1.bias.grad.numpy(), -1)
	this_grad[cumel_list[5]:cumel_list[6]] = np.reshape(this_net.fc1.weight.grad.numpy(), -1)
	this_grad[cumel_list[6]:cumel_list[7]] = np.reshape(this_net.fc2.bias.grad.numpy(), -1)
	this_grad[cumel_list[7]:cumel_list[8]] = np.reshape(this_net.fc2.weight.grad.numpy(), -1)
	this_grad[cumel_list[8]:cumel_list[9]] = np.reshape(this_net.fc3.bias.grad.numpy(), -1)
	this_grad[cumel_list[9]:cumel_list[10]] = np.reshape(this_net.fc3.weight.grad.numpy(), -1)
	return this_grad
size_list = get_size_list(net)
cumel_list = get_cumel_list(net)
num_unknowns = cumel_list[-1]
#msg(f'size_list = {size_list}')
#msg(f'cumel_list = {cumel_list}')
msg(f'num_unknowns = {num_unknowns}')

# Initialize our "x".
# DRaburn 2023-02-16:
#  There may be a better way to do this, but, I don't know what it is.
#  First, we'll create a "shared_vecX" and populate it.
#  Then, we make the weights (and biases) in the net point to our "shared_vecX".
shared_vecX = init_x_from_net(net)
net.conv1.bias.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[0]:cumel_list[1]],size_list[0]))
net.conv1.weight.data = torch.from_numpy(np.reshape(shared_vecX[cumel_list[1]:cumel_list[2]],size_list[1]))
net.conv2.bias.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[2]:cumel_list[3]],size_list[2]))
net.conv2.weight.data = torch.from_numpy(np.reshape(shared_vecX[cumel_list[3]:cumel_list[4]],size_list[3]))
net.fc1.bias.data     = torch.from_numpy(np.reshape(shared_vecX[cumel_list[4]:cumel_list[5]],size_list[4]))
net.fc1.weight.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[5]:cumel_list[6]],size_list[5]))
net.fc2.bias.data     = torch.from_numpy(np.reshape(shared_vecX[cumel_list[6]:cumel_list[7]],size_list[6]))
net.fc2.weight.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[7]:cumel_list[8]],size_list[7]))
net.fc3.bias.data     = torch.from_numpy(np.reshape(shared_vecX[cumel_list[8]:cumel_list[9]],size_list[8]))
net.fc3.weight.data   = torch.from_numpy(np.reshape(shared_vecX[cumel_list[9]:cumel_list[10]],size_list[9]))
# I'm not sure what happens to original memory in net, and I don't particularly care.
if (''!=fname_x0):
	msg(f'Reading x0 from disk using dtype = "{dtype_x0}".')
	vecX0 = np.fromfile(fname_x0,dtype=dtype_x0)
	shared_vecX[:] = vecX0[:]
else:
	vecX0 = shared_vecX.copy()
msg(f'vecX0[0] = {vecX0[0]:0.18E}')
msg(f'||vecX0|| = {np.linalg.norm(vecX0):0.18E}')

# Initialize our gradient.
# DRaburn 2023-02-16:
#  There may be a better way to do this, but, I don't know what it is.
#  This approach is basically the same as for shared_vecX, but,
#   it seems we need to do at least one gradient calculation first.
torch_optim_SGD.zero_grad()
for batch_index, batch_data in enumerate(trainloader, 0):
	batch_inputs, batch_labels = batch_data
	batch_outputs = net(batch_inputs)
	batch_loss = loss_criterion(batch_outputs, batch_labels)
	batch_loss.backward()
	break
shared_vecG = init_g_from_net(net)
net.conv1.bias.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[0]:cumel_list[1]],size_list[0]))
net.conv1.weight.grad = torch.from_numpy(np.reshape(shared_vecG[cumel_list[1]:cumel_list[2]],size_list[1]))
net.conv2.bias.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[2]:cumel_list[3]],size_list[2]))
net.conv2.weight.grad = torch.from_numpy(np.reshape(shared_vecG[cumel_list[3]:cumel_list[4]],size_list[3]))
net.fc1.bias.grad     = torch.from_numpy(np.reshape(shared_vecG[cumel_list[4]:cumel_list[5]],size_list[4]))
net.fc1.weight.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[5]:cumel_list[6]],size_list[5]))
net.fc2.bias.grad     = torch.from_numpy(np.reshape(shared_vecG[cumel_list[6]:cumel_list[7]],size_list[6]))
net.fc2.weight.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[7]:cumel_list[8]],size_list[7]))
net.fc3.bias.grad     = torch.from_numpy(np.reshape(shared_vecG[cumel_list[8]:cumel_list[9]],size_list[8]))
net.fc3.weight.grad   = torch.from_numpy(np.reshape(shared_vecG[cumel_list[9]:cumel_list[10]],size_list[9]))

# Initialize momentum.
# This one is easy, since we don't need to interface with Torch.
vecP = np.zeros(num_unknowns, dtype=np.float32)
if (''!=fname_p0):
	msg(f'Reading p0 from disk using dtype = "{dtype_p0}".')
	vecP0 = np.fromfile(fname_p0,dtype=dtype_p0)
	vecP[:] = vecP0[:]
else:
	vecP0 = vecP.copy()
msg(f'vecP0[0] = {vecP0[0]:0.18E}')
msg(f'||vecP0|| = {np.linalg.norm(vecP0):0.18E}')

# Initialize records.
record_matX = np.zeros(( num_unknowns, max_num_records ))
record_matG = np.zeros(( num_unknowns, max_num_records ))
record_rvcF = np.zeros(( 1, max_num_records ))
num_records = 0

# Initialize inter-interval-info and QNJ param.
prev_f = -1.0 # Negative indicates no prev.
prev_fPred = -1.0
tr_size = 0.0
qnj_prm = qnj.prm()

# Main loop.
msg('Finished initialization.')
msgtime()
msg('Starting main loop...')
print('')
print('[')
for epoch_index in range(max_num_epochs):
	vecXSeed = shared_vecX.copy()
	vecPSeed = vecP.copy()
	batch_count = 0
	avg_f = 0.0
	avg_vecX = np.zeros(num_unknowns)
	avg_vecG = np.zeros(num_unknowns)
	avg_fSq = 0.0
	avg_vecXSq = np.zeros(num_unknowns)
	avg_vecGSq = np.zeros(num_unknowns)
	for batch_index, batch_data in enumerate(trainloader, 0):
		# Calculate f (loss) and gradient.
		torch_optim_SGD.zero_grad()
		batch_inputs, batch_labels = batch_data
		batch_outputs = net(batch_inputs)
		batch_loss = loss_criterion(batch_outputs, batch_labels)
		batch_loss.backward()
		batch_f = batch_loss.item()
		
		# Grab data.
		batch_count += 1
		avg_f += batch_f
		avg_vecX[:] += shared_vecX[:]
		avg_vecG[:] += shared_vecG[:]
		avg_fSq += (batch_f**2)
		avg_vecXSq[:] += (shared_vecX[:]**2)
		avg_vecGSq[:] += (shared_vecG[:]**2)
		
		# Take step.
		vecP[:] = (momentum_coefficient * vecP[:]) - (learning_rate * shared_vecG[:])
		shared_vecX[:] += vecP[:]
	# End batch loop.
	
	# Calc stuff over epoch.
	vecXHarvest = shared_vecX.copy()
	vecPHarvest = vecP.copy()
	avg_f /= batch_count
	avg_vecX[:] /= batch_count
	avg_vecG[:] /= batch_count
	avg_fSq /= batch_count
	avg_vecXSq[:] /= batch_count
	avg_vecGSq[:] /= batch_count
	var_f = danutil.var( avg_f, avg_fSq )
	var_x = danutil.var( avg_vecX, avg_vecXSq )
	avg_d = np.linalg.norm( avg_vecX - vecX0 )
	avg_g = np.linalg.norm( avg_vecG )
	var_g = danutil.var( avg_vecG, avg_vecGSq )
	hmsx = np.linalg.norm( vecXHarvest - vecXSeed )
	hmsp = np.linalg.norm( vecPHarvest - vecPSeed )
	
	# Add record.
	# DRaburn 2023-02-16: Always rolling may be wasteful, but POITROME.
	record_matX = np.roll( record_matX, 1 )
	record_matG = np.roll( record_matG, 1 )
	record_rvcF = np.roll( record_rvcF, 1 )
	# DRaburn 2023-02-23:
	# I had intended for roll() to do the following:
	#record_matX[:,1:] = record_matX[:,:-1]
	#record_matG[:,1:] = record_matG[:,:-1]
	#record_rvcF[0,1:] = record_rvcF[0,:-1]
	# What it actually does is different, but the incorrect part gets overwritten anyway,
	# so it doesn't impact the calculations.
	if ( num_records < max_num_records ):
		num_records += 1
	record_matX[:,0] = avg_vecX[:]
	record_matG[:,0] = avg_vecG[:]
	record_rvcF[0,0] = avg_f
	
	# Update trust region.
	tr_lo = 0.1*hmsx
	#tr_hi = 10.0*hmsx
	tr_hi = 100.0*hmsx
	#tr_decel = 0.0
	tr_decel = 0.5
	tr_accel = 1.2
	if (prev_f < 0.0):
		msg(f'Useless(?) TR: {tr_size} {tr_lo}, {tr_hi}')
		tr_size = tr_lo
	elif (prev_fPred < 0.0):
		msg(f'(Re)init TR: {tr_size} {tr_lo}, {tr_hi}')
		tr_size = tr_lo
	elif (avg_f > prev_f):
		msg(f'Decel TR: {tr_size} {tr_lo}, {tr_hi}')
		tr_size *= tr_decel
	#elif (avg_f < (prev_f+prev_fPred)/2.0):
	elif (avg_f < prev_fPred):
		msg(f'Accel TR: {tr_size} {tr_lo}, {tr_hi}')
		# NOTE: "fPred" does nto consider that we start lower and go lower after landing.
		tr_size *= tr_accel
	tr_size = np.median([ tr_size, tr_lo, tr_hi ])
	# End update trust region.
	
	# Calculate jump.
	qnjCode, vecXNext, vecPNext, sizeK, gammaRat, fPred = qnj.calcJump(
	  vecXHarvest,
	  vecPHarvest,
	  record_matX[:,0:num_records],
	  record_matG[:,0:num_records],
	  record_rvcF[:,0:num_records],
	  tr_size,
	  qnj_prm )
	# End calculate jump.
	
	# Report.
	print(f'[', end='')
	print(f' {time.time()-start_time:10.3f} {epoch_index+1:5d}', end='')
	print(f' ', end='')
	print(f'  {avg_f:12.6E} {var_f:12.6E}', end='')
	print(f' ', end='')
	print(f'  {avg_d:12.6E} {var_x:12.6E}', end='')
	print(f'  {avg_g:12.6E} {var_g:12.6E}', end='')
	print(f'  ', end='')
	print(f' {hmsx:12.6E} {np.linalg.norm(vecPSeed):12.6E} {np.linalg.norm(vecPHarvest):12.6E} {tr_size:12.6E}', end='')
	print(f'  ', end='')
	print(f' {np.linalg.norm(vecXNext-vecXHarvest):12.6E} {np.linalg.norm(vecPNext):12.6E}', end='')
	print(f' {qnjCode:3d} {sizeK:3d} {gammaRat:12.6E} {fPred:13.6E}', end='' )
	print(f' ]')
	
	# Save progress
	if ( (epoch_index+1)%100 == 0 ):
		shared_vecX.tofile(f'out_vecX_{epoch_index+1:06d}.np')
		vecP.tofile(f'out_vecP_{epoch_index+1:06d}.np')
	
	# Move to next step.
	prev_f = avg_f
	prev_fPred = fPred
	shared_vecX[:] = vecXNext[:]
	vecP[:] = vecPNext[:]
# End epoch loop.

print('];')
print('')
msg('Finished main loop.')

# Exit.
msgtime()
msg( "Goodbye!" )
