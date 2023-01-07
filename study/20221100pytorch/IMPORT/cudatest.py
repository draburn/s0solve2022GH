import torch
x = torch.rand(5, 3)
print(x)

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")
