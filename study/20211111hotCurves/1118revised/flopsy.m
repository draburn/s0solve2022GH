tic();
x = 0.9
z = 500
for k=1:z
for m=1:z
for n=1:z
	x = 3.9*x*(1.0-x);
end
end
end
echo__x = x
toc()
