myclear;
%setprngstates(90855792); %Use with 2x2, dat.ord = 3 to get 3 loc min.
%
sizeX = 2;
sizeF = 2;
numVals = 10;
%
dat.c0.ary = randn(sizeF,1);
dat.c(1).ary = randn(sizeF,sizeX);
dat.c(2).ary = randn(sizeF,sizeX,sizeX)/2;
dat.c(3).ary = randn(sizeF,sizeX,sizeX,sizeX)/6;
dat.c(4).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX)/24;
dat.c(5).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX)/120;
dat.c(6).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX,sizeX)/720;
dat.c(7).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX,sizeX,sizeX)/5040;
dat.ord = 7;
%
funchF = @(x)(studyfunc0510_eval(x,dat));
%
assert(2==sizeX);
funchFNorm = @(x)( sqrt(sum(funchF(x).^2,1)) );
%funchZ = @(x,y)( funchFNorm([x;y]) );
funchZ = @(x,y)( asinh(funchFNorm([x;y])*1)/1 );
contourfunch(funchZ);
axis equal;
colormap(0.5+0.5*jet(1000));
grid on;
%
matX = randn(sizeX,numVals);
matF = funchF(matX);
