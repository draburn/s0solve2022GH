myclear;
setprngstates(90855792);
%
sizeX = 2;
sizeF = 2;
numVals = 10;
%
dat.c0.ary = randn(sizeF,1);
dat.c(1).ary = randn(sizeF,sizeX);
dat.c(2).ary = randn(sizeF,sizeX,sizeX);
dat.c(3).ary = randn(sizeF,sizeX,sizeX,sizeX);
dat.ord = 3;
%
funchF = @(x)(studyfunc0510_eval(x,dat));
%
assert(2==sizeX);
funchFNorm = @(x)( sqrt(sum(funchF(x).^2,1)) );
%funchZ = @(x,y)( funchFNorm([x;y]) );
funchZ = @(x,y)( asinh(funchFNorm([x;y])*1)/1 );
%contourfunch(funchZ);
[ zMin, zMax, zAvg, zVar ] = contourfunch(funchZ,[-1.5,2,-1,1.5],51,51,20,true,true)
axis equal;
colormap(0.5+0.5*jet(1000));
grid on;
%
matX = randn(sizeX,numVals);
matF = funchF(matX);
