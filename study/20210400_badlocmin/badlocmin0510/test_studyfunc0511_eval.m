myclear;
%setprngstates(95180704);
%
sizeX = 2;
sizeF = 2;
%
dat.numerDat.c0.ary = randn(sizeF,1);
dat.numerDat.c(1).ary = randn(sizeF,sizeX);
dat.numerDat.c(2).ary = randn(sizeF,sizeX,sizeX)/2;
dat.numerDat.c(3).ary = randn(sizeF,sizeX,sizeX,sizeX)/6;
dat.numerDat.c(4).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX)/24;
dat.numerDat.c(5).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX)/120;
dat.numerDat.c(6).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX,sizeX)/720;
dat.numerDat.c(7).ary = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX,sizeX,sizeX)/5040;
dat.numerDat.ord = 7;
denomScale = 1E-3;
dat.denomDat.c0.ary = randn(sizeF,1)*sqrt(denomScale);
dat.denomDat.c(1).ary = randn(sizeF,sizeX)*sqrt(denomScale);
dat.denomDat.c(2).ary = randn(sizeF,sizeX,sizeX)*(sqrt(denomScale)/2);
dat.denomDat.c(3).ary = randn(sizeF,sizeX,sizeX,sizeX)*(sqrt(denomScale)/6);
dat.denomDat.ord = 2;
dat.vecXRoot = randn(sizeX,1);
%
funchF = @(x)(studyfunc0511_eval(x,dat));
%
assert(2==sizeX);
funchFNorm = @(x)( sqrt(sum(funchF(x).^2,1)) );
%funchZ = @(x,y)( funchFNorm([x;y]) );
vizC = 1E0;
funchZ = @(x,y)( asinh(funchFNorm([x;y])*vizC)/vizC );
[ zMin, zMax, zAvg, zVar ] = contourfunch(funchZ,5*[-1,1,-1,1])
axis equal;
colormap(0.5+0.5*jet(1000));
grid on;
%
numVals = 10;
matX = randn(sizeX,numVals);
matF = funchF(matX);
