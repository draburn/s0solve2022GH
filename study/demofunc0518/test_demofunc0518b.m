myclear;
%setprngstates(11835232);
setprngstates();
sizeX = 2;
sizeF = 2;
funcPrm.vecF0 = randn(sizeF,1);
funcPrm.matJ0 = randn(sizeF,sizeX);
funcPrm.polyOrder = 5;
funcPrm.c2 = randn(sizeF,sizeX,sizeX)/20;
funcPrm.c3 = randn(sizeF,sizeX,sizeX,sizeX)/100;
funcPrm.c4 = randn(sizeF,sizeX,sizeX,sizeX,sizeX)/500;
funcPrm.c5 = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX)/2000;
ax = 4*[-5,5,-5,5];
nx = 101;
ny = 101;
nc = 30;
%
funcPrm.expPrm.numTerms = 0*round(exp(2+randn));
for n=1:funcPrm.expPrm.numTerms
	funcPrm.expPrm.matC(:,n) = randn(sizeF,1);
	funcPrm.expPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.expPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecXScale(1,n) = abs(randn);
end
%
funcPrm.vecF0Denom = ones(2,1);
funcPrm.denomExpPrm.numTerms = 0*round(exp(2+randn));
for n=1:funcPrm.denomExpPrm.numTerms
	funcPrm.denomExpPrm.matC(:,n) = abs(randn(sizeF,1));
	funcPrm.denomExpPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.denomExpPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.denomExpPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.denomExpPrm.rvecXScale(1,n) = abs(randn);
end
%
vecXRoot = zeros(2,1);
vecFPreRoot = demofunc0518_eval(vecXRoot,funcPrm);
%
funchF = @(x)( demofunc0518_eval(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
%funchF([1;1])
%
funchF1XY = @(x,y)( funchF([x;y])(1,:) );
funchF2XY = @(x,y)( funchF([x;y])(2,:) );
funchFNXY = @(x,y)( sqrt(sum(funchF([x;y]).^2,1)) );
%
vizScl = 10.0;
vizFunch = @(x)(abs(asinh(x*vizScl)/vizScl));
%vizFunch = @(x)( double(abs(x)<0.05) );
cMap = 0.6 +0.4*jet;
cMap(1,:) *= 0.4;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
%ax = [-1.569,-1.568,-2.785,-2.784];
%ax = [3.3727,3.374,-0.7037,-0.7034]
%ax = [1.8,3,3.7,4.5];
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchF1XY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc );
colormap(cMap);
axis equal;
grid on;
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchF2XY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc );
colormap(cMap);
axis equal;
grid on;
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchFNXY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc )
colormap(cMap);
axis equal;
grid on;
%
toc
