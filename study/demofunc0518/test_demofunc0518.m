myclear;
%setprngstates();
sizeX = 2;
sizeF = 2;
funcPrm.vecF0 = zeros(sizeF,1);
funcPrm.matJ0 = eye(sizeF,sizeX);
funcPrm.polyOrder = 5;
funcPrm.c2 = randn(sizeF,sizeX,sizeX)/2;
funcPrm.c3 = randn(sizeF,sizeX,sizeX,sizeX)/6;
funcPrm.c4 = randn(sizeF,sizeX,sizeX,sizeX,sizeX)/24;
funcPrm.c5 = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX)/120;
%
funcPrm.expPrm.numTerms = round(exp(2+randn));
for n=1:funcPrm.expPrm.numTerms
	funcPrm.expPrm.matC(:,n) = randn(sizeF,1);
	funcPrm.expPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.expPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecXScale(1,n) = abs(randn);
end
%
funcPrm.vecF0Denom = ones(2,1);
funcPrm.denomExpPrm.numTerms = round(exp(2+randn));
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
vizScl = 1.0;
vizFunch = @(x)(abs(asinh(x*vizScl)/vizScl));
%vizFunch = @(x)( double(abs(x)<0.05) );
cMap = 0.6 +0.4*jet;
cMap(1,:) *= 0.4;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
ax = [-5,5,-5,5];
%ax = [-1.569,-1.568,-2.785,-2.784];
%ax = [3.3727,3.374,-0.7037,-0.7034]
%ax = [1.8,3,3.7,4.5];
nx = 201;
ny = 201;
nc = 20;
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

%vecXA = [4;0];
%vecXA = [-1.6;-2.8];
vecXA = [-1.56846;-2.78451];
vecFA = funchF(vecXA)
vecFAHat = vecFA/sqrt(sum(vecFA.^2,1));
matSolvedPart = eye(2,2)-(vecFAHat*(vecFAHat'));
funchFSolvedPart = @(x)( matSolvedPart * funchF(x) );
funchFSolvedPartNormXY = @(x,y)( sqrt(sum(funchFSolvedPart([x;y]).^2,1)) );
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchFSolvedPartNormXY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc );
colormap(cMap);
axis equal;
grid on;
%
%
%vecXA = [3.3734;-0.7035];
vecXA = [2;4];
vecFA = funchF(vecXA)
vecFAHat = vecFA/sqrt(sum(vecFA.^2,1));
matSolvedPart = eye(2,2)-(vecFAHat*(vecFAHat'));
funchFSolvedPart = @(x)( matSolvedPart * funchF(x) );
funchFSolvedPartNormXY = @(x,y)( sqrt(sum(funchFSolvedPart([x;y]).^2,1)) );
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchFSolvedPartNormXY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc );
colormap(cMap);
axis equal;
grid on;
