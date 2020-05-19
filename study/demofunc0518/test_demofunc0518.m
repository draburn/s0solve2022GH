myclear;
sizeX = 2;
sizeF = 2;
funcPrm.vecF0 = zeros(sizeF,1);
funcPrm.matJ0 = eye(sizeF,sizeX);
funcPrm.polyOrder = 5;
funcPrm.c2 = randn(sizeF,sizeX,sizeX);
funcPrm.c3 = randn(sizeF,sizeX,sizeX,sizeX);
funcPrm.c4 = randn(sizeF,sizeX,sizeX,sizeX,sizeX);
funcPrm.c5 = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX);
%
funcPrm.expPrm.numTerms = 0;
for n=1:0
	funcPrm.expPrm.numTerms++;
	funcPrm.expPrm.matC(:,n) = randn(sizeF,1);
	funcPrm.expPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.expPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecXScale(1,n) = abs(randn);
end
%
funcPrm.vecF0Denom = ones(2,1);
funcPrm.denomExpPrm.numTerms = 0;
for n=1:0
	funcPrm.denomExpPrm.numTerms++;
	funcPrm.denomExpPrm.matC(:,n) = abs(randn(sizeF,1));
	funcPrm.denomExpPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.denomExpPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.denomExpPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.denomExpPrm.rvecXScale(1,n) = abs(randn);
end
%
vecXRoot = [0;0];
vecFPreRoot = demofunc0518_eval(vecXRoot,funcPrm);
%
funchF = @(x)( demofunc0518_eval(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
funchF([1;1])
%
funchF1XY = @(x,y)( funchF([x;y])(1,:) );
funchF2XY = @(x,y)( funchF([x;y])(2,:) );
funchFNXY = @(x,y)( sqrt(sum(funchF([x;y]).^2,1)) );
%
vizScl = 1.0;
vizFunch = @(x)(asinh(x));
funchVizXY = @(x,y)( vizFunch(vizScl*funchFNXY(x,y))/vizScl );
%contourfunch( funchVizXY, [1.725,1.727,-0.785,-0.783] )
contourfunch( funchVizXY, 3*[-1,1,-1,1] );
cMap = 0.6 +0.4*jet;
cMap(1,:) *= 0.4;
cMap(end,:) *= 0.4;
cMap(1,:) += 0.6;
colormap(cMap);
axis equal;
grid on;
toc
