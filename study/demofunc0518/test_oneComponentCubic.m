myclear;
%setprngstates(11835232);
%setprngstates();
sizeX = 2;
sizeF = 2;
ax = [-0.3,1.2,-0.5,0.5];
%ax = [-1.5,5,-0.5,0.5];
nx = 51;
ny = 51;
nc = 20;
%
funcPrm.matM = eye(sizeF,sizeX);
funcPrm.vecV = [1;0];
funcPrm.vecA = [2;0];
funcPrm.cQuad = -1.5;
%
vecXRoot = [-0.1;0];
vecFPreRoot = oneComponentCubic(vecXRoot,funcPrm);
%
funchF = @(x)( oneComponentCubic(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
%funchF([1;1])
%
funchF1XY = @(x,y)( funchF([x;y])(1,:) );
funchF2XY = @(x,y)( funchF([x;y])(2,:) );
funchFNXY = @(x,y)( sqrt(sum(funchF([x;y]).^2,1)) );
%
figIndex++; figure(figIndex);
rvecX = linspace(ax(1),ax(2),1001);
plot( rvecX, funchF1XY(rvecX,0*rvecX), 'o-' );
grid on;
%return;
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
%axis equal;
grid on;
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchF2XY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc );
colormap(cMap);
%axis equal;
grid on;
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)( vizFunch(funchFNXY(x,y)) );
contourfunch( funchVizXY, ax, nx, ny, nc )
colormap(cMap);
%axis equal;
grid on;
%
toc
