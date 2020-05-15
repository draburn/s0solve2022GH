myclear;
%
ax = 5*[-1,1,-1,1];
numXVals = 101;
numYVals = 101;
xVals = linspace(ax(1),ax(2),numXVals);
yVals = linspace(ax(3),ax(4),numYVals);
[ gridX, gridY ] = ndgrid( xVals, yVals );
funchFX = @(x,y)( 1 + x.^2 + 1.5*x.*y );
funchFY = @(x,y)( y.^2 + x.*y );
%
gridFX = funchFX(gridX,gridY);
gridFY = funchFY(gridX,gridY);
gridFN = sqrt(gridFX.^2+gridFY.^2);
min(min(gridFN))
%
figIndx = 0;
cMap = 0.6 + 0.4*jet;
cMap(1,:) *= 0.4;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
numSlices = 20;
%
figIndex++; figure(figIndex);
contourf(gridX,gridY,asinh(gridFN),numSlices);
colormap(cMap)
grid on;
