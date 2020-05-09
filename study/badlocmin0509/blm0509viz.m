blm0509init;
thisFile = "blm0509viz";
%
numX1Vals = 201;
numX2Vals = 201;
%x1Vals = linspace(-2,2,numX1Vals);
%x2Vals = linspace(-2,2,numX2Vals);
%x1Vals = linspace(-0.9,0.1,numX1Vals);
%x2Vals = linspace(-0.1,0.7,numX2Vals);
x1Vals = linspace(-0.55,0.05,numX1Vals);
x2Vals = linspace(-0.05,0.35,numX2Vals);
%
[ gridX1, gridX2 ] = ndgrid( x1Vals, x2Vals );
matX = [ ...
  reshape( gridX1, [1,numX1Vals*numX2Vals] );
  reshape( gridX2, [1,numX1Vals*numX2Vals] ) ];
matF = funchF(matX);
gridF1 = reshape( matF(1,:), [numX1Vals,numX2Vals] );
gridF2 = reshape( matF(2,:), [numX1Vals,numX2Vals] );
gridOmega = 0.5*(gridF1.^2+gridF2.^2);
%
%funchViz = @(f)( sign(f).*( abs(f).^(1.0/2.0) ) );
%funchViz = @(f)( ( abs(f).^(1.0/2.0) ) );
funchViz = @(f)( abs(asinh(f*1000.0)/1000.0) );
cmap = colormap(jet(1000));
cmap(1,:) = 0.7;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(gridF1) );
colormap(jet(1000));
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(gridF2) );
colormap(jet(1000));
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(sqrt(gridOmega)), 20 );
colormap(cmap);
grid on;
