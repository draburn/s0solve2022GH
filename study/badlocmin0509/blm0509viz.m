blm0509init;
thisFile = "blm0509viz";
%
numX1Vals = 201;
numX2Vals = 201;
x1Vals = linspace(-2,2,numX1Vals);
x2Vals = linspace(-2,2,numX2Vals);
x1Vals = linspace(-0.6,0.05,numX1Vals);
x2Vals = linspace(-0.02,0.05,numX2Vals);
%x1Vals = linspace(-0.5,-0.45,numX1Vals);
%x2Vals = linspace(0.006,0.014,numX2Vals);
%
x1Vals = linspace(-0.7,0.02,numX1Vals);
x2Vals = linspace(-0.1,0.02,numX2Vals);
x1Vals = linspace(-0.6,-0.4,numX1Vals);
x2Vals = linspace(-0.06,-0.02,numX2Vals);
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
funchViz = @(f)( abs(asinh(f*100.0)/100.0) );
%funchViz = @(f)( abs(asinh(f*1000.0)/1000.0) );
%funchViz = @(f)( abs(asinh(f*10000.0)/10000.0) );
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
contourf( gridX1, gridX2, funchViz(sqrt(gridOmega)) );
colormap(cmap);
grid on;
