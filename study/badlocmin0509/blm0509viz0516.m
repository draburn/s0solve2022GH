blm0509init0516;
%
numX1Vals = 201;
numX2Vals = 201;
numLevels = 20;
x1Vals = linspace(-2,2,numX1Vals);
x2Vals = linspace(-2,2,numX2Vals);
%x1Vals = linspace(-0.65,0.05,numX1Vals);
%x2Vals = linspace(-0.1,0.05,numX2Vals);
%x1Vals = linspace(-0.55,-0.45,numX1Vals);
%x2Vals = linspace(-0.055,-0.035,numX2Vals);
%x1Vals = linspace(-0.505,-0.495,numX1Vals);
%x2Vals = linspace(-0.047,-0.044,numX2Vals);
%%%matC = [ -0.00999950121796285; 0.99995000373788323 ]; %Columnspace Jm.
%
vecXM = [-0.5;-0.045];
vecFM = funchF(vecXM);
vecFMHat = vecFM/sqrt(sum(vecFM.^2));
%
matJM = fdjaco( funchF, vecXM );
[ matUM, matSM, matVM ] = svd( matJM );
matC = matUM(:,1)
%
%x1Vals = linspace(-0.7,0.02,numX1Vals);
%x2Vals = linspace(-0.1,0.02,numX2Vals);
%x1Vals = linspace(-0.6,-0.4,numX1Vals);
%x2Vals = linspace(-0.06,-0.02,numX2Vals);
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
%Psuedo-root curve...
%gridPRC = reshape( (matUM(:,1)')*matF, [numX1Vals,numX2Vals] );
gridPRC = reshape( sqrt(sum((matF-vecFMHat*(vecFMHat'*matF)).^2,1)), [numX1Vals,numX2Vals] );
%
%funchViz = @(f)( sign(f).*( abs(f).^(1.0/2.0) ) );
%funchViz = @(f)( ( abs(f).^(1.0/2.0) ) );
funchViz = @(f)( abs(asinh(f*100.0)/100.0) );
%funchViz = @(f)( abs(asinh(f*1000.0)/1000.0) );
%funchViz = @(f)( abs(asinh(f*10000.0)/10000.0) );
cMap = 0.6 + 0.4*jet(1000);
cMap(1,:) *= 0.6;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(gridF1), numLevels );
%contourf( gridX1, gridX2, funchViz(gridPRC) );
colormap(cMap);
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(gridF2), numLevels );
colormap(cMap);
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(sqrt(gridOmega)), numLevels );
colormap(cMap);
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX1, gridX2, funchViz(gridPRC), numLevels );
colormap(cMap);
grid on;
