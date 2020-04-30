myclear;
thisFile = "blm0429";
%
sizeN = 2;
%
if (1)
	matA = eye(sizeN,sizeN);
	matB = eye(sizeN,sizeN);
	vecX0 = zeros(sizeN,1);
	indexNNonLin = 1;
else
	matA = randn(sizeN,sizeN);
	matB = randn(sizeN,sizeN);
	vecX0 = randn(sizeN,1);
	indexNNonLin = ceil((10*eps) + ((sizeN-(20*eps))*rand) );
end
%
c1 = 1.0;
cx = 0.1;
c0 = (c1^2)*(1.0+cx)/4.0;
%
numX1Vals = 201;
numX2Vals = 101;
x1Vals = linspace(-1,0.5,numX1Vals);
x2Vals = linspace(-0.1,0.1,numX2Vals);
[ gridX1, gridX2 ] = ndgrid( x1Vals, x2Vals );
matX = [ ...
  reshape( gridX1, [1,numX1Vals*numX2Vals] );
  reshape( gridX2, [1,numX1Vals*numX2Vals] ) ];
%vecX = randn(2,1)
%blm0429_func_alt( vecX, matA, matB, vecX0, indexNNonLin, c1, c0 )
%blm0429_func( vecX, matA, matB, vecX0, indexNNonLin, c1, c0 )
%return
matF = blm0429_func( matX, matA, matB, vecX0, indexNNonLin, c1, c0 );
gridF1 = reshape( matF(1,:), [numX1Vals,numX2Vals] );
gridF2 = reshape( matF(2,:), [numX1Vals,numX2Vals] );
gridOmega = 0.5*(gridF1.^2+gridF2.^2);
%
funchViz = @(f)( sign(f).*( abs(f).^(1.0/2.0) ) );
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
contourf( gridX1, gridX2, funchViz(sqrt(sqrt(gridOmega))), 20 );
colormap(cmap);
grid on;
