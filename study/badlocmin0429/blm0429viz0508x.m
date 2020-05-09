myclear;
setprngstates(0);
thisFile = "blm0429";
%
sizeN = 2;
%
if (0)
	matA = eye(sizeN,sizeN);
	matB = eye(sizeN,sizeN);
	vecX0 = zeros(sizeN,1);
	indexNNonLin = 1;
elseif (0)
	matA = eye(sizeN,sizeN) + 0.1*randn(sizeN,sizeN);
	matB = eye(sizeN,sizeN) + 0.1*randn(sizeN,sizeN);
	vecX0 = randn(sizeN,1);
	indexNNonLin = 1;
elseif (1)
	matA = eye(sizeN,sizeN) + 0.5*randn(sizeN,sizeN);
	matB = eye(sizeN,sizeN);% + 0.5*randn(sizeN,sizeN);
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
funchF = @(x)(blm0429_func( x, matA, matB, vecX0, indexNNonLin, c1, c0 ));
%
if (0)
%vecXM = [ -0.47277667109; 0.0 ]
%vecXM = matB*[ -0.472776671; 0.0 ]
vecXM = [-0.473;0.003];
%vecXM = [-0.47276;0.003];
vecFM = funchF(vecXM);
fdjaco_prm.epsFD = 1E-6;
fdjaco_prm.fdOrder = 2;
matJM = fdjaco( funchF, vecXM, fdjaco_prm );
if (0)
	vecXTest = vecXM + 1e-8*[1;1];
	vecFTest1 = funchF(vecXTest);
	vecFTest2 = vecFM + (matJM*[1;1]*1e-8);
	norm(vecFTest2-vecFTest1)
end
[matUM,matSM,matVM] = svd(matJM);
diag(matSM)
vecVNull = matVM(:,end);
numLambdaVals = 10001;
%rvecLambda = linspace(-10,10,numLambdaVals);
rvecLambda = linspace(-1,1,numLambdaVals);
%rvecLambda = linspace(-0.47278,-0.47276,numLambdaVals);
matX = vecVNull*rvecLambda + repmat(vecXM,[1,numLambdaVals]);
matF = funchF(matX);
rvecFNorm = sqrt(sum(matF.^2,1));
semilogy(rvecLambda,rvecFNorm,'o-');
grid on;
return;
end
%
numX1Vals = 201;
numX2Vals = 201;
%
x1Vals = linspace(-0.8,0.2,numX1Vals);
x2Vals = linspace(-0.02,0.02,numX2Vals);
%
%x1Vals = linspace(-0.6,-0.4,numX1Vals);
%x2Vals = linspace(-0.0,0.01,numX2Vals);
%
%x1Vals = linspace(-0.474,-0.472,numX1Vals);
%x2Vals = linspace(0.0028,0.0032,numX2Vals);
%
%x1Vals = linspace(-0.4729,-0.4727,numX1Vals);
%x2Vals = linspace(-0.00001,0.00001,numX2Vals);
%
[ gridX1, gridX2 ] = ndgrid( x1Vals, x2Vals );
matX = [ ...
  reshape( gridX1, [1,numX1Vals*numX2Vals] );
  reshape( gridX2, [1,numX1Vals*numX2Vals] ) ];
%vecX = randn(2,1)
%blm0429_func_alt( vecX, matA, matB, vecX0, indexNNonLin, c1, c0 )
%blm0429_func( vecX, matA, matB, vecX0, indexNNonLin, c1, c0 )
%return
%matF = blm0429_func( matX, matA, matB, vecX0, indexNNonLin, c1, c0 );
matF = funchF(matX);
gridF1 = reshape( matF(1,:), [numX1Vals,numX2Vals] );
gridF2 = reshape( matF(2,:), [numX1Vals,numX2Vals] );
gridOmega = 0.5*(gridF1.^2+gridF2.^2);
%
%funchViz = @(f)( sign(f).*( abs(f).^(1.0/2.0) ) );
funchViz = @(f)( ( abs(f).^(1.0/2.0) ) );
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
