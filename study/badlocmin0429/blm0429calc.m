myclear;
%setprngstates(0);
thisFile = "blm0429";
%
sizeN = 5;
%
if (0)
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
funchF = @(dummy_matX)( blm0429_func_alt( dummy_matX, matA, matB, vecX0, indexNNonLin, c1, c0 ) );
%
vecX1 = randn(sizeN,1);
vecX2 = randn(sizeN,1);
vecX3 = randn(sizeN,1);
vecX4 = randn(sizeN,1);
vecF1 = funchF(vecX1);
vecF2 = funchF(vecX2);
vecF3 = funchF(vecX3);
vecF4 = funchF(vecX4);
matJ1 = fdjaco( funchF, vecX1 );
matJ2 = fdjaco( funchF, vecX2 );
matJ3 = fdjaco( funchF, vecX3 );
matJ4 = fdjaco( funchF, vecX4 );
toc();
%
matD1 = matJ1-matJ4;
matD2 = matJ2-matJ4;
matD3 = matJ3-matJ4;
vecATemp1 = matD1*randn(sizeN,1);
vecATemp2 = matD2*randn(sizeN,1);
vecATemp3 = matD3*randn(sizeN,1);
% They can differ by a sign.
%vecAHat1 = vecATemp1 / norm(vecATemp1)
%vecAHat2 = vecATemp2 / norm(vecATemp2)
%vecAHat3 = vecATemp3 / norm(vecATemp3)
vecATempS = vecATemp1 + vecATemp2 + vecATemp3;
vecAHat = vecATempS / norm(vecATempS)
%
vecBTemp1 = (matD1')*randn(sizeN,1);
vecBTemp2 = (matD2')*randn(sizeN,1);
vecBTemp3 = (matD3')*randn(sizeN,1);
% They can differ by a sign.
%vecBHat1 = vecBTemp1 / norm(vecBTemp1)
%vecBHat2 = vecBTemp2 / norm(vecBTemp2)
%vecBHat3 = vecBTemp3 / norm(vecBTemp3)
vecBTempS = vecBTemp1 + vecBTemp2 + vecBTemp3;
vecBHat = vecBTempS / norm(vecBTempS)
%
matABHat = vecAHat * (vecBHat');
%
%dphip1 = matD1 ./ matABHat
%dphip1 = sum(sum(matD1 .* matABHat)) / ( eps + sum(sum(matABHat.^2)) )
dphip1 = (vecAHat')*matD1*vecBHat
dphip2 = (vecAHat')*matD2*vecBHat
dphip3 = (vecAHat')*matD3*vecBHat
%
% STUFF BELOW HERE BREAKS DEFINITIONS IN NOTES PRIOR TO 2020.04.30.
y1 = vecBHat'*vecX1
y2 = vecBHat'*vecX2
y3 = vecBHat'*vecX3
y4 = vecBHat'*vecX4
% Take phip = phiC2*y + phiC3*y^2...
% => phi = phiC2*y^2/2 + phiC3*y^3/3.
matDY = [ ...
  y1-y4, (y1^2)-(y4^2); ...
  y2-y4, (y2^2)-(y4^2) ];
vecDPhiP = [ dphip1; dphip2 ];
vecPhiC = matDY \ vecDPhiP
%matDY = [ ...
%  y1-y4, (y1^2)-(y4^2); ...
%  y3-y4, (y3^2)-(y4^2) ];
%vecDPhiP = [ dphip1; dphip3 ];
%vecPhiC = matDY \ vecDPhiP
%
phi1 = vecPhiC(1)*(y1^2)/2.0 + vecPhiC(2)*(y1^3)/3.0
phi1p = vecPhiC(1)*y1 + vecPhiC(2)*(y1^2);
matM = matJ1 - phi1p*matABHat
%
rvecS = linspace(-1000, 1000, 1000001 );
%rvecS = linspace(-404.33,-404.32,10001);
%matXGuess = repmat( ((matM\vecF1)-vecX1), size(rvecS) ) - (matM\vecAHat)*rvecS;
matXGuess = repmat( vecX1-(matM\vecF1), size(rvecS) ) + (matM\vecAHat)*rvecS;
matFGuess = funchF(matXGuess);
rvecOmegaGuess = 0.5*sum(matFGuess.^2,1);
[ omegaMin, indexOfMin ] = min(rvecOmegaGuess)
vecFBestGuess = matFGuess(:,indexOfMin)
vecXBestGuess = matXGuess(:,indexOfMin)
%
%vecDelta = -matM \ (vecF1 - (vecAHat*phi1));
%vecXGuess = vecX1 + vecDelta
echo__vecX0 = vecX0
vecDiff = vecXBestGuess - vecX0
return;

vecPhiC = matPsi \ [ phip1; phip2; phip3 ]
matPsiIntegrated = [ ...
  psi1, (psi1^2)/2.0, (psi1^3)/3.0; ...
  psi2, (psi2^2)/2.0, (psi2^3)/3.0; ...
  psi3, (psi3^2)/2.0, (psi3^3)/3.0 ];
vecPhi = matPsiIntegrated * vecPhiC
%
% Take phi4 = 0. Okay?
matM1 = matJ1 - vecPhi(1)*matABHat
matM2 = matJ2 - vecPhi(2)*matABHat
matM3 = matJ3 - vecPhi(3)*matABHat
matM4 = matJ4

return;
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
