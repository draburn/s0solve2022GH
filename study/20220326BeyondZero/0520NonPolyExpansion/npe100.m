% "Non-Polynomial Expansion"
% Somebody must have done this before???
% For a matrix M = H + mu*C, pick some mu0, define M0 = H + mu0*C, d = mu - mu0,
%  pick some parameter a...
%   M = ( 1 + a * d ) * M0 - d * ( a*M0 - c ).
clear;
numFigs = 0;
setprngstates(0); sizeX = 5; sizeF = sizeX; sizeB = sizeX; sxexp = 0.0; sfexp = 0.0; sbexp = 0.0; bexp = 0.0; jexp = 0.0; x0exp = 0.0;
%
matSX = diag(exp(sxexp*randn(sizeX,1)));
matSF = diag(exp(sfexp*randn(sizeF,1)));
matSB = diag(exp(sbexp*randn(sizeB,1)));
matJ = matSF*( randn(sizeF,sizeX).*exp(jexp*randn(sizeF,sizeX)) )/matSX;
matB = matSB*( randn(sizeB,sizeX).*exp(bexp*randn(sizeB,sizeX)) ); % Do not divide by matSX, to represent scaling wrong.
vecX0 = matSX*( randn(sizeX,1).*exp(x0exp*randn(sizeX,1)) );
%
vecF = matJ*vecX0;
vecG = matJ'*vecF;
matH = matJ'*matJ;
matC = matB'*matB;
hScale = norm(diag(matH));
cScale = norm(diag(matC));
muScale = hScale/cScale;
matR = chol( matH );
bMax = norm( matB * ( matR \ ( matR' \ vecG ) ) );
%
muSize = 101;
foo = linspace( 1.0, 0.0, muSize );
muVals = (1.0-foo.^2).^2;
for n=1:muSize
	matR = chol( matH + muVals(n)*muScale*matC );
	bOfMuVals(n) = norm( matB * ( matR \ ( matR' \ vecG ) ) );
endfor
nuSize = 101;
foo = linspace( 1.0, 0.0, nuSize );
nuVals = 1.0-(1.0-foo.^2).^2;
for n=1:nuSize
	matR = chol( nuVals(n)*matH + muScale*matC );
	bOfNuVals(n) = nuVals(n) * norm( matB * ( matR \ ( matR' \ vecG ) ) );
endfor
%
muValsAug = [ muVals(2:end), 1.0./nuVals(1:end-1) ];
nuValsAug = [ 1.0./muVals(2:end), nuVals(1:end-1) ];
vuValsAug = [ muVals(2:end), 2.0-nuVals(1:end-1) ];
%
numFigs++; figure(numFigs);
plot( ...
  muVals, bOfMuVals, 'o-', ...
  2.0-nuVals, bOfNuVals, 'o-' );
grid on;
%
%
%
mu0 = 0.5;
matR = chol( matH + mu0*muScale*matC );
vecY0 = matR\(matR'\(-vecG));
b0 = norm(matB*vecY0);
m = sumsq(matR\(matC*vecY0))/(b0^2);
alpha = m;
%
x = muValsAug - mu0;
opax = 1.0 + alpha*x;
bModelOfMuValsAug = b0*( (1.0./opax) + (alpha-m)*x./(opax.^2) );
msk = isfinite(bModelOfMuValsAug)&(bModelOfMuValsAug>=0)&(bModelOfMuValsAug<=bMax);
%
numFigs++; figure(numFigs);
plot( ...
  muVals, bOfMuVals, 'o-', ...
  2.0-nuVals, bOfNuVals, 'o-', ...
  vuValsAug(msk), bModelOfMuValsAug(msk), 'x-' );
grid on;

