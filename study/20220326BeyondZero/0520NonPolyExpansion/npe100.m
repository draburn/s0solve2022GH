% "Non-Polynomial Expansion"
% Somebody must have done this before???
% For a matrix M = H + mu*C, pick some mu0, define M0 = H + mu0*C, d = mu - mu0,
%  pick some parameter a...
%   M = ( 1 + a * d ) * M0 - d * ( a*M0 - c ).
clear;
numFigs = 0;
setprngstates(0); sizeX = 5; sizeF = sizeX; sizeB = sizeX; sxexp = 0.0; sfexp = 0.0; sbexp = 0.0; bexp = 0.0; jexp = 0.0; x0exp = 0.0;
%setprngstates(0); sizeX = 100; sizeF = sizeX; sizeB = sizeX; sxexp = 2.0; sfexp = 1.0; sbexp = 1.0; bexp = 1.0; jexp = 1.0; x0exp = 1.0;
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
s = hScale/cScale
%s = 1.0e-4
matR = chol( matH );
bMax = norm( matB * ( matR \ ( matR' \ vecG ) ) );
clear matR;
%
funchVuOfMu = @(mu)( mu + (mu>1.0).*( 2.0 - 1.0./mu - mu ) );
%
muSize = 101;
foo = linspace( 1.0, 0.0, muSize );
muVals = (1.0-foo.^2).^2;
for n=1:muSize
	matR = chol( matH + muVals(n)*s*matC );
	bOfMuVals(n) = norm( matB * ( matR \ ( matR' \ vecG ) ) );
endfor
clear matR;
nuSize = 101;
foo = linspace( 1.0, 0.0, nuSize );
nuVals = 1.0-(1.0-foo.^2).^2;
for n=1:nuSize
	matR = chol( nuVals(n)*matH + s*matC );
	bOfNuVals(n) = nuVals(n) * norm( matB * ( matR \ ( matR' \ vecG ) ) );
endfor
clear matR;
%
muPts = [ muVals(2:end), 1.0./nuVals(1:end-1) ]';
nuPts = [ 1.0./muVals(2:end), nuVals(1:end-1) ]';
vuPts = [ muVals(2:end), 2.0-nuVals(1:end-1) ]';
%
numFigs++; figure(numFigs);
plot( ...
  muVals, bOfMuVals, 'o-', ...
  2.0-nuVals, bOfNuVals, 'o-' );
grid on;
%
%
% Look at impact of alpha.
mu0 = 0.5;
matR0 = chol( matH + mu0*s*matC );
vecY0 = matR0\(matR0'\(-vecG));
b0 = norm(matB*vecY0);
m = s*sumsq(matR0\(matC*vecY0))/(b0^2);
xPts = muPts - mu0;
%
alphaVals = [ 0.0, 0.1*m, 0.4*m, 0.5*m, 0.6*m, 0.9*m, 1.0*m, 1.1*m ];
sizeAlpha = length(alphaVals);
for n=1:sizeAlpha
	alpha = alphaVals(n);
	opaxPts = 1.0 + alpha*xPts;
	bModelPts = b0*( (1.0./opaxPts) + (alpha-m)*xPts./(opaxPts.^2) );
	mskPts = isfinite(bModelPts)&(bModelPts>=0)&(bModelPts<=bMax);
	%
	bModelPtsOfAlphaVals(:,n) = bModelPts;
	mskPtsOfAlphaVals(:,n) = mskPts;
endfor
%
numFigs++; figure(numFigs);
plot( ...
  muVals, bOfMuVals, 'o-', ...
  2.0-nuVals, bOfNuVals, 'o-' );
hold on;
for n=1:sizeAlpha
	plot( vuPts(mskPtsOfAlphaVals(:,n)), bModelPtsOfAlphaVals(mskPtsOfAlphaVals(:,n),n), 'x-' );
endfor
hold off;
grid on;
%
%
% Pick alpha to match mu1...
%
numFigs++; figure(numFigs);
plot( ...
  muVals, bOfMuVals, 'o-', 'linewidth', 3, ...
  2.0-nuVals, bOfNuVals, 'o-', 'linewidth', 3 );
grid on;
mu1 = 2.0;
matR1 = chol( matH + mu1*s*matC );
vecY1 = matR1\(matR1'\(-vecG));
b1 = norm(matB*vecY1)
f1 = b1/b0;
x1 = mu1-mu0;
c2 = f1*(x1^2);
c1 = 2.0*x1*(f1-1.0);
c0 = f1-1.0+(m*x1);
discrim = c1^2 - 4.0*c0*c2
if ( discrim >= 0.0 )
	alphaP = ( -c1 + sqrt(discrim) ) / (2.0*c2)
	alphaM = ( -c1 - sqrt(discrim) ) / (2.0*c2)
	%
	alpha = alphaP;
	opaxPts = 1.0 + alpha*xPts;
	bModelPts = b0*( (1.0./opaxPts) + (alpha-m)*xPts./(opaxPts.^2) );
	mskPts = isfinite(bModelPts)&(bModelPts>=0)&(bModelPts<=bMax);
	bModelPtsOfAlphaP = bModelPts;
	mskPtsOfAlphaP = mskPts;
	%
	alpha = alphaM;
	opaxPts = 1.0 + alpha*xPts;
	bModelPts = b0*( (1.0./opaxPts) + (alpha-m)*xPts./(opaxPts.^2) );
	mskPts = isfinite(bModelPts)&(bModelPts>=0)&(bModelPts<=bMax);
	bModelPtsOfAlphaM = bModelPts;
	mskPtsOfAlphaM = mskPts;
	%
	hold on;
	plot( vuPts(mskPtsOfAlphaP), bModelPtsOfAlphaP(mskPtsOfAlphaP), 'p-', 'linewidth', 5 );
	plot( vuPts(mskPtsOfAlphaM), bModelPtsOfAlphaM(mskPtsOfAlphaM), '-', 'linewidth', 4 );
	plot( funchVuOfMu(mu0), b0, 'o', 'linewidth', 6, 'markersize', 15 );
	plot( funchVuOfMu(mu1), b1, 'o', 'linewidth', 6, 'markersize', 15);
	hold off;
endif
