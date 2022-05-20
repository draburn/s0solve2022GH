% "Non-Polynomial Expansion"
% Somebody must have done this before???
% For a matrix M = H + mu*C, pick some mu0, define M0 = H + mu0*C, d = mu - mu0,
%  pick some parameter a...
%   M = ( 1 + a * d ) * M0 - d * ( a*M0 - c ).
clear;
numFigs = 0;
c = 1.0;
x0 = 1.0;
f00 = 1.0/(1.0+(c*x0));

%

numPts = 101;
alphaVals = (c/f00)*[ 1.0, 0.0, 0.1, 0.5, 0.9, 1.1 ];
vecDeltaX = linspace(-0.2,2.0,numPts)';
numAlphaVals = length(alphaVals);
for n=1:numAlphaVals;
	alpha = alphaVals(n);
	vecDenom = (1.0 + alpha*vecDeltaX)*f00;
	vecNumer = vecDeltaX*(alpha-(c/f00));
	r = vecNumer./vecDenom;
	f0(:,n) = 1.0./vecDenom;
	f1(:,n) = f0(:,n) + r.*f0(:,n);
	f2(:,n) = f1(:,n) + r.*r.*f0(:,n);
endfor
%
numFigs++; figure(numFigs);
plot( 
  x0+vecDeltaX, f0(:,1), 'k-', "linewidth", 5, ...
  x0+vecDeltaX, f0(:,2), 'o-', "linewidth", 2, ...
  x0+vecDeltaX, f0(:,3), '+-', "linewidth", 2, ...
  x0+vecDeltaX, f0(:,4), 's-', "linewidth", 2, ...
  x0+vecDeltaX, f0(:,5), 'v-', "linewidth", 2, ...
  x0+vecDeltaX, f0(:,6), '^-', "linewidth", 2 );
grid on;
%
numFigs++; figure(numFigs);
plot( 
  x0+vecDeltaX, f1(:,1), 'k-', "linewidth", 5, ...
  x0+vecDeltaX, f1(:,2), 'o-', "linewidth", 2, ...
  x0+vecDeltaX, f1(:,3), '+-', "linewidth", 2, ...
  x0+vecDeltaX, f1(:,4), 's-', "linewidth", 2, ...
  x0+vecDeltaX, f1(:,5), 'v-', "linewidth", 2, ...
  x0+vecDeltaX, f1(:,6), '^-', "linewidth", 2 );
grid on;
%
numFigs++; figure(numFigs);
plot( 
  x0+vecDeltaX, f2(:,1), 'k-', "linewidth", 5, ...
  x0+vecDeltaX, f2(:,2), 'o-', "linewidth", 2, ...
  x0+vecDeltaX, f2(:,3), '+-', "linewidth", 2, ...
  x0+vecDeltaX, f2(:,4), 's-', "linewidth", 2, ...
  x0+vecDeltaX, f2(:,5), 'v-', "linewidth", 2, ...
  x0+vecDeltaX, f2(:,6), '^-', "linewidth", 2 );
grid on;
