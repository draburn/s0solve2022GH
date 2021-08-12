clear;
commondefs;
thisFile = "interpo";
numFigs = 0;
setprngstates(0);
%
funchF = @(x)( x.^2 );
numPts = 10;
xVals = randn(numPts,1);
%
sort(xVals);
fVals = funchF(xVals);
%
r1 = 0.1;
r0 = 0.01;
fVals .*= 1.0 + r1*randn(size(fVals));
fVals += r0*randn(size(fVals));
%
%
%
mod_numPts = 50;
modXVals = linspace( min(xVals), max(xVals), mod_numPts );
modX0 = 0.0001;
matW = zeros(mod_numPts,numPts);
for n=1:numPts
	matW(:,n) = 1.0./( modX0^2 + (modXVals-xVals(n)).^2 );
end
for n=1:mod_numPts
	matW(n,:) /= (eps+sum(abs(matW(n,:))));
end
modFVals = matW*fVals;
%
%
%
xLo = min([ min(xVals) ]);
xHi = max([ max(xVals) ]);
viz_numPts = 1000;
viz_xVals = linspace( xLo, xHi, viz_numPts );
%
numFigs++; figure(numFigs);
plot( ...
  modXVals, modFVals, 'o-', 'linewidth', 2, ...
  viz_xVals, funchF(viz_xVals), 'k-', ...
  xVals, fVals, 'ko', 'linewidth', 2, 'markersize', 25, ...
  xVals, fVals, 'k+', 'linewidth', 2, 'markersize', 25 );
grid on;