% Try a linear fit.
clear;
commondefs;
thisFile = "interpo";
numFigs = 0;
setprngstates(0);
%
funchF = @(x)( x.^4 );
numPts = 11;
%xVals = randn(numPts,1);
%xVals = linspace(-5,5,numPts)';
xVals = [ -5.0, -3.0, -2.0, -2.0, -2.0, 3.0, 3.0, 2.0, 1.0, 1.0, 4.0 ]';
%
xVals = sort(xVals);
fVals = funchF(xVals);
%
wVals0 = ones(size(xVals));
for n=1:numPts
	wVals0(n) = sum((abs(xVals-xVals(n))).^2);
end
%
r1 = 0.1;
r0 = 0.01;
fVals .*= 1.0 + r1*randn(size(fVals));
fVals += r0*randn(size(fVals));
%
%
%
modNumPts = 500;
modXLo = min(xVals);
modXHi = max(xVals);
modXVals = linspace( modXLo, modXHi, modNumPts );
modX0 = 0.1;
matW = zeros(modNumPts,numPts);
for n=1:numPts
	%matW(:,n) = wVals0(n)./( modX0^2 + (modXVals-xVals(n)).^2 );
	matW(:,n) = wVals0(n)*exp( -0.5*abs( (modXVals-xVals(n))/modX0 ) );
end
for n=1:modNumPts
	matW(n,:) /= (eps+sum(abs(matW(n,:))));
end
%modFVals = matW*fVals;
%%%matX = [ ones(size(xVals)), xVals, xVals.^2 ];
matX = [ ones(size(xVals)), xVals ];
vecF = fVals;
for n=1:modNumPts
	matD = diag(sqrt(matW(n,:)));
	vecC = (matD*matX)\(matD*vecF);
	%%%modFVals(n) = vecC(1) + vecC(2)*modXVals(n) + vecC(3)*modXVals(n)^2;
	%%%modHVals(n) = ( vecC(2) + 2.0*vecC(3)*modXVals(n) ) / abs(2.0*vecC(3));
	% Above modH does not include deriv due to weights.
	%
	modFVals(n) = vecC(1) + vecC(2)*modXVals(n);
	modHVals(n) = 0.0;
end
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
%
dfVals = diff(modFVals)./diff(modXVals);
cxVals = cent(modXVals);
ddfVals = diff(dfVals)./diff(cxVals);
ccxVals = cent(cxVals);
hVals = cent(dfVals)./(eps*max(abs(ddfVals)) + abs(ddfVals));
sortAbsModHVals = sort(abs(hVals));
hCap = sortAbsModHVals(1+round(2.0*(modNumPts-1)/3.0));
numFigs++; figure(numFigs);
plot( ...
  modXVals, cap(modHVals,-hCap,hCap), 'x-', 'linewidth', 2, ...
  ccxVals, cap(hVals,-hCap,hCap), 'o-', 'linewidth', 2, ...
  xVals, 0*xVals, 'ko-' );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  xVals, wVals0, 'o-', 'linewidth', 2 );
grid on;