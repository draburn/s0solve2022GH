clear;
%
setprngstates
%setprngstates(0);
%setprngstates(67938960);
numYTerms = 4;
yScale = 50.0;
ya0_orig = 0.0;
ya_orig = abs(randn(1,numYTerms));
ya0 = ya0_orig/(ya0_orig+sum(ya_orig));
ya = ya_orig/(ya0_orig+sum(ya_orig));
yf = 2*pi*randn(1,numYTerms)/yScale;
yd = 2*pi*rand(1,numYTerms);
funch_y = @(x)( ya0*x + ...
 ya(1)*( x + cos(yf(1)*x+yd(1))/yf(1) ) + ...
 ya(2)*( x + cos(yf(2)*x+yd(2))/yf(2) ) + ...
 ya(3)*( x + cos(yf(3)*x+yd(3))/yf(3) ) + ...
 ya(4)*( x + cos(yf(4)*x+yd(4))/yf(4) ) );
%
numFTerms = 5;
fScale = 1.0;
fa = randn(1,numFTerms);
ff = 2*pi*randn(1,numFTerms)/fScale;
fd = 2*pi*rand(1,numFTerms);
funch_f = @(x)( ...
  fa(1)*cos(ff(1)*x+fd(1)) + ...
  fa(2)*cos(ff(2)*x+fd(2)) + ...
  fa(3)*cos(ff(3)*x+fd(3)) + ...
  fa(4)*cos(ff(4)*x+fd(4)) + ...
  fa(5)*cos(ff(5)*x+fd(5)) );
%
numFigs = 0;
%
xVals = linspace(-500,500,10001);
numFigs++; figure(numFigs);
plot( cent(xVals), diff(funch_y(xVals))./diff(xVals), 'o-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, xVals, 'kx-', xVals, funch_y(xVals), 'mx-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
%
yVals = funch_y(xVals);
dydxVals = diff(yVals)./diff(xVals);
[ foo, iOfMin ] = min(dydxVals);
xOfMin = (xVals(iOfMin)+xVals(iOfMin+1))/2.0;
[ foo, iOfMax ] = max(dydxVals);
xOfMax = (xVals(iOfMax)+xVals(iOfMax+1))/2.0;
%
xLo = xOfMin-0.05*(xOfMax-xOfMin);
xHi = xOfMin+0.05*(xOfMax-xOfMin);
xVals = linspace(xLo,xHi,10001);
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
%
xLo = xOfMax-0.05*(xOfMax-xOfMin);
xHi = xOfMax+0.05*(xOfMax-xOfMin);
xVals = linspace(xLo,xHi,10001);
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
%
xLo = xOfMax+xOfMin-0.05*(xOfMax-xOfMin);
xHi = xOfMax+xOfMax+0.05*(xOfMax-xOfMin);
xVals = linspace(xLo,xHi,10001);
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
