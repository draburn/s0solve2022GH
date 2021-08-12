clear;
commondefs;
thisFile = "demo";
numFigs = 0;
setprngstates(0);
%
bigF0 = randn()
bigF1 = randn()
p = abs(randn())
s = randn()
funchF = @(x)( bigF0 + bigF1 * abs( x - s ).^p );
%
numPts = 5;
xVals = randn(1,numPts)
%
xVals = sort(xVals);
fVals = funchF(xVals);
%
r1 = 0;%0.01;
r0 = 0;%0.001;
fVals .*= 1.0 + r1*randn(size(fVals));
fVals += r0*randn(size(fVals));
%
xLo = min([ s, min(xVals) ]);
xHi = max([ s, max(xVals) ]);
viz_numPts = 1000;
viz_xVals = linspace( xLo, xHi, viz_numPts );
%
numFigs++; figure(numFigs);
plot( ...
  xVals, fVals, 'ko', 'linewidth', 2, 'markersize', 20, ...
  viz_xVals, funchF(viz_xVals), 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
%
df = diff(fVals)./diff(xVals);
cx = cent(xVals);
ddf = diff(df)./diff(cx);
ccx = cent(cx);
h = cent(df)./abs(ddf);
%
numFigs++; figure(numFigs);
plot( ...
  ccx, h, 'o-' );
grid on;
xlabel( "x" );
ylabel( "h" );
title( "h vs x" );
