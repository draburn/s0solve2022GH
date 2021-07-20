clear;
commondefs;
%
numPts=10;
x=sort(randn(1,numPts));
f=abs(randn(1,numPts))+0.1*x;
%f = (x+0.1*randn(1,numPts)).^3;
xNext=findGoodCand(x,f);
plot( ...
  [min(x), max(x)], [0.0,0.0], 'k-', ...
  x, f, 'o-', 'markersize', 10, ...
  xNext*[1,1], [0,0], 'rx-', 'markersize',15);
grid on;
return;
%
sizeX = 1;
sizeF = sizeX;
seedPrm = demoFunc0101_genSeedPrm("xhard");
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
funchF = @(vecXDummy)( demoFunc0101_eval( 1E4*vecXDummy, funcPrm ) );
%
xVals = linspace(-5,5,101);
fVals = funchF(xVals);
%
plot( xVals, fVals, 'o-' );
grid on;
