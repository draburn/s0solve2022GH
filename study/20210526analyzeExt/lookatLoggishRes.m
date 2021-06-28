clear;
commondefs;
thisFile = "lookatLoggishRes";
%setprngstates(2);
setprngstates();
numFigs = 0;
%
modelOrder = 2;
epsFD = sqrt(sqrt(eps));
vecP = randn(modelOrder+2,1);
vecP(end) *= 0.1;
switch (modelOrder)
case 0
	funchF = @(x)( vecP(1) + vecP(2)*x );
case 1
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 );
case 2
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 + vecP(4)*x.^3 );
case 3
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 + vecP(4)*x.^3 + vecP(5)*x.^4 );
otherwise
	error(["Invalid value of modelOrder (", num2str(modelOrder), ")."]);
end
funchDF  = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
%
numPts = modelOrder+1;
xVals = sort(randn(1,numPts));
fVals = funchF(xVals);
%
vecX = xVals';
matX = ones(numPts,1);
for m=1:modelOrder
	matX = [ matX, vecX.^m ];
end
vecF = fVals';
vecC = matX\vecF;
switch (modelOrder)
case 0
	funchFModel = @(x)( vecC(1) );
case 1
	funchFModel = @(x)( vecC(1) + vecC(2)*x );
case 2
	funchFModel = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
case 3
	funchFModel = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 + vecC(4)*x.^3 );
otherwise
	error(["Invalid value of modelOrder (", num2str(modelOrder), ")."]);
end
funchDFModel  = @(x)( (funchFModel(x+epsFD)-funchFModel(x-epsFD))/(2.0*epsFD) );
funchDDFModel = @(x)( (funchFModel(x+epsFD)+funchFModel(x-epsFD)-2*funchFModel(x))/(epsFD^2) );
%
viz_numPts = 1000;
viz_xVals = linspace(min(xVals),max(xVals),viz_numPts);
viz_fVals   = funchF(   viz_xVals );
viz_dfVals  = funchDF(  viz_xVals );
viz_ddfVals = funchDDF( viz_xVals );
viz_fModelVals   = funchFModel(   viz_xVals );
viz_dfModelVals  = funchDFModel(  viz_xVals );
viz_ddfModelVals = funchDDFModel( viz_xVals );
%
xAvg = sum(xVals)/numPts;
xMatchDDFforMO2 = xAvg;
xSqAvg = sum(xVals.^2)/numPts;
xVarSq = xSqAvg - (xAvg^2);
assert( xVarSq > 0.0 );
xVar = sqrt(xVarSq);
xMatchDFforMO2m = xAvg - xVar/sqrt(2.0);
xMatchDFforMO2p = xAvg + xVar/sqrt(2.0);
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_fVals, '^-', ...
  viz_xVals, viz_fModelVals, 'v-', ...
  xVals, fVals, 'ko', 'linewidth', 4, 'markersize', 25 );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_dfVals, '^-', ...
  viz_xVals, viz_dfModelVals, 'v-', ...
  xMatchDFforMO2m, funchDFModel(xMatchDFforMO2m), 'cs', 'linewidth', 4, 'markersize', 25,
  xMatchDFforMO2p, funchDFModel(xMatchDFforMO2p), 'cs', 'linewidth', 4, 'markersize', 25 );
grid on;
xlabel( "x" );
ylabel( "f'" );
title( "f' vs x" );
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_ddfVals, '^-', ...
  viz_xVals, viz_ddfModelVals, 'v-', ...
  xMatchDDFforMO2, funchDDFModel(xMatchDDFforMO2), 'cs', 'linewidth', 4, 'markersize', 25 );
grid on;
xlabel( "x" );
ylabel( "f''" );
title( "f'' vs x" );
%
numFigs++; figure(numFigs);
%plot( ...
%  viz_xVals, viz_dfVals./(0.001*max(abs(viz_dfVals)) + abs(viz_ddfVals)), '^-', ...
%  viz_xVals, viz_dfModelVals./(0.001*max(abs(viz_dfModelVals)) + abs(viz_ddfModelVals)), 'v-' );
plot( ...
  viz_xVals, viz_dfVals./(0*max(abs(viz_dfVals)) + abs(viz_ddfVals)), '^-', ...
  viz_xVals, viz_dfModelVals./(0*max(abs(viz_dfModelVals)) + abs(viz_ddfModelVals)), 'v-' );
grid on;
xlabel( "x" );
ylabel( "f'/|f''|" );
title( "f'/|f''| vs x" );
ax0 = axis();
ax0(3) = max([ax0(3),-2.0]);
ax0(4) = min([ax0(4),+2.0]);
axis(ax0);
%
return;
