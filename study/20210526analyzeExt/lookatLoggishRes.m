clear;
commondefs;
thisFile = "lookatLoggishRes";
%setprngstates(3);
setprngstates(33754416);
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
case 4
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 + vecP(4)*x.^3 + vecP(5)*x.^4 + vecP(6)*x.^5 );
otherwise
	error(["Invalid value of modelOrder (", num2str(modelOrder), ")."]);
end

if (1)
	msg( thisFile, __LINE__, "OVER-RIDING FUNCTION." );
	%funchF = @(x)( exp(-1.0./(x.^2)) );
	%funchF = @(x)( 1E-4 + abs(x-0.0).^20.5 );
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
end

funchDF  = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
%
numPts = modelOrder+1;
xVals = sort(randn(1,numPts));
if (1)
	msg( thisFile, __LINE__, "OVER-RIDING X-VALUES." );
	%xVals = [-2.1288  -1.1109   1.4579];
	% Zero is at 0.173.
	% [-2.1288  -1.1109   0.173    1.4579] = [ 5.3305e+06   8.6366e+00   1.0000e-04   2.2721e+03 ];
	% ptwise min is 3. drop 1.
	%xVals = [-1.1109    0.173   1.4579];
	% Zero at -0.464.
	%[-1.1109  -0.464   0.173   1.4579] = [8.6366e+00   1.0015e-04   1.0000e-04   2.2721e+03];
	% ptwise min is 3. drop 1.
	%xVals = [-0.464   0.173   1.4579];
	% zero at -0.145
	%[-0.464   -0.145  0.173   1.4579] = [ 1.0015e-04   1.0000e-04   1.0000e-04   2.2721e+03];
	% = 1.00145753368791e-04   1.00000000000006e-04   1.00000000000240e-04   2.27205377617069e+03
	% ptwise min is 2. drop 4.
	%xVals = [-0.464   -0.145  0.173 ];
	% zero at 0.014;
	% [-0.464   -0.145  0.014  0.173 ] =
	% 1.00145753368791e-04   1.00000000000006e-04   1.00000000000000e-04   1.00000000000240e-04
	% 2 is ptwise min.
	% And, it's pretty clear the ext val is 1E-4.
	% And, this happens simply with quad model?!?!
	%
	%
	%xVals = [-1.5,-1.0,-0.5]; %Next is -0.85 as pt 3. min is 3. Drop 1.
	%xVals = [-1.0,-0.85,-0.5]; %Next is -0.917 as pt 2. min is 3. Drop 1.
	%xVals = [-0.917,-0.85,-0.5]; %Next is -0.8835 as pt 2. min is 2. Drop 4.
	%xVals = [-0.917, -0.8835, -0.85]; %Next is -0.8823 as pt 3. min is 3. Drop 1.
	%xVals = [-0.8835, -0.8823, -0.85]; %Next is -0.8828...
	%%%
	% Originally looked like bad pts,
	% with pt on right being continually used.
	% Maybe this only happens some of the time?
	%xVals = [ -0.9, -0.85, -0.5]; %-0.875 as 2. Below tol.
	% Try again...
	xVals = [ -1.4, -0.82, -0.8 ]; %-0.81! On wrong side!
	xVals = [ -1.4, -0.82, -0.815 ]; %-0.821! Barely on other side!
	% Look at this case with correct points!
end

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
case 4
	funchFModel = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 + vecC(4)*x.^3 + vecC(5)*x.^4 );
otherwise
	error(["Invalid value of modelOrder (", num2str(modelOrder), ")."]);
end
funchDFModel  = @(x)( (funchFModel(x+epsFD)-funchFModel(x-epsFD))/(2.0*epsFD) );
funchDDFModel = @(x)( (funchFModel(x+epsFD)+funchFModel(x-epsFD)-2*funchFModel(x))/(epsFD^2) );
%
viz_numPts = 1000;
viz_xVals = linspace(-1.6,1.6,viz_numPts);
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
  viz_xVals, viz_dfModelVals./(0*max(abs(viz_dfModelVals)) + abs(viz_ddfModelVals)), 'v-', ...
  viz_xVals, 0*viz_xVals, 'k-' );
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
