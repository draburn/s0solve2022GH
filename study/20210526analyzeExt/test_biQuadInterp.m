clear;
commondefs;
thisFile = "test_biQuadInterp";
numFigs = 4;
%
%
caseNum = 1;
switch (caseNum)
case 0
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x_secret = -0.890874387338734; % This is the left abs min.
	xVals0 = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
case 1
	bigA = 0.0;
	bigB = 1.0;
	bigX = pi/2.0;
	bigP = 4.0;
	funchF = @(x)( bigA + bigB*abs(x-bigX).^bigP );
	xVals0 = [0.5,1.0,2.2];
	x_secret = bigX;
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
%
% Lets re-do...
xVals1 = [];
for n=1:15
	xVals = sort([ xVals0, xVals1 ]);
	fVals = funchF(xVals);
	xNext = biQuadInterp( xVals, fVals );
	fNext = funchF(xNext);
	msg( thisFile, __LINE__, sprintf( ...
	  "%3d,   %12.8f,   %10.3e", n, xNext, fNext ) );
	xVals1 = [ xVals1, xNext ];
end
%funchFOrig = funchF;
%funchF = @(x)( funchFOrig(x).^(1.0/10.0) );
fVals0 = funchF(xVals0);
fVals1 = funchF(xVals1);

%epsFD = sqrt(sqrt(eps));
%funchDF = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
%funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
viz_numPts = 5000;
%viz_xVals = linspace(-1.6,1.6,viz_numPts);
viz_xVals = linspace(min(xVals1),max(xVals1),viz_numPts);

numFigs++; figure(numFigs);
plot( ...
  viz_xVals, funchF(viz_xVals), 'o-', ...
  xVals1, fVals1, 'ko-', 'linewidth', 2, 'markersize', 25, ...
  xNext, fNext, 'g*', 'linewidth', 4, 'markersize', 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );

numFigs++; figure(numFigs);
semilogy( ...
  fVals1, 'o-', ...
  abs(xVals1-x_secret), 'x-' );
grid on;
xlabel( "n" );
ylabel( "f" );
title( "f vs n" );
