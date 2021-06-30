clear;
commondefs;
thisFile = "test_biQuadInterp";
numFigs = 0;
%
%
caseNum = 0;
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
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
%
% Lets re-do...
xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
for n=1:15
	fVals = funchF(xVals);
	xNext = biQuadInterp( xVals, fVals );
	fNext = funchF(xNext);
	msg( thisFile, __LINE__, sprintf( ...
	  "%3d,   %12.8f,   %10.3e", n, xNext, fNext ) );
	xVals = sort([ xNext, xVals ]);
end
fVals = funchF(xVals);

epsFD = sqrt(sqrt(eps));
funchDF = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );

viz_numPts = 5000;
viz_xVals = linspace(-1.6,1.6,viz_numPts);

numFigs++; figure(numFigs);
plot( ...
  viz_xVals, funchF(viz_xVals), 'o-', ...
  xVals, fVals, 'ko', 'linewidth', 2, 'markersize', 25, ...
  xNext, fNext, 'g*', 'linewidth', 4, 'markersize', 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
