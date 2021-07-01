clear;
commondefs;
thisFile = "lookatLoggishRes2";
numFigs = 0;
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
%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843 ]);
xVals = sort([ -0.831698038563058, xVals ]);
xVals = sort([ -0.897620231966365, xVals ]);
xVals = sort([ -0.867992251154230, xVals ]);
xVals = sort([ -0.88284017185230, xVals ]);
xVals = sort([ -0.887137697860992, xVals ]); % xB goes past data.
xVals = sort([ -0.886078607798801, xVals ]);
xVals = sort([ -0.890138678468405, xVals ]); % xB goes way past data!
xVals = sort([ -0.889487776317036, xVals ]); % must balance...
%xVals = sort([ -0.893393189225250, xVals ]);
%xVals = sort([ -0.890922389166615, xVals ]); % xB goes past...
%xVals = sort([ -0.890775831691018, xVals ]); % must balance...
%xVals = sort([ -0.891655176544600, xVals ]);
%xVals = sort([ -0.890874186962029, xVals ]); %1E-16.
%
%
%
numPts = size(xVals,2);
fVals = funchF(xVals);
%
[ foo, nC ] = min(fVals);
n = nC;
vecX = xVals(n-1:n+1)';
vecF = fVals(n-1:n+1)';
matX = [ ones(3,1), vecX, vecX.^2 ];
vecC = matX \ vecF;
funchFC = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
funchDFC = @(x)( vecC(2) + 2.0*vecC(3)*x );
xExtC = -vecC(2)/(2.0*vecC(3));
xAvgC = sum(vecX)/3.0;
%
if (0)
	% L = "Left", R = "Right"...
	nL = max([ nC-1, 2 ]);
	nR = min([ nC+1, numPts-1 ]);
else
	% L = "loser", R = "real"...
	if ( xExtC > xVals(nC) )
		nL = max([ nC-1, 2 ]);
		nR = min([ nC+1, numPts-1 ]);
	else
		nR = max([ nC-1, 2 ]);
		nL = min([ nC+1, numPts-1 ]);
	end
end
%
%
n = nL;
vecX = xVals(n-1:n+1)';
vecF = fVals(n-1:n+1)';
matX = [ ones(3,1), vecX, vecX.^2 ];
vecC = matX \ vecF;
funchFL = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
funchDFL = @(x)( vecC(2) + 2.0*vecC(3)*x );
xExtL = -vecC(2)/(2.0*vecC(3));
xAvgL = sum(vecX)/3.0;
%
n = nR;
vecX = xVals(n-1:n+1)';
vecF = fVals(n-1:n+1)';
matX = [ ones(3,1), vecX, vecX.^2 ];
vecC = matX \ vecF;
funchFR = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
funchDFR = @(x)( vecC(2) + 2.0*vecC(3)*x );
xExtR = -vecC(2)/(2.0*vecC(3));
xAvgR = sum(vecX)/3.0;
%
%
%
xVals = xVals(nC-2:nC+2);
fVals = funchF(xVals);
%
%
%
epsFD = sqrt(sqrt(eps));
funchDF = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
%
viz_numPts = 5000;
viz_xVals = linspace(min(xVals),max(xVals),viz_numPts);
%
viz_fValAll = [ funchF(viz_xVals), ...
  funchFL(viz_xVals), funchFC(viz_xVals), funchFR(viz_xVals), ...
  0.0 ];
viz_fValMin = min(viz_fValAll);
viz_fValMax = max(viz_fValAll);
%
%
%
numFigs++; figure(numFigs);
plot( ...
  xVals, fVals, 'ko', 'linewidth', 2, 'markersize', 25, ...
  viz_xVals, funchF(viz_xVals),  'o-', "color", [0.3,0.3,0.3], "markersize", 4, ...
  viz_xVals, funchFL(viz_xVals), '^-', "color", [0.8,0.0,0.0], "markersize", 4, ...
  viz_xVals, funchFC(viz_xVals), 's-', "color", [0.0,0.6,0.0], "markersize", 4, ...
  viz_xVals, funchFR(viz_xVals), 'v-', "color", [0.0,0.0,0.9], "markersize", 4, ...
  xExtL*[1,1], [viz_fValMin,viz_fValMax], '^-', "linewidth", 2, "markersize", 25, "color", [0.9,0.3,0.3], ...
  xExtC*[1,1], [viz_fValMin,viz_fValMax], 's-', "linewidth", 2, "markersize", 25, "color", [0.3,0.8,0.3], ...
  xExtR*[1,1], [viz_fValMin,viz_fValMax], 'v-', "linewidth", 2, "markersize", 25, "color", [0.3,0.3,1.0], ...
  x_secret*[1,1], [viz_fValMin,viz_fValMax], '*-', 'linewidth', 2, 'markersize', 25, "color", [0.8,0.6,0.0], ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
%
%
%
viz_dfValAll = [ funchDF(viz_xVals), ...
  funchDFL(viz_xVals), funchDFC(viz_xVals), funchDFR(viz_xVals), ...
  0.0 ];
viz_dfValMin = min(viz_dfValAll);
viz_dfValMax = max(viz_dfValAll);
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, funchDF(viz_xVals),  'o-', "color", [0.3,0.3,0.3], "markersize", 4, ...
  viz_xVals, funchDFL(viz_xVals), '^-', "color", [0.8,0.0,0.0], "markersize", 4, ...
  viz_xVals, funchDFC(viz_xVals), 's-', "color", [0.0,0.6,0.0], "markersize", 4, ...
  viz_xVals, funchDFR(viz_xVals), 'v-', "color", [0.0,0.0,0.9], "markersize", 4, ...
  xExtL*[1,1], [viz_dfValMin,viz_dfValMax], '^-', "linewidth", 2, "markersize", 25, "color", [0.9,0.3,0.3], ...
  xExtC*[1,1], [viz_dfValMin,viz_dfValMax], 's-', "linewidth", 2, "markersize", 25, "color", [0.3,0.8,0.3], ...
  xExtR*[1,1], [viz_dfValMin,viz_dfValMax], 'v-', "linewidth", 2, "markersize", 25, "color", [0.3,0.3,1.0], ...
  [xAvgC,xAvgR], [funchDFC(xAvgC),funchDFR(xAvgR)], 'x-', "linewidth", 2, "markersize", 25, "color", [0.3,0.8,0.9], ...
  x_secret*[1,1], [viz_dfValMin,viz_dfValMax], '*-', 'linewidth', 2, 'markersize', 25, "color", [0.8,0.6,0.0], ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "df" );
title( "df vs x" );
