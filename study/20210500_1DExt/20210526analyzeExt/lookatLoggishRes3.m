clear;
commondefs;
thisFile = "lookatLoggishRes3";
numFigs = 0;
%
caseNum = 1;
switch (caseNum)
case 1
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x_secret = -0.890874387338734; % This is the left abs min.
	%%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
	xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843 ]);
	xVals = sort([ -0.831698038563058, xVals ]);
	xVals = sort([ -0.897620231966365, xVals ]);
	xVals = sort([ -0.867992251154230, xVals ]);
	xVals = sort([ -0.88284017185230, xVals ]);
	xVals = sort([ -0.887137697860992, xVals ]); % xB goes past data.
	xVals = sort([ -0.886078607798801, xVals ]);
	xVals = sort([ -0.890138678468405, xVals ]); % xB goes way past data!
	xVals = sort([ -0.889487776317036, xVals ]); % must balance...
	xVals = sort([ -0.893393189225250, xVals ]);
	%xVals = sort([ -0.890922389166615, xVals ]); % xB goes past...
	%xVals = sort([ -0.890775831691018, xVals ]); % must balance...
	%xVals = sort([ -0.891655176544600, xVals ]);
	%xVals = sort([ -0.890874186962029, xVals ]); %1E-16.
case 2
	x_secret = 0.8;
	funchF = @(x)( abs(x-x_secret).^6.5 );
	xVals = x_secret+1e-1*sort([ -5, -2, 1, 4, 7 ]);
case 3
	x_secret = 0.8;
	funchF = @(x)( abs(x-x_secret).^7.3 );
	%%%xVals = x_secret+linspace(-1.0,1.0,8)-0.3;
	xVals = x_secret+linspace(-1.0,1.0,12)-0.3;
case 4
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x_secret = -0.890874387338734; % This is the left abs min.
	xVals = x_secret + 1e-3*linspace(-1.0,1.0,8);
case 5
	funchF = @(x)( x + 1.0./x );
	x_secret = 1.0;
	xVals = 0.01+linspace(0.0,3.0,11);
case 6
	secret_bigA = 1.0;
	secret_bigB = 1.0;
	secret_bigC = 1.0;
	secret_bigS = 1.0;
	secret_bigP = 2.1;
	secret_foo = secret_bigB/(secret_bigC*secret_bigP);
	x_secret = secret_bigS - sign(secret_foo)*abs(secret_foo)^(1.0/(secret_bigP-1.0));
	funchF = @(x)( secret_bigA + secret_bigB*x + secret_bigC*abs(x-secret_bigS).^secret_bigP );
	xVals = x_secret+0.01+linspace(-1.0,1.0,11);
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
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
xAvg = sum(vecX)/3.0;
xSqAvg = sum(vecX.^2)/3.0;
xVar = sqrt( xSqAvg - xAvg^2 );
funchFC = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
funchDFC = @(x)( vecC(2) + 2.0*vecC(3)*x );
funchHC = @(x)( ( vecC(2) + 2.0*vecC(3)*x ) ./ abs(2.0*vecX(3)) );
xExtC = -vecC(2)/(2.0*vecC(3));
xAvgC = xAvg;
xVarC = xVar;
xExtCRes = xExtC-x_secret

%if ( xExtC > xVals(nC) )
%	extFit_xVals = xVals(nC-1:nC+2);
%	extFit_fVals = fVals(nC-1:nC+2);
%else
%	extFit_xVals = xVals(nC-2:nC+1);
%	extFit_fVals = fVals(nC-2:nC+1);
%end
%extFit_xVals = xVals(nC-1:nC+3);
%extFit_fVals = fVals(nC-1:nC+3);
extFit_xVals = xVals(nC-2:nC+2);
extFit_fVals = fVals(nC-2:nC+2);
%extFit_xVals = xVals(nC:nC+3);
%extFit_fVals = fVals(nC:nC+3);
%extFit_xVals = xVals(nC-3:nC);
%extFit_fVals = fVals(nC-3:nC);
%extFit_xVals = xVals(nC+1:nC+4);
%extFit_fVals = fVals(nC+1:nC+4);
extFit_wVals = 1.0./sqrt(sqrt(eps)+extFit_fVals-fVals(nC))
useNewExtFitModel = true;
extFitPrm = [];
%extFitPrm.verbLev = VERBLEV__COPIOUS;
extFitPrm.iterLimit = 100;
if (useNewExtFitModel)
%extFitDat = extFit( xVals(nC), 3.0, extFit_xVals, extFit_fVals, [], extFitPrm );
msg( thisFile, __LINE__, "..." );
extFitDat = extFit( xExtC, 3.0, extFit_xVals, extFit_fVals, extFit_wVals, extFitPrm );
msg( thisFile, __LINE__, "..." );
%extFitDat = extFit( secret_bigS, secret_bigP, extFit_xVals, extFit_fVals, extFit_wVals, extFitPrm );
[ extFit_omega, extFit_rho, extFit_bigA, extFit_bigB, extFit_bigC ] = extFit_calcOmega( ...
  extFit_xVals, extFit_fVals, extFitDat.bigX, extFitDat.bigP, extFit_wVals );
msg( thisFile, __LINE__, "..." );
echo__extFit_omega = extFit_omega
echo__extFit_rho = extFit_rho
echo__extFit_bigX = extFitDat.bigX
echo__extFit_bigP = extFitDat.bigP
echo__extFit_bigA = extFit_bigA
echo__extFit_bigB = extFit_bigB
echo__extFit_bigC = extFit_bigC
funchFExtFit = @(x)( extFit_bigA + extFit_bigB*x ...
  + extFit_bigC*abs(x-extFitDat.bigX).^extFitDat.bigP );
extFit_shiftyX = extFit_bigB / ( extFit_bigC * extFitDat.bigP )
xExtFit = extFitDat.bigX - sign(extFit_shiftyX)*abs(extFit_shiftyX)^(1.0/(extFitDat.bigP-1.0))
else
%extFitDat = extFit( xVals(nC), 3.0, extFit_xVals, extFit_fVals, [], extFitPrm );
extFitDat = extFit( xExtC, 2.0, extFit_xVals, extFit_fVals, extFit_wVals, extFitPrm );
[ extFit_omega, extFit_rho, extFit_bigA, extFit_bigB ] = extFit_calcOmega( ...
  extFit_xVals, extFit_fVals, extFitDat.bigX, extFitDat.bigP, extFit_wVals );
echo__extFit_omega = extFit_omega
echo__extFit_rho = extFit_rho
echo__extFit_bigA = extFit_bigA
echo__extFit_bigB = extFit_bigB
funchFExtFit = @(x)( extFit_bigA + extFit_bigB*abs(x-extFitDat.bigX).^extFitDat.bigP );
xExtFit = extFitDat.bigX
end
xExtFitRes = xExtFit-x_secret

%
extFitVizPrm = [];
extFitVizPrm.wVals = extFit_wVals;
extFit_deltaX = sqrt(eps)*(xVals(nC+1)-xVals(nC-1));
numFigs++; figure(numFigs);
extFit_viz( extFit_xVals, extFit_fVals, ...
  min(xVals), max(xVals), 1.5, 2.5 );
%  secret_bigS-1e0, secret_bigS+1e0, 1.70, 2.20, extFitVizPrm );
%  xVals(nC-1)+extFit_deltaX, xVals(nC+1)-extFit_deltaX, 1.5, 3.0 );
%hold on;
%plot( x_secret*[1,1], [1.0,10.0], 'r-' );
%hold off;


	msg( thisFile, __LINE__, "" );
	msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
	msg( thisFile, __LINE__, sprintf( "QuadC  residual = %g.", xExtC-x_secret ) );
	msg( thisFile, __LINE__, sprintf( "ExtFit residual = %g.", xExtFit-x_secret ) );
	msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	msg( thisFile, __LINE__, "" );

%
if (0)
	% L = "loser", R = "real"...
	if ( xExtC > xVals(nC) )
		nL = max([ nC-1, 2 ]);
		nR = min([ nC+1, numPts-1 ]);
	else
		nR = max([ nC-1, 2 ]);
		nL = min([ nC+1, numPts-1 ]);
	end
else
	% L = "Left", R = "Right"...
	nL = max([ nC-1, 2 ]);
	nR = min([ nC+1, numPts-1 ]);
end
%
%
n = nL;
vecX = xVals(n-1:n+1)';
vecF = fVals(n-1:n+1)';
matX = [ ones(3,1), vecX, vecX.^2 ];
vecC = matX \ vecF;
xAvg = sum(vecX)/3.0;
xSqAvg = sum(vecX.^2)/3.0;
xVar = sqrt( xSqAvg - xAvg^2 );
funchFL = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
funchDFL = @(x)( vecC(2) + 2.0*vecC(3)*x );
funchHL = @(x)( ( vecC(2) + 2.0*vecC(3)*x ) ./ abs(2.0*vecX(3)) );
xExtL = -vecC(2)/(2.0*vecC(3));
xAvgL = xAvg;
xVarL = xVar;
%
n = nR;
vecX = xVals(n-1:n+1)';
vecF = fVals(n-1:n+1)';
matX = [ ones(3,1), vecX, vecX.^2 ];
vecC = matX \ vecF;
xAvg = sum(vecX)/3.0;
xSqAvg = sum(vecX.^2)/3.0;
xVar = sqrt( xSqAvg - xAvg^2 );
funchFR = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
funchDFR = @(x)( vecC(2) + 2.0*vecC(3)*x );
funchHR = @(x)( ( vecC(2) + 2.0*vecC(3)*x ) ./ abs(2.0*vecX(3)) );
xExtR = -vecC(2)/(2.0*vecC(3));
xAvgR = xAvg;
xVarR = xVar;
%
%
%
%%%xVals = xVals(nC-2:nC+2);
%%%xVals = xVals(6:end-6);
fVals = funchF(xVals);
%
%
%
epsFD = sqrt(sqrt(eps));
funchDF = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
funchH = @(x)( funchDF(x)./abs(funchDDF(x)) );
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
  viz_xVals, funchFExtFit(viz_xVals), "x-", "color", [0.7,0.0,0.8], "markersize", 4, ...
  xExtL*[1,1], [viz_fValMin,viz_fValMax], '^-', "linewidth", 2, "markersize", 25, "color", [0.9,0.3,0.3], ...
  xExtC*[1,1], [viz_fValMin,viz_fValMax], 's-', "linewidth", 2, "markersize", 25, "color", [0.3,0.8,0.3], ...
  xExtR*[1,1], [viz_fValMin,viz_fValMax], 'v-', "linewidth", 2, "markersize", 25, "color", [0.3,0.3,1.0], ...
  xExtFit*[1,1], [viz_fValMin,viz_fValMax], 'x-', "linewidth", 2, "markersize", 25, "color", [0.9,0.3,1.0], ...
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
cxVals = cent(xVals);
dfVals = diff(fVals)./diff(xVals);
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, funchDF(viz_xVals),  'o-', "color", [0.3,0.3,0.3], "markersize", 4, ...
  viz_xVals, funchDFL(viz_xVals), '^-', "color", [0.8,0.0,0.0], "markersize", 4, ...
  viz_xVals, funchDFC(viz_xVals), 's-', "color", [0.0,0.6,0.0], "markersize", 4, ...
  viz_xVals, funchDFR(viz_xVals), 'v-', "color", [0.0,0.0,0.9], "markersize", 4, ...
  cxVals, dfVals,'+-', "color", [0.5,0.0,0.6], "markersize", 10, "linewidth", 1, ...
  xExtL*[1,1], [viz_dfValMin,viz_dfValMax], '^-', "linewidth", 2, "markersize", 25, "color", [0.9,0.3,0.3], ...
  xExtC*[1,1], [viz_dfValMin,viz_dfValMax], 's-', "linewidth", 2, "markersize", 25, "color", [0.3,0.8,0.3], ...
  xExtR*[1,1], [viz_dfValMin,viz_dfValMax], 'v-', "linewidth", 2, "markersize", 25, "color", [0.3,0.3,1.0], ...
  [xAvgL+xVarL,xAvgR-xVarR], [funchDFL(xAvgL+xVarL),funchDFR(xAvgR-xVarR)], 'x-', "linewidth", 5, "markersize", 30, "color", [0.3,0.8,0.9], ...
  x_secret*[1,1], [viz_dfValMin,viz_dfValMax], '*-', 'linewidth', 2, 'markersize', 25, "color", [0.8,0.6,0.0], ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "df" );
title( "df vs x" );
%
%
%
viz_hValAll = [ funchH(viz_xVals), ...
  funchHL(viz_xVals), funchHC(viz_xVals), funchHR(viz_xVals), ...
  0.0 ];
viz_hValMin = min(viz_hValAll);
viz_hValMax = max(viz_hValAll);
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, funchH(viz_xVals),  'o-', "color", [0.3,0.3,0.3], "markersize", 4, ...
  viz_xVals, funchHL(viz_xVals), '^-', "color", [0.8,0.0,0.0], "markersize", 4, ...
  viz_xVals, funchHC(viz_xVals), 's-', "color", [0.0,0.6,0.0], "markersize", 4, ...
  viz_xVals, funchHR(viz_xVals), 'v-', "color", [0.0,0.0,0.9], "markersize", 4, ...
  xExtL*[1,1], [viz_hValMin,viz_hValMax], '^-', "linewidth", 2, "markersize", 25, "color", [0.9,0.3,0.3], ...
  xExtC*[1,1], [viz_hValMin,viz_hValMax], 's-', "linewidth", 2, "markersize", 25, "color", [0.3,0.8,0.3], ...
  xExtR*[1,1], [viz_hValMin,viz_hValMax], 'v-', "linewidth", 2, "markersize", 25, "color", [0.3,0.3,1.0], ...
  [xAvgL-xVarL,xAvgR+xVarR], [funchHL(xAvgL-xVarL),funchHR(xAvgR+xVarR)], 'x-', "linewidth", 5, "markersize", 30, "color", [0.3,0.8,0.9], ...
  x_secret*[1,1], [viz_hValMin,viz_hValMax], '*-', 'linewidth', 2, 'markersize', 25, "color", [0.8,0.6,0.0], ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "h" );
title( "h vs x" );
