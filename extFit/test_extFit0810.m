clear;
commondefs;
thisFile = "test_extFit";
numFigs = 0;
%
tic();
hadBad = false;
numTests = 100;
seed0 = 0%mod( round(now*1E11), 1E8 )
%
% These cases have issue using p_secret = 1.0 + abs(randn());
%  because p is close to unity.
%setprngstates(14767083)  % Fails mildly on test 1.
%setprngstates(328669); % p is nearly 1.0.
%setprngstates(49002736); % Fails hard.
setprngstates(91175170); % Has problems. Also after 4th test.
for n=1:numTests
	%setprngstates(seed0+n);
	bigF0_secret = randn();
	bigF1_secret = randn();
	s_secret = randn();
	p_secret = 1.5 + abs(randn());
	funchF = @(x)( bigF0_secret + bigF1_secret * abs( x - s_secret ).^p_secret );
	%
	numPts = 4+abs(round(4*randn()));
	if ( randn()>0.0 )
		r0 = randn();
		xVals = r0+sign(r0)*abs(randn(1,numPts));
	else
		xVals = randn(1,numPts);
	end
	fVals = funchF(xVals);
	%
	[ s, p, bigF0, bigF1, retCode, datOut ] = extFit( xVals, fVals );
	if ( retCode ~= RETCODE__SUCCESS )
		msg( thisFile, __LINE__, "" );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, sprintf( "extFit returned %s on test %d.", retcode2str(retCode), n ) );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, "" );
		hadBad = true;
		break;
	end
	if ( abs(s-s_secret) > eps025*(max(xVals)-min(xVals)) ...
	  || abs(bigF0-bigF0_secret) > eps025 )
		msg( thisFile, __LINE__, sprintf( "Fit for test %d is not great (%g,%g;%g,%g).", ...
		   n, s_secret, s, bigF0_secret, bigF0 ) );
	end
	funchFModel = @(x)( bigF0 + bigF1*abs(x-s).^p );
	fModelVals = funchFModel(xVals);
	rhoVals = fModelVals - fVals;
	myOmega = 0.5*sum(rhoVals.^2);
	myTol = 0.5*eps075*sum(fVals.^2);
	if ( myOmega > myTol )
		msg( thisFile, __LINE__, "" );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, sprintf( "Bad fit (%g vs %g) on test %d.", myOmega, myTol, n ) );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, "" );
		hadBad = true;
		break;
	end
	if ( abs(s-s_secret) > 0.1*(max(xVals)-min(xVals)) ...
	  || abs(p-p_secret) > 0.1 ...
	  || abs(bigF0-bigF0_secret) > 0.1 ...
	  || abs(bigF1-bigF1_secret) > 0.1 )
		msg( thisFile, __LINE__, "" );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, sprintf( "Incorrect fit on test %d.", n ) );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, "" );
		hadBad = true;
		break;
	end
end
if (~hadBad)
	msg( thisFile, __LINE__, sprintf( "Found a good result for all %d tests.", numTests ) );
	toc();
	return;
end
toc();
%
msg( thisFile, __LINE__, sprintf( "Redoing test %d....", n ) );
msg( thisFile, __LINE__, sprintf( "numPts = %d.", numPts ) );
%echo__xVals = xVals
%echo__fVals = fVals
wVals = [];
prm = [];
prm.verbLev = VERBLEV__COPIOUS;
prm.prm_findFit.verbLev = VERBLEV__COPIOUS;
prm.prm_findFit.prm_findStep.verbLev = VERBLEV__COPIOUS;
prm.prm_genSimpleFit = VERBLEV__COPIOUS;
[ s, p, bigF0, bigF1, retCode, datOut ] = extFit( xVals, fVals, wVals, prm );
msg( thisFile, __LINE__, "" );
msg( thisFile, __LINE__, "***" );
msg( thisFile, __LINE__, sprintf( "extFit returned %s.", retcode2str(retCode) ) );
msg( thisFile, __LINE__, "" );
echo__xVals = xVals
echo__fVals = fVals
msg( thisFile, __LINE__, sprintf( "numPts = %d.", numPts ) );
viz_numPts = 1000;
viz_xLo = min([ min(xVals), s_secret, s ]);
viz_xHi = max([ max(xVals), s_secret, s ]);
viz_xVals = linspace(viz_xLo,viz_xHi,viz_numPts);
if ( retCode ~= RETCODE__SUCCESS )
	msg( thisFile, __LINE__, sprintf( "s_secret:     %11.3e.", s_secret ) );
	msg( thisFile, __LINE__, sprintf( "p_secret:     %11.3e.", p_secret ) );
	msg( thisFile, __LINE__, sprintf( "bigF0_secret: %11.3e.", bigF0_secret ) );
	msg( thisFile, __LINE__, sprintf( "bigF1_secret: %11.3e.", bigF1_secret ) );
	numFigs++; figure(numFigs);
	plot( ...
	  viz_xVals, funchF(viz_xVals), 'k-', ...
	  s_secret, funchF(s_secret), 'r+', 'linewidth', 2, 'markersize', 25, ...
	  s_secret, bigF0_secret, 'ro', 'linewidth', 2, 'markersize', 25, ...
	  xVals, funchF(xVals), 'kx', 'linewidth', 3, 'markersize', 20, ...
	  xVals, fVals, 'ko', 'linewidth', 3, 'markersize', 20 );
	grid on;
	legend( ...
	  "actual f", ...
	  "actual s", ...
	  "actual s", ...
	  "data pts", ...
	  "data pts", ...
	  "location", "northWest" );
	xlabel( "x" );
	ylabel( "f" );
	title( "f vs x" );
	return;
end
funchFModel = @(x)( bigF0 + bigF1*abs(x-s).^p );
fModelVals = funchFModel(xVals);
rhoVals = fModelVals - fVals;
myOmega = 0.5*sum(rhoVals.^2);
myTol = 0.5*eps100*sum(fVals.^2);
%echo__fModelVals = fModelVals
%echo__rhoVals = rhoVals
msg( thisFile, __LINE__, sprintf( "quant:   true value,   fit result,   residual." ) );
msg( thisFile, __LINE__, sprintf( "s:      %11.3e,  %11.3e,  %11.3e.", s_secret,     s,     s-s_secret ) );
msg( thisFile, __LINE__, sprintf( "p:      %11.3e,  %11.3e,  %11.3e.", p_secret,     p,     p-p_secret ) );
msg( thisFile, __LINE__, sprintf( "bigF0:  %11.3e,  %11.3e,  %11.3e.", bigF0_secret, bigF0, bigF0-bigF0_secret ) );
msg( thisFile, __LINE__, sprintf( "bigF1:  %11.3e,  %11.3e,  %11.3e.", bigF1_secret, bigF1, bigF1-bigF1_secret ) );
msg( thisFile, __LINE__, sprintf( "total residual:  %11.3e  ( vs %11.3e,  %11.3e ).", myOmega, myTol, myOmega/myTol ) );
msg( thisFile, __LINE__, "***" );
msg( thisFile, __LINE__, "" );
%
viz_fModelVals = funchFModel(viz_xVals);
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_fModelVals, 'o-', 'linewidth', 2, ...
  s, funchFModel(s), 'g+', 'linewidth', 2, 'markersize', 25, ...
  s, bigF0, 'go', 'linewidth', 2, 'markersize', 25, ...
  viz_xVals, funchF(viz_xVals), 'k-', ...
  s_secret, funchF(s_secret), 'r+', 'linewidth', 2, 'markersize', 25, ...
  s_secret, bigF0_secret, 'ro', 'linewidth', 2, 'markersize', 25, ...
  xVals, funchF(xVals), 'kx', 'linewidth', 3, 'markersize', 20, ...
  xVals, fVals, 'ko', 'linewidth', 3, 'markersize', 20 );
grid on;
legend( ...
  "fModel", ...
  "sModel", ...
  "sModel", ...
  "actual f", ...
  "actual s", ...
  "actual s", ...
  "data pts", ...
  "data pts", ...
  "location", "northWest" );
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
