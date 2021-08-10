clear;
commondefs;
thisFile = "test_extFit";
numFigs = 0;
%
hadBad = false;
numTests = 100;
for n=1:numTests
	setprngstates(n);
	bigF0_secret = randn();
	bigF1_secret = randn();
	s_secret = randn();
	p_secret = 1.0 + abs(randn());
	funchF = @(x)( bigF0_secret + bigF1_secret * abs( x - s_secret ).^p_secret );
	%
	numPts = 5+abs(round(4*randn()));
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
	funchFModel = @(x)( bigF0 + bigF1*abs(x-s).^p );
	fModelVals = funchFModel(xVals);
	rhoVals = fModelVals - fVals;
	myOmega = 0.5*sum(rhoVals.^2);
	myTol = 0.5*eps100*sum(fVals.^2);
	if ( myOmega > myTol )
		msg( thisFile, __LINE__, "" );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, sprintf( "Bad fit (%g vs %g) on test %d.", myOmega, myTol, n ) );
		msg( thisFile, __LINE__, "***" );
		msg( thisFile, __LINE__, "" );
		hadBad = true;
		break;
	end
	if ( ~fleq(s,s_secret,eps025) || ~fleq(bigF0,bigF0_secret,eps025)...
	  || ~fleq(p,p_secret,eps025) || ~fleq(bigF1,bigF1_secret,eps025) );
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
	return;
end
%
msg( thisFile, __LINE__, sprintf( "Redoing test %d....", n ) );
msg( thisFile, __LINE__, sprintf( "numPts = %d.", numPts ) );
%echo__xVals = xVals
%echo__fVals = fVals
wVals = [];
prm = [];
prm.verbLev = VERBLEV__COPIOUS;
prm.prm_findFit.verbLev = VERBLEV__COPIOUS;
[ s, p, bigF0, bigF1, retCode, datOut ] = extFit( xVals, fVals, wVals, prm );
msg( thisFile, __LINE__, sprintf( "extFit returned %s.", retcode2str(retCode) ) );
if ( retCode ~= RETCODE__SUCCESS )
	return;
end
msg( thisFile, __LINE__, "" );
msg( thisFile, __LINE__, "***" );
funchFModel = @(x)( bigF0 + bigF1*abs(x-s).^p );
fModelVals = funchFModel(xVals);
rhoVals = fModelVals - fVals;
myOmega = 0.5*sum(rhoVals.^2);
myTol = 0.5*eps100*sum(fVals.^2);
%echo__fModelVals = fModelVals
%echo__rhoVals = rhoVals
msg( thisFile, __LINE__, sprintf( "quant:   true value,   fit result,   residual." ) );
msg( thisFile, __LINE__, sprintf( "p:      %11.3e,  %11.3e,  %11.3e.", p_secret,     p,     p-p_secret ) );
msg( thisFile, __LINE__, sprintf( "s:      %11.3e,  %11.3e,  %11.3e.", s_secret,     s,     s-s_secret ) );
msg( thisFile, __LINE__, sprintf( "bigF0:  %11.3e,  %11.3e,  %11.3e.", bigF0_secret, bigF0, bigF0-bigF0_secret ) );
msg( thisFile, __LINE__, sprintf( "bigF1:  %11.3e,  %11.3e,  %11.3e.", bigF1_secret, bigF1, bigF1-bigF1_secret ) );
msg( thisFile, __LINE__, sprintf( "total residual:  %11.3e  ( vs %11.3e,  %11.3e ).", myOmega, myTol, myOmega/myTol ) );
%
viz_numPts = 1000;
viz_xVals = linspace(min(xVals),max(xVals),viz_numPts);
viz_fModelVals = funchFModel(viz_xVals);
plot( ...
  viz_xVals, viz_fModelVals, 'o-', 'linewidth', 2, ...
  xVals, fVals, 'kx', 'linewidth', 3, 'markersize', 20 );
grid on;
legend( ...
  "fModel", ...
  "original points", ...
  "location", "north" );
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
