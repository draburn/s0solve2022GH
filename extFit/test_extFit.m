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
	if ( ~fleq(s,s_secret,eps025) || ~ fleq(bigF0,bigF0_secret,eps025)...
	  || ~fleq(p,p_secret,eps025) || ~ fleq(bigF1,bigF1_secret,eps025) );
		msg( thisFile, __LINE__, sprintf( "Bad results on test %d.", n ) );
		hadBad = true;
		break;
	end
end
if (~hadBad)
	msg( thisFile, __LINE__, sprintf( "Found a good result for all %d tests.", numTests ) );
	return;
end
