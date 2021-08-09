function [ s, p, bigF0, bigF1, retCode, datOut ] = extFit__internal( xVals, fVals, wVals, prm=[] )
	commondefs;
	thisFile = "extFit__internal";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	datOut = [];
	%
	sMin = mygetfield( prm, "sMin", [] );
	sMax = mygetfield( prm, "sMax", [] );
	pMin = mygetfield( prm, "pMin", [] );
	pMax = mygetfield( prm, "pMax", [] );
	%
	prm_genConstants = mygetfield( prm, "prm_genConstants", [] );
	[ wVals, sMin, sMax, pMin, pMax, retCode, prm_genParams ] = extFit__genConstants( ...
	  xVals, fVals, wVals, sMin, sMax, pMin, pMax, prm_genConstants );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__genConstants() failed for initial guess (%d).", retCode ) );
		return;
	end
	if (doChecks)
		if (~isempty(sMin))
			assert( isrealscalar(sMin) );
		end
		if (~isempty(sMax))
			assert( isrealscalar(sMax) );
		end
		if ( ~isempty(sMin) && ~isempty(sMax) )
			assert( sMin < sMax );
		end
		if (~isempty(pMin))
			assert( isrealscalar(pMin) );
		end
		if (~isempty(pMax))
			assert( isrealscalar(pMax) );
		end
		if ( ~isempty(pMin) && ~isempty(pMax) )
			assert( pMin < pMax );
		end
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	%
	%
	if (1)
		msg( thisFile, __LINE__, "HACK! Overiding sMin and Max..." );
		sMin = [];
		sMax = [];
	end
	%
	if ( verbLev >= VERBLEV__COPIOUS )
		msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
		echo__xVals = xVals
		echo__fVals = fVals
		echo__wVals = wVals
		echo__sMin = sMin
		echo__sMax = sMax
		echo__pMin = pMin
		echo__pMax = pMax
		plot( xVals, fVals, 'o-' );
		grid on;
		msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	end
	%
	%
	%
	prm_genInitialGuess = mygetfield( prm, "prm_genInitialGuess", [] );
	prm_genInitialGuess.sMin = mygetfield( prm_genInitialGuess, "sMin", sMin );
	prm_genInitialGuess.sMax = mygetfield( prm_genInitialGuess, "sMax", sMax );
	prm_genInitialGuess.pMin = mygetfield( prm_genInitialGuess, "pMin", pMin );
	prm_genInitialGuess.pMax = mygetfield( prm_genInitialGuess, "pMax", pMax );
	[ s0, p0, retCode, datOut_genInitialGuess ] = extFit__genInitialGuess( ...
	  xVals, fVals, wVals, prm_genInitialGuess );
	datOut.datOut_genInitialGuess = datOut_genInitialGuess;
	if ( retCode ~= RETCODE__SUCCESS )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "WARNING: __genInitialGuess returned %s.", retcode2str(retCode) ) );
		s = s0;
		p = p0;
		return;
	end
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "s0 = %g.", s0 ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "p0 = %g.", p0 ) );
	if (doChecks)
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		if (~isempty(sMin))
			assert( sMin <= s0 );
		end
		if (~isempty(sMax))
			assert( s0 <= sMax );
		end
		if (~isempty(pMin))
			assert( pMin <= p0 );
		end
		if (~isempty(pMax))
			assert( p0 <= pMax );
		end
	end
	%
	prm_findFit = mygetfield( prm, "prm_findFit", [] );
	prm_findFit.sMin = mygetfield( prm_findFit, "sMin", sMin );
	prm_findFit.sMax = mygetfield( prm_findFit, "sMax", sMax );
	prm_findFit.pMin = mygetfield( prm_findFit, "pMin", pMin );
	prm_findFit.pMax = mygetfield( prm_findFit, "pMax", pMax );
	[ s, p, bigF0, bigF1, retCode, datOut_findFit ] = extFit__findFit( s0, p0, xVals, fVals, wVals, prm_findFit );
	datOut.datOut_findFit = datOut_findFit;
	if ( retCode ~= RETCODE__SUCCESS
	  && retCode ~= RETCODE__IMPOSED_STOP
	  && retCode ~= RETCODE__IMPOSED_STOP )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "WARNING: __findFit returned %s.", retcode2str(retCode) ) );
	end
	if (doChecks)
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
		if (~isempty(sMin))
			assert( sMin <= s );
		end
		if (~isempty(sMax))
			assert( s <= sMax );
		end
		if (~isempty(pMin))
			assert( pMin <= p );
		end
		if (~isempty(pMax))
			assert( p <= pMax );
		end
	end
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "s = %g.", s ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "p = %g.", p ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "bigF0 = %g.", bigF0 ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "bigF1 = %g.", bigF1 ) );
return;
end
