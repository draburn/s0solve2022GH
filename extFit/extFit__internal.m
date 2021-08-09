function [ s, p, retCode, datOut ] = extFit__internal( xVals, fVals, wVals, prm=[] )
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
	prm_genInitialGuess = mygetfield( prm, "prm_genInitialGuess", [] );
	prm_genInitialGuess.sMin = mygetfield( prm_genInitialGuess, "sMin", sMin );
	prm_genInitialGuess.sMax = mygetfield( prm_genInitialGuess, "sMax", sMax );
	prm_genInitialGuess.pMin = mygetfield( prm_genInitialGuess, "pMin", pMin );
	prm_genInitialGuess.pMax = mygetfield( prm_genInitialGuess, "pMax", pMax );
	[ s0, p0, retCode, datOut_genInitialGuess ] = extFit__genInitialGuess( xVals, fVals, wVals, prm_genInitialGuess );
	datOut.datOut_genInitialGuess = datOut_genInitialGuess;
	if ( retCode ~= RETCODE__SUCCESS )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "WARNING: __genInitialGuess returned %s.", retcode2str(retCode) ) );
		s = s0;
		p = p0;
		return;
	end
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
	else
	end
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "s0 = %g.", s0 ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "p0 = %g.", p0 ) );
	%
	prm_findFit = mygetfield( prm, "prm_findFit", [] );
	prm_findFit.sMin = mygetfield( prm_findFit, "sMin", sMin );
	prm_findFit.sMax = mygetfield( prm_findFit, "sMax", sMax );
	prm_findFit.pMin = mygetfield( prm_findFit, "pMin", pMin );
	prm_findFit.pMax = mygetfield( prm_findFit, "pMax", pMax );
	[ s, p, retCode, datOut_findFit ] = extFit__findFit( s0, p0, xVals, fVals, wVals, prm_findFit );
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
		assert( 0.0 < p );
	end
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "s = %g.", s ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "p = %g.", p ) );
return;
end
