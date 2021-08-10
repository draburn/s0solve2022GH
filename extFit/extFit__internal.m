function [ s, p, bigF0, bigF1, retCode, datOut ] = extFit__internal( xVals, fVals, wVals, prm=[] )
	commondefs;
	thisFile = "extFit__internal";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	datOut = [];
	%
	%
	%
	prm_genSimpleFit = mygetfield( prm, "prm_genSimpleFit", [] );
	[ wVals, s0, p0, sMin, sMax, pMin, pMax, retCode_genSimpleFit, datOut_genSimpleFit ] = ...
	  extFit__genSimpleFit( xVals, fVals, prm_genSimpleFit )
	datOut.datOut_genSimpleFit = datOut_genSimpleFit;
	if ( RETCODE__SUCCESS ~= retCode_genSimpleFit )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__genSimpleFit() returned %s.", retcode2str(retCode_genSimpleFit) ) );
		return;
	end
	%
	s0 = mygetfield( prm, "s0", s0 );
	p0 = mygetfield( prm, "p0", p0 );
	sMin = mygetfield( prm, "sMin", sMin );
	sMax = mygetfield( prm, "sMax", sMax );
	pMin = mygetfield( prm, "pMin", pMin );
	pMax = mygetfield( prm, "pMax", pMax );
	%
	numPts = size(xVals,2);
	if (doChecks)
		assert( 3 <= numPts );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( max(xVals) > min(xVals) );
		assert( max(fVals) > min(fVals) );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
		%
		assert( isrealscalar(s0) );
		if (~isempty(sMin))
			assert( s0 >= sMin );
		end
		if (~isempty(sMax))
			assert( s0 <= sMax );
		end
		%
		assert( isrealscalar(p0) );
		if (~isempty(pMin))
			assert( p0 >= pMin );
		end
		if (~isempty(pMax))
			assert( p0 <= pMax );
		end
	end
	%
	%
	%
	prm_findFit = mygetfield( prm, "prm_findFit", [] );
	prm_findFit.sMin = mygetfield( prm_findFit, "sMin", sMin );
	prm_findFit.sMax = mygetfield( prm_findFit, "sMax", sMax );
	prm_findFit.pMin = mygetfield( prm_findFit, "pMin", pMin );
	prm_findFit.pMax = mygetfield( prm_findFit, "pMax", pMax );
	[ s, p, bigF0, bigF1, retCode, datOut_findFit ] = extFit__findFit( ...
	  s0, p0, xVals, fVals, wVals, prm_findFit );
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
