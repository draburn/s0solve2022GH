function [ omega, bigF0, bigF1 ] = extFit_getOmega( s, p, xVals, fVals, nFit, wVals=[], prm=[] )
	thisFile = "extFit_getOmega";
	doChecks_default = true;
	%
	if ( mygetfield( prm, "doChecks", doChecks_default ) )
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
		numPts = size(xVals,2);
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealscalar(nFit) );
		assert( fleq(nFit,round(nFit)) );
		assert( 1 <= nFit );
		assert( nFit <= numPts );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
	end
	%
	xFit = xVals(nFit);
	fFit = fVals(nFit);
	gFit = abs( xFit - s ).^p;
	gVals = abs( xVals - s).^p;
	cVals = gVals - gFit;
	dVals = fVals - fFit;
	%
	if ( isempty(wVals) )
		bigF1 = sum( cVals .* dVals ) / sum( cVals .* cVals );
		bigF0 = fFit - bigF1*gFit;
		omega = 0.5 * sum( ( bigF0 + bigF1*gVals - fVals ).^2 );
		return;
	end
	%
	if ( mygetfield( prm, "doChecks", doChecks_default ) )
		assert( isrealarray(wVals,[1,numPts]) );
		noWValsAreNegative = (0==sum(wVals<0.0));
		assert( noWValsAreNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	bigF1 = sum( wVals .* cVals .* dVals ) / sum( wVals .* cVals .* cVals );
	bigF0 = fFit - bigF1*gFit;
	omega = 0.5 * sum( wVals .* ( bigF0 + bigF1*gVals - fVals ).^2 );
return;
end
