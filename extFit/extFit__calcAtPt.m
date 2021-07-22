function [ rhoVals, bigF0, bigF1, omega ] = extFit__calcAtPt( ...
  s, p, xVals, fVals, nExactFit, wVals=[], prm=[] )
	%
	thisFile = "extFit__calcAtPt";
	doChecks = mygetfield( prm, "doChecks", true );
	%
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	if ( doChecks )
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
		assert( 0.0 < p );
		numPts = size(xVals,2);
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isposintscalar(nExactFit) );
		assert( 1 <= nExactFit );
		assert( nExactFit <= numPts );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	% Least-squares fit to f = F0 + F1 * | x - s |^p,
	% subject to an exact match at nExactFit.
	%
	yVals = abs( xVals - s ) .^p;
	dfVals = fVals - fVals(nExactFit);
	dyVals = yVals - yVals(nExactFit);
	%
	bigF1 = sum( wVals .* dyVals .* dfVals ) / sum( wVals .* dyVals.^2 );
	bigF0 = fVals(nExactFit) - bigF1*yVals(nExactFit);
	rhoVals = bigF0 + bigF1*yVals - fVals;
	omega = 0.5*sum( wVals .* rhoVals.^2 );
return;
end
