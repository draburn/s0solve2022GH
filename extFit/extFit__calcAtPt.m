function [ rhoVals, bigF0, bigF1, omega ] = extFit__calcAtPt(
  s, p, xVals, fVals, wVals, prm=[] )
	%
	thisFile = "extFit__calcAtPt";
	doChecks = mygetfield( prm, "doChecks", true );
	%
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	numPts = size(xVals,2);
	if ( doChecks )
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
		assert( 0.0 < p );
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	% Least-squares fit to f = F0 + F1 * | x - s |^p.
	yVals = abs( xVals - s ).^p;
	if ( doChecks )
		assert( isrealarray(yVals,[1,numPts]) );
		% 47^209 = Inf, for example.
	end
	vecF = fVals';
	matY = [ ones(numPts,1), yVals' ];
	matW = diag(sqrt(wVals));
	vecC = (matW*matY)\(matW*vecF);
	bigF0 = vecC(1);
	bigF1 = vecC(2);
	rhoVals = bigF0 + bigF1*yVals - fVals;
	omega = 0.5*sum( wVals .* rhoVals.^2 );
return;
end
	%
	% DRaburn 2021.08.02.
	% This is an older version that had one point fit exactly.
	%
	% Least-squares fit to f = F0 + F1 * | x - s |^p,
	% subject to an exact match at nExactFit.
	%
	%yVals = abs( xVals - s ) .^p;
	%dfVals = fVals - fVals(nExactFit);
	%dyVals = yVals - yVals(nExactFit);
	%%
	%bigF1 = sum( wVals .* dyVals .* dfVals ) / sum( wVals .* dyVals.^2 );
	%bigF0 = fVals(nExactFit) - bigF1*yVals(nExactFit);
	%rhoVals = bigF0 + bigF1*yVals - fVals;
	%omega = 0.5*sum( wVals .* rhoVals.^2 );
%return;
%end
