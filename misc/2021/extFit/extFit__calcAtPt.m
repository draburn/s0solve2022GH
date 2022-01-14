function [ rhoVals, bigF0, bigF1, omega, retCode ] = extFit__calcAtPt(
  s, p, xVals, fVals, wVals, prm=[] )
	%
	commondefs;
	thisFile = "extFit__calcAtPt";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	doChecks = mygetfield( prm, "doChecks", false );
	%
	numPts = size(xVals,2);
	if ( doChecks )
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
		%assert( 0.0 < p );
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
	rhoVals = zeros(1,numPts);
	bigF0 = 0.0;
	bigF1 = 0.0;
	omega = -1.0;
	%
	% Least-squares fit to f = F0 + F1 * | x - s |^p.
	yVals = abs( xVals - s ).^p;
	if ( ~isrealarray(yVals,[1,numPts]) )
		% 47^209 = Inf, for example.
		% Or, maybe p is negative.
		if ( verbLev >= VERBLEV__COPIOUS )
			msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
			echo__yVals = yVals
			msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
		end
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	vecF = fVals';
	matY = [ ones(numPts,1), yVals' ];
	matW = diag(sqrt(wVals));
	matA = matW*matY;
	if ( eps >= rcond(matA'*matA) )
		if ( verbLev >= VERBLEV__COPIOUS )
			msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
			echo__matY = matY
			echo__matW = matW
			echo__matA = matA
			msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
		end
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	vecC = matA\(matW*vecF);
	bigF0 = vecC(1);
	bigF1 = vecC(2);
	rhoVals = bigF0 + bigF1*yVals - fVals;
	omega = 0.5*sum( wVals .* rhoVals.^2 );
	retCode = RETCODE__SUCCESS;
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
