function [ s0, p0, retCode, datOut ] = extFit__genInitialGuess( xVals, fVals, wVals, prm=[] )
	commondefs;
	thisFile = "extFit__genInitialGuess";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	datOut = [];
	%
	%
	%
	if (1)
		msg( thisFile, __LINE__, "HACK! Forcing initial guess..." );
		p0 = 2.0;
		%s0 = xVals(2)+0.0236742;
		s0 = xVals(2)+0.0236743;
		retCode = RETCODE__SUCCESS;
		return
	end
	if (1)
		msg( thisFile, __LINE__, "HACK! Overiding wVals for initial guess..." );
		wVals(:) = 0.0;
		wVals(1:3) = 1.0;
	end
	%
	%
	%
	p0 = mygetfield( prm, "p0", 2.0 );
	assert( isrealscalar(p0) );
	assert( 0.0 < p0 );
	s0 = mygetfield( prm, "s0", [] );
	if ( ~isempty(s0) )
		assert( isrealscalar(s0) );
		retCode = RETCODE__SUCCESS;
		return;
	end
	%
	numPts = size(xVals,2);
	xMin = min(xVals);
	xMax = max(xVals);
	if ( doChecks )
		assert( 3 <= numPts );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealarray(wVals,[1,numPts]) );
		haveAtLeast3UniqueXVals = ( 0.0 < max(abs((xVals-xMin).*(xVals-xMax))) );
		assert( haveAtLeast3UniqueXVals );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	% Generate quadratic fit.
	vecX = xVals';
	vecF = fVals';
	matX = [ ones(numPts,1), vecX, vecX.^2 ];
	matW = diag(sqrt(wVals));
	vecC = (matW*matX)\(matW*vecF);
	if ( 0.0 == abs(vecC(3)) )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: Points are linear." );
		p = 0.0;
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	if ( abs(vecC(3)) <= eps075*abs(vecC(2)) )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: Points are nearly linear." );
	end
	s0 = -0.5*vecC(2)/vecC(3);
	assert( isrealscalar(s0) );
	if ( s0 < xMin - (xMax-xMin)/eps050 ...
	  || s0 > xMax + (xMax-xMin)/eps050 )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: s0 is very far from xVals" );
	end
	%
	retCode = RETCODE__SUCCESS;
return;
end
%
% DRaburn 2021.08.05:
% Another good way to approach extFit would be to look at f' / f'';
%  this is particularly true for generating the initial guess.
% However, from what I've seen, the calculate value near the exterma
%  tends to be inaccurate -- unless we go out of our way to use a small epislon,
%  which is not part of the "greedy" framework I'm foregrounding.
% And, well, KISS.
