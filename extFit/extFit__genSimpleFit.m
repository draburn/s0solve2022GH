function [ wVals, s0, p0, sMin, sMax, pMin, pMax, retCode, datOut ] = extFit__genSimpleFit( xVals, fVals, prm=[] )
	%
	commondefs;
	thisFile = "extFit__genSimpleFit";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	%
	numPts = size(xVals,2);
	if (doChecks)
		assert( 3 <= numPts );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( max(xVals) > min(xVals) );
		assert( max(fVals) > min(fVals) );
	end
	%
	wVals = ones(size(xVals)); % May be changed later.
	wVals /= numPts;
	s0 = ( max(xVals) + min(xVals) ) / 2.0; % May be changed later.
	p0 = 2.0;
	sMin = [];
	sMax = [];
	pMin = 0.0;
	pMax = [];
	retCode = RETCODE__SUCCESS; % Allow return at any time.
	datOut = [];
	%
	%
	%
	% First, consider a quad fit to EVERYTHING.
	vecX = xVals';
	vecF = fVals';
	matX = [ ones(numPts,1), vecX, vecX.^2 ];
	if ( rcond(matX'*matX) > eps300 )
		vecC = matX \ vecF;
		if ( abs(vecC(3)) > eps075 * abs(vecC(2)) )
			s0 = -vecC(2)/(2.0*vecC(3)); % Overwrite s0!
		end
	end
	%
	%
	% Look for a clean pt-wise ext...
	% Bail if any issues come up.
	%
	deltaF0 = fVals(2)-fVals(1);
	if ( 0.0 == deltaF0 )
		return;
	end
	%
	n = 2;
	while (1)
		if ( n>=numPts )
			% Didn't find a change in slope.
			return;
		end
		if ( deltaF0*(fVals(n+1)-fVals(n)) <= 0.0 )
			% Found first change in slope.
			break;
		end
		n++;
	end
	nC = n;
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "nC = %d.", nC ) );
	while (1)
		if ( n>=numPts )
			% Didn't find a second change in slope.
			break;
		end
		if ( deltaF0*(fVals(n+1)-fVals(n)) >= 0.0 )
			% Found a second change in slope.
			return;
		end
		n++;
	end
	assert( 2 <= nC );
	assert( nC <= numPts-1 );
	if ( xVals(nC-1) >= xVals(nC) || xVals(nC) >= xVals(nC+1) )
		% We're noping on this case.
		return;
	end
	%
	% So, looks like we have a clean pt-wise ext at nC!
	vecX = xVals(nC-1:nC+1)';
	vecF = fVals(nC-1:nC+1)';
	matX = [ ones(3,1), vecX, vecX.^2 ];
	vecC = matX \ vecF;
	%
	if ( abs(vecC(3)) < eps075*abs(vecC(2)) )
		return;
	end
	s0 = -vecC(2)/(2.0*vecC(3)); % Overwrite s0!
	bigF0 = vecC(1) - (vecC(2)^2)/(4.0*vecC(3));
	%
	gVals = fVals - bigF0;
	if (doChecks)
		assert( xVals(nC-1) <= s0 );
		assert( s0 <= xVals(nC+1) );
		numFValsOnWrongSideOfExt = sum( gVals(1)*gVals < 0.0 );
		assert( 0 == numFValsOnWrongSideOfExt );
	end
	fScale1 = max([ abs(bigF0-fVals(nC)), ...
	  min([ abs(fVals(nC-1)-fVals(nC)), abs(fVals(nC+1)-fVals(nC)) ]) ]);
	wVals = 1.0./( abs(gVals) + fScale1 + eps*max(abs(fVals)) );
	wVals /= sum(wVals);
	%
return;
end
%
