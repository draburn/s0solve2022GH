function [ wVals, sMin, sMax, pMin, pMax, retCode, datOut ] = extFit__genConstants( ...
  xVals, fVals, wVals, sMin, sMax, pMin, pMax, prm=[] )
	%
	% If we have a single pt-wise local ext, set default sMin and/or sMax based on that.
	% Allow for the (real-world extremely unlikely) case of two exactly equal values.
	commondefs;
	thisFile = "extFit__genDefaults";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	%
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	numPts = size(xVals,2);
	if (doChecks)
		assert( 3 <= numPts );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		if (~isempty(wVals))
			assert( isrealarray(wVals,[1,numPts]) );
			noWValIsNegative = (0==sum(wVals<0.0));
			assert( noWValIsNegative );
			atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
			assert( atLeastOneWValIsPositive );
		end
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
	end
	%
	% Determine if we have a clear single pt-wise ext;
	%  this is (currently) the only case in which we'll apply specific constants.
	% Let's not worry about what happens in the (real-world unlikely) case
	%  where two values are exactly equal.
	% If we encounter an exception to this expectation, just return.
	% In order to do that safely, we need all of our default return values set.
	wValsWASEmpty = isempty(wVals);
	if (wValsWASEmpty)
		wVals = ones(size(xVals)); % Unless changed later.
		wVals /= sum(wVals);
	end
	retCode = RETCODE__SUCCESS; % In a sense.
	%
	%
	% So, here goes!
	s0 = fVals(2)-fVals(1);
	if ( 0.0 == s0 )
		return;
	end
	%
	n = 2;
	while (1)
		if ( n>=numPts )
			% Didn't find a change in slope.
			return;
		end
		if ( s0*(fVals(n+1)-fVals(n)) <= 0.0 )
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
		if ( s0*(fVals(n+1)-fVals(n)) >= 0.0 )
			% Found a second change in slope.
			return;
		end
		n++;
	end
	assert( 2 <= nC );
	assert( nC <= numPts-1 );
	%
	% So, looks like we have a clean pt-wise ext at nC!
	%
	if ( isempty(sMin) && isempty(sMax) )
		sMin = xVals(nC-1);
		sMax = xVals(nC+1);
	elseif ( isempty(sMin) )
		sMin = xVals(nC-1); % Unless...
		if ( sMin >= sMax )
			sMin = [];
		end
	elseif ( isempty(sMax) )
		sMax = xVals(nC+1); % Unless...
		if ( sMax <= sMin )
			sMax = [];
		end
	end
	%
	if ( isempty(pMin) && isempty(pMax) )
		pMin = 0.0;
	elseif ( isempty(pMin) )
		pMin = 0.0; % Unless...
		if ( pMin >= pMax )
			pMin = [];
		end
	elseif ( isempty(pMax) )
		% Nothing to do.
	end
	%
	if (wValsWASEmpty)
		msg( thisFile, __LINE__, "SET wVals HERE!" );
		%wVals(4) = 0.0;
		%wVals /= sum(wVals);
		%assert(0);
	end
	%
return;
end
%
