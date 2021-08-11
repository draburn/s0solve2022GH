function [ s, p, bigF0, bigF1, retCode, datOut ] = extFit( xVals, fVals, wVals=[], prm=[] )
	commondefs;
	thisFile = "extFit";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	datOut = [];
	%
	% prm is passed down directly to __internal,
	%  but _genInitialGuess and _findFit are needed below that.
	% datOut is similarly passed back up.
	%
	%
	numPts = size(xVals,2);
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	if (~isempty(wVals) )
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		wSum = sum(wVals);
		assert( wSum > 0.0 );
	end
	%
	xMax = max(xVals);
	xMin = min(xVals);
	fMax = max(fVals);
	fMin = min(fVals);
	if ( xMax == xMin )
		s = (xMax+xMin)/2.0;
		p = 0.0;
		msg_warn( verbLev, thisFile, __LINE__, ...
		  "WARNING: xVals are constant." );
		retCode = RETCODE__BAD_INPUT;
		datOut = [];
		return;
	end
	if ( fMax == fMin )
		s = (xMax+xMin)/2.0;
		p = 0.0;
		msg_warn( verbLev, thisFile, __LINE__, ...
		  "WARNING: fVals are constant." );
		retCode = RETCODE__BAD_INPUT;
		datOut = [];
		return;
	end
	%
	if ( xMax <= xMin + eps050*(abs(xMax)+abs(xMin)) )
		msg_warn( verbLev, thisFile, __LINE__, ...
		  "WARNING: Fractional variation of xVals is small." );
	end
	if ( fMax <= fMin + eps050*(abs(fMax)+abs(fMin)) )
		msg_warn( verbLev, thisFile, __LINE__, ...
		  "WARNING: Fractional variation of fVals is small." );
	end
	%
	% Sort and normalize data.
	[ xVals_sorted, sortingIndexes ] = sort( xVals );
	fVals_sorted = fVals(sortingIndexes);
	if (~isempty(wVals))
		wVals_sorted = wVals(sortingIndexes);
	else
		wVals_sorted = [];
	end
	%
	xVals_normalized = ( xVals_sorted - xMin ) / ( xMax - xMin );
	fVals_normalized = ( fVals_sorted - fMin ) / ( fMax - fMin );
	if (~isempty(wVals_sorted))
		wVals_normalized = wVals_sorted / wSum;
	else
		wVals_normalized = [];
	end
	%
	% Use retCode and prm directly with __internal.
	[ s_normalized, p_normalized, bigF0_normalized, bigF1_normalized, retCode, datOut ] = ...
	  extFit__internal( xVals_normalized, fVals_normalized, wVals_normalized, prm );
	if ( retCode == RETCODE__SUCCESS
	  || retCode == RETCODE__IMPOSED_STOP
	  || retCode == RETCODE__ALGORITHM_BREAKDOWN )
		% These are the three return codes we'll treat as
		% more or less reasonable results.
		assert( isrealscalar(s_normalized) );
		assert( isrealscalar(p_normalized) );
		assert( isrealscalar(bigF0_normalized) );
		assert( isrealscalar(bigF1_normalized) );
		s = xMin + s_normalized*(xMax-xMin);
		p = p_normalized;
		bigF0 = fMin + bigF0_normalized*(fMax-fMin);
		bigF1 = bigF1_normalized*(fMax-fMin)/((xMax-xMin)^p);
		msg_copious( verbLev, thisFile, __LINE__, sprintf( "s = %g.", s ) );
		msg_copious( verbLev, thisFile, __LINE__, sprintf( "p = %g.", p ) );
		msg_copious( verbLev, thisFile, __LINE__, sprintf( "bigF0 = %g.", bigF0 ) );
		msg_copious( verbLev, thisFile, __LINE__, sprintf( "bigF1 = %g.", bigF1 ) );
		retCode = RETCODE__SUCCESS;
	else
		s = [];
		p = [];
		bigF0 = [];
		bigF1 = [];
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
	end
return;
end
