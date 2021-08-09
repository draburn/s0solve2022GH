function [ s, p, retCode, datOut ] = extFit( xVals, fVals, wVals=[], prm=[] )
	commondefs;
	thisFile = "extFit";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	datOut = [];
	%
	% The predominant hierarchy is:
	%  extFit ---> extFit__internal ---> extFit__genInitalGuess
	%                                \-> extFit__findFit
	% prm is passed down directly to __internal,
	%  but _genInitialGuess and _findFit are needed below that.
	% datOut is similarly passed back up.
	%
	msg( thisFile, __LINE__, "To-do..." );
	msg( thisFile, __LINE__, "  * Set s/p min/max when it makes sense to do so." );
	msg( thisFile, __LINE__, "  * Set wVals when it makes sense to do so." );
	msg( thisFile, __LINE__, "  * Broader testing." );
	msg( thisFile, __LINE__, "  * Replace 'doChecks' with varLev." );
	msg( thisFile, __LINE__, "Recently done..." );
	msg( thisFile, __LINE__, "  * Check normalization." );
	%
	numPts = size(xVals,2);
	if ( isempty(wVals) )
		wVals = ones(1,numPts);
		wVals /= sum(wVals);
	end
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	assert( isrealarray(wVals,[1,numPts]) );
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
	noWValIsNegative = (0==sum(wVals<0.0));
	assert( noWValIsNegative );
	wSum = sum(wVals);
	assert( wSum > 0.0 );
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
	wVals_sorted = wVals(sortingIndexes);
	%
	xVals_normalized = ( xVals_sorted - xMin ) / ( xMax - xMin );
	fVals_normalized = ( fVals_sorted - fMin ) / ( fMax - fMin );
	wVals_normalized = wVals_sorted / wSum;
	%
	% Use retCode and prm directly with __internal.
	[ s_normalized, p_normalized, retCode, datOut ] = ...
	  extFit__internal( xVals_normalized, fVals_normalized, wVals_normalized, prm );
	if ( retCode == RETCODE__SUCCESS
	  || retCode == RETCODE__IMPOSED_STOP
	  || retCode == RETCODE__ALGORITHM_BREAKDOWN )
		% These are the three return codes we'll treat as
		% more or less reasonable results.
		assert( isrealscalar(s_normalized) );
		assert( isrealscalar(p_normalized) );
		s = xMin + s_normalized*(xMax-xMin);
		p = p_normalized;
		[ rhoVals_normalized, ...
		  bigF0_normalized, ...
		  bigF1_normalized, ...
		  omega_normalized ] = ...
		 extFit__calcAtPt( ...
		   s_normalized, ...
		   p_normalized, ...
		   xVals_normalized, ...
		   fVals_normalized, ...
		   wVals_normalized, ...
		   prm );
		datOut.rhoVals = fMin + rhoVals_normalized*(fMax-fMin);
		datOut.bigF0 = fMin + bigF0_normalized*(fMax-fMin);
		datOut.bigF1 = bigF1_normalized*(fMax-fMin)/((xMax-xMin)^p);
		datOut.omega = omega_normalized * wSum;
	else
		s = (xMax+xMin)/2.0;
		p = 0.0;
	end
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "s = %g.", s ) );
	msg_copious( verbLev, thisFile, __LINE__, sprintf( "p = %g.", p ) );
	retCode = RETCODE__SUCCESS;
return;
end
