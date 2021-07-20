thisFile = "groot1d__getXNew";

numPts = size(xVals_raw,2);
assert( fevalCount == numPts );
if ( 0 == numPts )
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'given - x1'." );
	xNew = x1;
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end
if ( 1 == numPts )
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'given - x2'." );
	xNew = x2;
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end

assert( 2 <= fevalCount );
xMin = min(xVals_raw);
xMax = max(xVals_raw);
fMin = min(fVals_raw);
fMax = max(fVals_raw);

% Handle typical case, including bounded.
xNew = findGoodCand( xVals_sorted, fVals_sorted );
if ( isrealscalar(xNew) )
	xCapHi = xVals_sorted(end) + (xVals_sorted(end)-xVals_sorted(1));
	xCapLo = xVals_sorted(1)   - (xVals_sorted(end)-xVals_sorted(1));
	xNew = cap( xNew, xCapLo, xCapHi );
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end

% Handle "cylce: expand left, expand right, bisect"?
% Handle "all values are the same".
%%% HOW SHOULD THIS COMPETE WITH LARGER-SCALE TYPICAL CASE???
m = mod(desperationIter,5); desperationIter++;
if (0==m)
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - expand right'." );
	xNew = xMax + 0.3*(xMax-xMin);
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
elseif (1==m)
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - expand left'." );
	xNew = xMin - 0.3*(xMax-xMin);
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end
%
msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - bisect largest gap'." );
[ largestGapSize, largestGapIndex ] = max(diff(xVals_sorted));
xNew = (xVals_sorted(largestGapIndex) + xVals_sorted(largestGapIndex+1))/2.0;
thisFile = "RETURNING FROM groot1d__getXNew";
return;
