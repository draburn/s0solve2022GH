thisFile = "groot1d__getXNew";

if ( 1 == fevalCount )
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'x2'." );
	xNew = x2;
	thisFile = "RETURNED FROM groot1d__getXNew";
	return;
end

assert( 2 <= fevalCount );
xMin = min(xVals);
xMax = max(xVals);
fMin = min(fVals);
fMax = max(fVals);

% Handle "bounded".
if ( fMin * fMax < 0.0 )
	n = 1;
	while ( fVals_sorted(n)*fVals_sorted(n+1) > 0.0 )
		n++;
	end
	assert( fVals_sorted(n)*fVals_sorted(n+1) < 0.0 );
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("n = %d.",n) );
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("boundedIter = %d.",boundedIter) );
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("xVals_sorted(n)   = %f.", xVals_sorted(n)  ) );
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("xVals_sorted(n+1) = %f.", xVals_sorted(n+1)) );
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("fVals_sorted(n)   = %f.", fVals_sorted(n)  ) );
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("fVals_sorted(n+1) = %f.", fVals_sorted(n+1)) );
	%
	% For simplicity, just alternate between lin and bisection.
	% More advanced methods would include using a cubic
	%  (after checking for good agreement with points further out),
	%  and using an "reduce distance between bounding points"
	%  that is smarter than just bisection.
	% But, all of this is irrelevant to higher dim anyway.
	m = mod(boundedIter,2); boundedIter++;
	if ( 0==m )
		msg_copious( verbLev, thisFile, __LINE__, "xNew via 'bounded - linear interpolation'." );
		xNew = xVals_sorted(n) - fVals_sorted(n) ...
		 * ( xVals_sorted(n)-xVals_sorted(n+1) ) ...
		 / ( fVals_sorted(n)-fVals_sorted(n+1) );
	else
		msg_copious( verbLev, thisFile, __LINE__, "xNew via 'bounded - bisection'." );
		xNew = ( xVals_sorted(n) + xVals_sorted(n+1) )/2.0;
	end
	%msg_copious( verbLev, thisFile, __LINE__, sprintf("xNew = %f.",xNew) );
	thisFile = "RETURNED FROM groot1d__getXNew";
	return;
end

% Handle typical case.
%%% PLACEHOLDER
%xNew = xVals(end) - fVals(end)*(xVals(end)-xVals(end-1))/(fVals(end)-fVals(end-1));
%return;
xNew = findGoodCand( xVals_sorted, fVals_sorted );
if ( isrealscalar(xNew) )
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'typical method'." );
	thisFile = "RETURNED FROM groot1d__getXNew";
	return;
end

% Handle "cylce: expand left, expand right, bisect"?
% Handle "all values are the same".
%%% HOW SHOULD THIS COMPETE WITH LARGER-SCALE TYPICAL CASE???
m = mod(desperationIter,5); desperationIter++;
if (0==m)
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - expand right'." );
	xNew = xMax + 0.3*(xMax-xMin);
	thisFile = "RETURNED FROM groot1d__getXNew";
	return;
elseif (1==m)
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - expand left'." );
	xNew = xMin - 0.3*(xMax-xMin);
	thisFile = "RETURNED FROM groot1d__getXNew";
	return;
end
%
msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - bisect largest gap'." );
[ largestGapSize, largestGapIndex ] = max(diff(xVals_sorted));
xNew = (xVals_sorted(largestGapIndex) + xVals_sorted(largestGapIndex+1))/2.0;
thisFile = "RETURNED FROM groot1d__getXNew";
return;
