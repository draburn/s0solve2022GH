thisFile = "groot1d__getXNew";

if ( 1 == fevalCount )
	xNew = x2;
	return;
end

assert( 2 <= fevalCount );


% Handle "bounded".
%%% TO-DO

% Handle typical case.
%%% PLACEHOLDER
xNew = xVals(end) - fVals(end)*(xVals(end)-xVals(end-1))/(fVals(end)-fVals(end-1));
return;

xMin = min(xVals);
xMax = max(xVals);
% Handle "cylce: expand left, expand right, bisect"?
% Handle "all values are the same".
m = mod(fevalCount,5);
if (0==m)
	xNew = 2.0*xMax - xMin;
	return;
elseif (1==m)
	xNew = 2.0*xMin - xMax;
	return;
end
%
[ largestGapSize, largestGapIndex ] = max(diff(xVals_sorted));
xNew = (xVals_sorted(largestGapIndex) + xVals_sorted(largestGapIndex+1))/2.0;
return;
