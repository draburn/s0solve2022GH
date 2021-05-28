% Ideally pre-compiled stuff...
commondefs; thisFile = "analyzeExt__init";

% Verbosity...
verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );

% Check data types...
numPts = size(xVals,2);
assert( isrealarray(xVals,[1,numPts]) );
assert( isrealarray(fVals,[1,numPts]) );

% Check unsupported cases...
assert( numPts >= 3 );
%
xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
assert(xValsAreStrictlyIncreasing);
%
fValsAreAllNonzero = (0==sum( 0.0 == fVals ));
assert( fValsAreAllNonzero );
%
fValsAllHaveSameSign = (0==sum( 0.0 >= sign(fVals(1)) * fVals ));
assert( fValsAllHaveSameSign );
%
gVals = abs(fVals);
[ foo, indexOfGMin ] = min( gVals );
clear foo;
if ( 1==indexOfGMin )
	if ( gVals(2)!=gVals(1) )
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	for n=2:numPts-1
	if ( gVals(n+1)<=gVals(n) )
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	end
	haveExactlyOneMinAndNoMax = true;
	haveUniquePtWiseMin = false; % This scenario is probably unlikely.
	indexOfGMin = 2;
elseif ( numPts==indexOfGMin )
	if ( gVals(numPts-1)!=gVals(numPts) )
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	% Assuming min() returns the first min, this path should be impossible.
	for n=1:numPts-2
	if ( gVals(n+1)>=gVals(n) )
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	end
	haveExactlyOneMinAndNoMax = true;
	haveUniquePtWiseMin = false; % This scenario is probably unlikely.
	indexOfGMin = numPts-1;
else
	for n=1:indexOfGMin-1
	if ( gVals(n+1)>=gVals(n) )
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	end
	if (  gVals(indexOfGMin) < gVals(indexOfGMin-1) ...
	  &&  gVals(indexOfGMin) < gVals(indexOfGMin+1)  )
		haveUniquePtWiseMin = true;
	elseif (  gVals(indexOfGMin) == gVals(indexOfGMin-1) ...
	     &&   gVals(indexOfGMin) <  gVals(indexOfGMin+1)  )
		haveUniquePtWiseMin = false; % This scenario is probably unlikely.
		% Let indexOfGMin be whichever.
	elseif (  gVals(indexOfGMin) <  gVals(indexOfGMin-1) ...
	     &&   gVals(indexOfGMin) == gVals(indexOfGMin+1)  )
		haveUniquePtWiseMin = false; % This scenario is probably unlikely.
		% Let indexOfGMin be whichever.
	else
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	for n=indexOfGMin+1:numPts-1
	if ( gVals(n+1)<=gVals(n) )
		haveExactlyOneMinAndNoMax = false;
		assert(haveExactlyOneMinAndNoMax);
	end
	end
end
% I intend no special handling for when haveUniquePtWiseMin is false,
% I just thought it was worth including.
%
diffGVals = diff(gVals);
