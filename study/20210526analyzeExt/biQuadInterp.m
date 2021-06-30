function [ xCand, meritCand, datOut ] = biQuadInterp( xVals, fVals, prm=[], datIn=[] )
	commondefs;
	thisFile = "biQuadInterp";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	datOut = [];
	%
	% Validate input...
	numPts = size(xVals,2);
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	assert( 3 <= numPts );
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert(xValsAreStrictlyIncreasing);
	fValsAreAllNonzero = (0==sum( 0.0 == fVals ));
	assert( fValsAreAllNonzero );
	signF = sign(fVals(1));
	gVals = signF * fVals;
	fValsAllHaveSameSign = (0==sum( 0.0 >= gVals ));
	assert( fValsAllHaveSameSign );
	%
	%
	% Identify point-wise minimum of g.
	[ gOfMin, nOfMin ] = min( gVals );
	% Allow exactly equal values of ptwise min.
	% Ignore this issue except where it would prohibit calculation:
	%  when the ptwise min is at the edge.
	if ( 1 == nOfMin )
	if ( gVals(2) == gVals(1) )
		nOfMin = 2;
	end
	end
	if ( numPts == nOfMin )
	if ( gVals(numPts-1) == gVals(numPts) )
		% I believe min() will never do this, but let's be safe.
		nOfMin = numPts-1;
	end
	end
	ptwiseAbsMinIsNotOnEdge = ( 2 <= nOfMin ) && ( numPts-1 >= nOfMin );
	assert( ptwiseAbsMinIsNotOnEdge );
	%
	maxRatioAllowed = mygetfield( prm, "maxRatioAllowed", 10.0 );
	ratioTarget = mygetfield( prm, "ratioTarget", 5.0 );
	assert( isrealscalar(maxRatioAllowed) );
	assert( 1.0 < maxRatioAllowed );
	assert( isrealscalar(ratioTarget) );
	assert( 1.0 < ratioTarget );
	assert( ratioTarget < maxRatioAllowed );
	if ( xVals(nOfMin-1) < xVals(nOfMin) - maxRatioAllowed*(xVals(nOfMin+1)-xVals(nOfMin)) )
		xCand = xVals(nOfMin) - ratioTarget*(xVals(nOfMin+1)-xVals(nOfMin));
		meritCand = -1.0;
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "Left-balancing ( %12.8f, [%12.8f], %12.8f, %12.8f ).", ...
		  xVals(nOfMin-1), xCand, xVals(nOfMin), xVals(nOfMin+1) ) );
		return;
	end
	if ( xVals(nOfMin+1) > xVals(nOfMin) + maxRatioAllowed*(xVals(nOfMin)-xVals(nOfMin-1)) )
		xCand = xVals(nOfMin) + ratioTarget*(xVals(nOfMin)-xVals(nOfMin-1));
		meritCand = -1.0;
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "Right-balancing ( %12.8f, %12.8f, [%12.8f], %12.8f ).", ...
		  xVals(nOfMin-1), xVals(nOfMin), xCand, xVals(nOfMin+1) ) );
		return;
	end
	%
	x0 = xVals(nOfMin);
	x1 = xVals(nOfMin+1)-xVals(nOfMin);
	yVals = (xVals-x0)/x1;
	%
	% Look at quadratic model about nOfMin, pt "A".
	vecYA = yVals(nOfMin-1:nOfMin+1)';
	vecGA = gVals(nOfMin-1:nOfMin+1)';
	matYA = [ ones(3,1), vecYA, vecYA.^2 ];
	vecCA = matYA \ vecGA;
	curvatureAIsPositive = (0.0<vecCA(3)); 
	assert( curvatureAIsPositive );
	yExtA = -vecCA(2)./(2.0*vecCA(3));
	extAIsInBounds = ( (yVals(nOfMin-1)<=yExtA) && (yExtA<=yVals(nOfMin+1)) );
	assert( extAIsInBounds );
	%
	% For now...
	xCand = x0 + (x1*yExtA);
	meritCand = -1.0;
	%
return;
end
