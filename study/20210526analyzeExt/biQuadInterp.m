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
	
	% IN CONTRAST WITH lookatLoggishRes2.m,
	% BI-QUAD DOES NOT HELP?!?!
	
	% Everything needs rewriting.
	% Also needs testing.
	
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
	% Add full inspection of |f| values,
	% ensuring exactly one local minimum (or an adjacent tie)?
	%
	maxRatioAllowed = mygetfield( prm, "maxRatioAllowed", 10.0 );
	ratioTarget = mygetfield( prm, "ratioTarget", 5.0 );
	assert( isrealscalar(maxRatioAllowed) );
	assert( 1.0 < maxRatioAllowed );
	assert( isrealscalar(ratioTarget) );
	assert( 1.0 <= ratioTarget );
	assert( ratioTarget < maxRatioAllowed );
	if (0)
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
	end
	%
	% Scale values to improve accuracy...?
	x0 = xVals(nOfMin);
	x1 = xVals(nOfMin+1)-xVals(nOfMin);
	yVals = (xVals-x0)/x1;
	%
	% Look at quadratic model about nOfMin, pt "A".
	% As long as the abs value of the ptwise abs min is less than or equal to
	%  its neighbors (and strictly less than at least one),
	%  all of the assertions should be satisfied.
	vecY = yVals(nOfMin-1:nOfMin+1)';
	vecG = gVals(nOfMin-1:nOfMin+1)';
	matY = [ ones(3,1), vecY, vecY.^2 ];
	vecC = matY \ vecG;
	curvatureIsPositive = (0.0<vecC(3));
	assert( curvatureIsPositive );
	yExt = -vecC(2)./(2.0*vecC(3));
	extIsInBounds = ( (yVals(nOfMin-1)<=yExt) && (yExt<=yVals(nOfMin+1)) );
	assert( extIsInBounds );
	%
	xCand = x0 + (x1*yExt);
	meritCand = -1.0;
	% May override below.
	
	%msg_copious( verbLev, thisFile, __LINE__, ...
	%  "Toggle this \"return\", balancing to see merit of bi-quad." );
	%return;
	
	%
	% If we have only 3 points and hit the ptwise ext exactly,
	%  we should instead eval nearby, like ~1/5 way through the interval.
	% But, if we have 4+ points in this case, we can check for cnvg.
	% Worry about that stuff later.
	%
	deltaXTol = mygetfield( prm, "deltaXTol", sqrt(eps) );
	if ( abs(yExt) < deltaXTol/x1 )
		msg_warn( verbLev, thisFile, __LINE__, "Hit extMin. Not sure what to do here!!" );
		return;
	elseif ( yExt < 0.0 )
		yIntervalHi = yVals(nOfMin);
		yIntervalLo = yVals(nOfMin-1);
		if ( nOfMin >= 3 )
			nOfPtB = nOfMin-1;
		else
			msg_copious( verbLev, thisFile, __LINE__, "Not enough points for model to left." );
			return;
		end
	else
		yIntervalLo = yVals(nOfMin);
		yIntervalHi = yVals(nOfMin+1);
		if ( nOfMin <= numPts-2 )
			nOfPtB = nOfMin+1;
		else
			msg_copious( verbLev, thisFile, __LINE__, "Not enough points for model to right." );
			return;
		end
	end
	%
	% We should check balance.
	if ( xVals(nOfPtB-1) < xVals(nOfPtB) - maxRatioAllowed*(xVals(nOfPtB+1)-xVals(nOfPtB)) )
		msg_copious( verbLev, thisFile, __LINE__, "Adjacent model is unbalanced on left." );
		return;
	end
	if ( xVals(nOfPtB+1) > xVals(nOfPtB) + maxRatioAllowed*(xVals(nOfPtB)-xVals(nOfPtB-1)) )
		msg_copious( verbLev, thisFile, __LINE__, "Adjacent model is unbalanced on right." );
		return;
	end
	%
	vecYB = yVals(nOfPtB-1:nOfPtB+1)';
	vecGB = gVals(nOfPtB-1:nOfPtB+1)';
	matYB = [ ones(3,1), vecYB, vecYB.^2 ];
	vecCB = matYB \ vecGB;
	if ( 0.0 >= vecCB(3) )
		msg_copious( verbLev, thisFile, __LINE__, "Adjacent model has bad curvature." );
		return;
	end
	yExtB = -vecCB(2)./(2.0*vecCB(3))
	
		%yExtB = cap( yExtB, yIntervalLo, yIntervalHi );
		%%yExtB = cap( yExtB, yVals(nOfMin-1), yVals(nOfMin+1) );
		%yExt = (yExt+yExtB)/2.0;
		%xCand = x0 + (x1*yExt);
		%return;
	
	
	
	
	if (0)
	%yExtB = cap( yExtB, yVals(nOfMin-1), yVals(nOfMin+1) ); % This doesn't work.
	if ( yExtB < yIntervalLo )
		msg_copious( verbLev, thisFile, __LINE__, "Capping adj to yIntervalLo." );
		yExtB = yIntervalLo;
	elseif ( yExtB > yIntervalHi )
		msg_copious( verbLev, thisFile, __LINE__, "Capping adj to yIntervalHi." );
		yExtB = yIntervalHi;
	end
	%yExtB = cap( yExtB, yIntervalLo, yIntervalHi );
	end
	if ( yExtB*(yExtB-yVals(nOfPtB)) < 0.0 )
	%if (1)
		% yExtB is between YVals(nOfPtB) and yVals(nOfMin).
		yExt = (yExt+yExtB)/2.0;
		xCand = x0 + (x1*yExt);
		meritCand = -1.0;
	else
		msg_copious( verbLev, thisFile, __LINE__, "Adjacent model is out of bounds." );
		return;
	end
	%
return;
end
