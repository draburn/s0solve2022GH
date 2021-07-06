	commondefs;
	thisFile = "extFinder__init";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	datOut = [];
	%
	haveCand = false;
	%
	numPts = size(xVals,2);
	assert( numPts >= 5 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	%
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	fValsAreAllNonzero = (0==sum( 0.0 == fVals ));
	assert( fValsAreAllNonzero );
	%
	fValsAllHaveSameSign = (0==sum( 0.0 >= sign(fVals(1)) * fVals ));
	assert( fValsAllHaveSameSign );
	gVals = abs(fVals);
	%
	%
	%
	% Find nOfPtWiseMin and validate g'.
	for n=1:numPts
		if ( n==numPts )
		error( "Bad absFVals: No internal point-wise min; last point is actual min." );
		end
		if ( gVals(n+1) >= gVals(n) )
			break;
		end
	end
	nOfPtWiseMin = n;
	if ( fleq( gVals(nOfPtWiseMin), gVals(nOfPtWiseMin+1), eps ) )
		% Unlikely in real-world scenarios, so, efficient handling is unimportant.
		for n=nOfPtWiseMin+1:numPts-1
		if ( gVals(n+1) <= gVals(n) )
			error( "Bad absFVals: Point-wise min is not unique and/or there is a ptwise local max." );
		end
		end
		xOfCand = ( xVals(nOfPtWiseMin+1) + xVals(nOfPtWiseMin) ) / 2.0;
		meritOfCand = -1.0;
		haveCand = true;
		thisFile = [ "RETURN from " thisFile ];
		return;
	end
	if ( 1 == nOfPtWiseMin )
		error( "Bad absFVals: No internal point-wise min; first point is actual min." );
	end
	for n=nOfPtWiseMin:numPts-1
	if ( gVals(n+1) <= gVals(n) )
		error( "Bad absFVals: Point-wise min is not unique and/or there is a ptwise local max." );
	end
	end
	%
	%
	% Validate g''.
	if ( nOfPtWiseMin >= 4 )
		n = 2;
		p = polyfit( xVals(n-1:n+1), gVals(n-1:n+1), 2 );
		msg( thisFile, __LINE__, sprintf( ...
		  "    [ %3d,   %g,   %g ]", n, xVals(n), p(2) ) );
		curvatureOnLeftIsOkay = abs(p(1)) > sqrt(eps)*abs(p(2));
		assert( curvatureOnLeftIsOkay );
		pToCompare = p;
		for n=3:nOfPtWiseMin-2
			p = polyfit( xVals(n-1:n+1), gVals(n-1:n+1), 2 );
			msg( thisFile, __LINE__, sprintf( ...
			  "    [ %3d,   %g,   %g ]", n, xVals(n), p(2) ) );
			curvatureOnLeftIsOkay = abs(p(1)) > sqrt(eps)*abs(p(2));
			assert( curvatureOnLeftIsOkay );
			curvatureOnLeftIsOkay = pToCompare(1)*p(1) > 0.0;
			assert( curvatureOnLeftIsOkay );
		end
	end
	if ( nOfPtWiseMin <= numPts-3 )
		n = numPts-1;
		p = polyfit( xVals(n-1:n+1), gVals(n-1:n+1), 2 );
		msg( thisFile, __LINE__, sprintf( ...
		  "    [ %3d,   %g,   %g ]", n, xVals(n), p(2) ) );
		curvatureOnRightIsOkay = abs(p(1)) > sqrt(eps)*abs(p(2));
		assert( curvatureOnRightIsOkay );
		curvatureRight = p(1);
		for n=nOfPtWiseMin+2:numPts-2
			p = polyfit( xVals(n-1:n+1), gVals(n-1:n+1), 2 );
			msg( thisFile, __LINE__, sprintf( ...
			  "    [ %3d,   %g,   %g ]", n, xVals(n), p(2) ) );
			curvatureOnRightIsOkay = abs(p(1)) > sqrt(eps)*abs(p(2));
			assert( curvatureOnRightIsOkay );
			curvatureOnRightIsOkay = pToCompare(1)*p(1) > 0.0;
			assert( curvatureOnRightIsOkay );
		end
	end
	%
	%
	xOfPtWiseMin = xVals(nOfPtWiseMin);
	fOfPtWiseMin = fVals(nOfPtWiseMin);
	gOfPtWiseMin = gVals(nOfPtWiseMin);
	if ( 3 <= nOfPtWiseMin )
		haveLeft = true;
	else
		haveLeft = false;
	end
	if ( numPts-2 >= nOfPtWiseMin )
		haveRight = true;
	else
		haveRight = false;
	end
	assert( haveLeft || haveRight );
	%
	%
	n = nOfPtWiseMin;
	bigDelta = xVals(n+1) - xVals(n-1);
	bigG0 = gVals(n);
	bigG1 = gVals(n+1) + gVals(n-1) - 2.0*gVals(n);
	%
	%
	thisFile = [ "RETURN from " thisFile ];
	return;
