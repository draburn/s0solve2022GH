function [ xCand, meritCand, datOut ] = quadSymInterp( xVals, fVals, prm=[], datIn=[] )
	commondefs;
	thisFile = "quadSymInterp";
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
	% Validate gVals...
	for n=1:numPts
		if ( n==numPts )
		error( "Bad absFVals: No internal ptwise min; last point is actual min." );
		end
		if ( gVals(n+1) >= gVals(n) )
			break;
		end
	end
	nOfMin = n;
	%
	if ( gVals(nOfMin+1) == gVals(nOfMin) )
		% Very unlike in real-world scenarios,
		% so, efficient handling is unimportant.
		for n=nOfMin+1:numPts-1
		if ( gVals(n+1) <= gVals(n) )
			error( "Bad absFVals: Ptwise min is not unique and/or there is a ptwise local max." );
		end
		end
		xCand = ( xVals(nOfMin+1) + xVals(nOfMin) ) / 2.0;
		meritCand = -1.0;
		return;
	end
	%
	if ( 1 == nOfMin )
		error( "Bad absFVals: No internal ptwise min; first point is actual min." );
	end
	%
	for n=nOfMin:numPts-1
	if ( gVals(n+1) <= gVals(n) )
		error( "Bad absFVals: Ptwise min is not unique and/or there is a ptwise local max." );
	end
	end
	%
	%
	% Do work...
	msg( thisFile, __LINE__, "THIS IS SUPER HACKISH." );
	n = nOfMin;
	if ( gVals(n-1) == gVals(n+1) )
		% Unlikely to happen in the wild.
		xCand = (xVals(n-1)+xVals(n+1))/2.0;
		meritCand = -1.0;
		return;
	elseif ( gVals(n-1) > gVals(n+1) )
		assert( n+3 <= numPts );
		assert( gVals(n-1) < gVals(n+3) ) % Not strictly necessary.
		nOfInterp = nOfMin+2;
		xOfTarget = xVals(n-1);
		gOfTarget = gVals(n-1);
		solutionSide = +1;
	else
		% gVals(n-1) < gVals(n+1)
		assert( n-3 >= 1 );
		assert( gVals(n-3) > gVals(n+1) ) % Not strictly necessary.
		nOfInterp = nOfMin-2;
		xOfTarget = xVals(n+1);
		gOfTarget = gVals(n+1);
		solutionSide = -1;
	end
	%	
	n = nOfInterp;
	vecX = xVals(n-1:n+1)';
	vecG = gVals(n-1:n+1)';
	matX = [ ones(3,1), vecX, vecX.^2 ];
	vecC = matX \ vecG;
	assert( 0.0 < vecC(3) );
	%xExtQuad = -vecC(2)/(2.0*vecC(3));
	discrim = vecC(2)^2 - 4.0 * vecC(3) * (vecC(1) - gOfTarget);
	assert( 0.0 < discrim );
	xSolu = ( -vecC(2) + solutionSide*sqrt(discrim) ) / (2.0*vecC(3) );
	xCand = (xSolu+xOfTarget)/2.0;
	meritCand = 0.0;
	return;
end
