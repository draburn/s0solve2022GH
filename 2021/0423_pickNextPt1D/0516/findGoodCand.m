function xCand = findGoodCand( xVals, fVals, prm = [] )
  	% Should-be-precompiled...
	thisFile = "findGoodCand.m";
	%
	% Check data types...
	numPts = size(xVals,2);
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	%
	% Check unsupported cases...
	assert( 2 <= numPts );
	%
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert(xValsAreStrictlyIncreasing);
	%
	fValsAreAllNonzero = (0==sum( 0.0 == fVals ));
	assert( fValsAreAllNonzero );
	%
	signF = sign(fVals(1));
	gVals = signF * fVals;
	fValsAllHaveSameSign = (0==sum( 0.0 >= gVals ));
	if (!fValsAllHaveSameSign)
		n = 1;
		while ( gVals(n+1) > 0.0 )
			n++;
		end
		x1 = xVals(n);
		x2 = xVals(n+1);
		g1 = gVals(n);
		g2 = gVals(n+1);
		assert( g1 > 0.0 );
		assert( g2 < 0.0 );
		xTemp = ( x1*g2 - x2*g1 ) / ( g2 - g1 );
		assert( xTemp > x1 );
		assert( xTemp < x2 );
		xCand = xTemp;
		return;
	end
	assert( fValsAllHaveSameSign );
	%
	%
	% Do work..
	if (1)
	%
	% 3-pt quad interp that has a root.
	% We'll take the extremum.
	if ( 3 <= numPts )
	for n=2:numPts-1
	if ( gVals(n-1) >= gVals(n) && gVals(n+1) >= gVals(n) )
	if ( gVals(n-1) != gVals(n+1) )
		matX = [ ones(3,1), xVals(n-1:n+1)', xVals(n-1:n+1).^2' ];
		vecG = gVals(n-1:n+1)';
		vecC = matX\vecG;
		a = vecC(3);
		b = vecC(2);
		c = vecC(1);
		assert( 0.0 < a );
		if ( b^2 >= 4.0*a*c )
			xTemp = -b / (2.0*a);
			if ( xVals(n-1) < xTemp && xTemp < xVals(n+1) )
				xCand = xTemp;
			end
				return;
		end
	end
	end
	end
	end
	%
	end
	%
	% 2-pt lin extrap... but internal.
	if ( 3<= numPts )
	for n=2:numPts-1
	if ( gVals(n) <= gVals(n-1) && gVals(n) <= gVals(n+1) )
	if ( gVals(n-1) != gVals(n+1) )
		x1 = xVals(n-1);
		x2 = xVals(n);
		x3 = xVals(n+1);
		g1 = gVals(n-1);
		g2 = gVals(n);
		g3 = gVals(n+1);
		xTemp = ( x1*g2 - x2*g1 ) / ( g2 - g1 );
		assert( xTemp > x2 );
		if ( xTemp < x3 )
			xCand = min([ xTemp, (x2+x3)/2.0 ]);
			return;
		end
		xTemp = ( x3*g2 - x2*g3 ) / ( g2 - g3 );
		assert( xTemp < x2 );
		if ( xTemp > x1 )
			xCand = max([ xTemp, (x2+x1)/2.0 ]);
			return;
		end
	end
	end
	end
	end
	%
	% 2-pt lin extrap.
	% Hypothetically, we could reject on basis of an apparent
	% horizontal asymptote above zero.
	if ( gVals(1) < gVals(2) )
		xCand = ( xVals(1)*gVals(2) - xVals(2)*gVals(1) ) / (gVals(2)-gVals(1));
		return;
	end
	if ( gVals(end) < gVals(end-1) )
		xCand = ( xVals(end)*gVals(end-1) - xVals(end-1)*gVals(end) ) / (gVals(end-1)-gVals(end));
		return;
	end
	%
	% Explore local pt-wise min, make sure we've found actual min.
	%error( "Not implemented!" );
	%
	msg( thisFile, __LINE__, sprintf( ...
	  "Using end-scale lin extrap with points ( %g, %g ) and ( %g, %g ).", ...
	  xVals(1), fVals(1), xVals(end), fVals(end) ) );
	xCand = ( xVals(end)*fVals(1) - xVals(1)*fVals(end) ) / (fVals(1)-fVals(end));
	msg( thisFile, __LINE__, sprintf( "New point is ( %g, ??? ).", xCand ) );
	return;
	%
	xCand = [];
	return;
	%
	%
	%
	%
	%
return;
end
