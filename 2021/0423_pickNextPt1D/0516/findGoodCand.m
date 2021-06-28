function xCand = findGoodCand( xVals, fVals, prm = [] )
  	% Should-be-precompiled...
	commondefs;
	thisFile = "findGoodCand.m";
	%
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
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
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "Bounded root between x = %f and %f.",xVals(n),xVals(n+1) ) );
		%
		% Consider two quadratic models, look at "merit" of each.
		% Negative merit would mean "do not use".
		% Only run the model if the point exists.
		% Perhaps require monotonicity in g?
		% Update 2021.05.19:
		%  Realized we can have an "other-side" extremum.
		%  So, we will not use quad in that scenario.
		%  In fact, let's require monotonicity in g.
		meritLeft = -1.0;
		meritRight = -1.0;
		prm_boundedQuad = mygetfield( prm, "prm_boundedQuad", [] );
		if ( n >= 2 )
		if ( gVals(n-1) > gVals(n) )
			[ xLeft, meritLeft ] = findGoodCand__boundedQuad( ...
			  xVals(n-1), xVals(n), xVals(n+1), ...
			  gVals(n-1), gVals(n), gVals(n+1), ...
			  prm_boundedQuad );
			assert( isrealscalar(meritLeft) );
			if ( meritLeft > 0.0 )
				assert( isrealscalar(xLeft) );
				assert( xVals(n) < xLeft );
				assert( xLeft < xVals(n+1) );
			end
		end
		end
		if ( n+2 <= numPts )
		if ( gVals(n+1) > gVals(n+2) )
			[ xRight, meritRight ] = findGoodCand__boundedQuad( ...
			  xVals(n), xVals(n+1), xVals(n+2), ...
			  gVals(n), gVals(n+1), gVals(n+2), ...
			  prm_boundedQuad );
			assert( isrealscalar(meritRight) );
			if ( meritRight > 0.0 )
				assert( isrealscalar(xRight) );
				assert( xVals(n) < xRight );
				assert( xRight < xVals(n+1) );
			end
		end
		end
		%
		if ( meritLeft < 0.0 )
			if ( meritRight < 0.0 )
				msg_copious( verbLev, thisFile, __LINE__, ...
				  "xNew via 'bounded - linear interpolation'." );
				xCand = xVals(n) - gVals(n) ...
				 * ( xVals(n)-xVals(n+1) ) ...
				 / ( gVals(n)-gVals(n+1) );
			else
				msg_copious( verbLev, thisFile, __LINE__, ...
				  "xNew via 'bounded - right-side quadratic interpolation'." );
				xCand = xRight;
			end
		else
			if ( meritRight < 0.0 )
				msg_copious( verbLev, thisFile, __LINE__, ...
				  "xNew via 'bounded - left-side quadratic interpolation'." );
				xCand = xLeft;
			else
				msg_copious( verbLev, thisFile, __LINE__, ...
				  "xNew via 'bounded - blended quadratic interpolation'." );
				xCand = ( meritLeft*xLeft + meritRight*xRight ) ...
				 / ( meritLeft + meritRight );
			end
		end
		return;
	end
	assert( fValsAllHaveSameSign );
	%
	if ( mygetfield(prm,"usefofprime",false) )
	% 3-pt reciprocal logarithmic derivative interp aka "fofprime".
	% This is good for glancing roots of order > 2.
	% But, we want to be strict about when we use it.
	% For simplicity, consider only one-sided models,
	% BUT, if purely one-sided interp, watch out for
	%  hitting the same point repeatedly!
	prm_fofprime = mygetfield( prm, "prm_fofprime", [] );
	[ xTemp, meritTemp ] = findGoodCand__fofprime( xVals, gVals, prm_fofprime );
	if ( 0.0 < meritTemp )
		msg_copious( verbLev, thisFile, __LINE__, ...
		  "  xNew via 'f over f prime'." );
		xCand = xTemp;
		return;
	end
	end
	%
	%
	% 3-pt quad interp that has a root.
	% We'll take the extremum.
	if ( 3 <= numPts )
	for n=2:numPts-1
	if ( gVals(n-1) >= gVals(n) && gVals(n+1) >= gVals(n) )
	if ( gVals(n-1) != gVals(n+1) )
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "Looking at pt-wise ext quadratic interpolation..." ) );
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "{ %f, %f, %f }, { %f, %f, %f }.", ...
		  xVals(n-1), xVals(n), xVals(n+1), ...
		  gVals(n-1), gVals(n), gVals(n+1) )  );
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
				msg_copious( verbLev, thisFile, __LINE__, ...
				  "  xNew via 'pt-wise extremum - quadratic interpolation'." );
				xCand = xTemp;
				return;
			end
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "  REJECTED this 'pt-wise extremum - quadratic interpolation'." );
		end
	end
	end
	end
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
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "xNew via 'left-side linear interpolation'." );
			return;
		end
		xTemp = ( x3*g2 - x2*g3 ) / ( g2 - g3 );
		assert( xTemp < x2 );
		if ( xTemp > x1 )
			xCand = max([ xTemp, (x2+x1)/2.0 ]);
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "xNew via 'right-side linear interpolation'." );
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
		msg_copious( verbLev, thisFile, __LINE__, ...
		  "xNew via 'rightward linear extrapolation'." );
		return;
	end
	if ( gVals(end) < gVals(end-1) )
		xCand = ( xVals(end)*gVals(end-1) - xVals(end-1)*gVals(end) ) / (gVals(end-1)-gVals(end));
		msg_copious( verbLev, thisFile, __LINE__, ...
		  "xNew via 'leftward linear extrapolation'." );
		return;
	end
	%
	% Explore local pt-wise min, make sure we've found actual min.
	%error( "Not implemented!" );
	%
	msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
	  "Using end-scale lin extrap with points ( %g, %g ) and ( %g, %g ).", ...
	  xVals(1), fVals(1), xVals(end), fVals(end) ) );
	xCand = ( xVals(end)*fVals(1) - xVals(1)*fVals(end) ) / (fVals(1)-fVals(end));
	msg_copious( verbLev, thisFile, __LINE__, ...
	  sprintf( "New point is ( %g, ??? ).", xCand ) );
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
