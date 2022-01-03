function x = getQuadGoodPt( a, b, c, tol=sqrt(eps) )
	warning( "This function is deprecated; use calcLinishRootOfQuad() instead." );
	thisFile = "getQuadGoodPt";
	% y = (a*(x^2)) + (b*x) + c;
	assert( isrealscalar(a) );
	assert( isrealscalar(b) );
	assert( isrealscalar(c) );
	if ( 0.0 == b )
		if ( a*c >= 0.0 )
			% We're at the extremum and there's no root.
			x = 0.0;
			return;
		end
		% Return the + root.
		msg( thisFile, __LINE__, "Warning: Quadratic has no 'good' point." );
		x = sqrt(-c/a);
		return;
	end
	discrim = (b^2) - (4.0*a*c);
	if ( discrim < 0.0 )
		% No zero; return extermum.
		x = -b/(2.0*a);
		return;
	end
	if ( abs(a*c) < tol*(b^2) )
		% Use near-linear model.
		x = -c*( 1.0 - ((a*c)/(b^2)) )/b;
		return;
	end
	% Use general quadratic model.
	x = b * ( sqrt( 1.0 - ((4.0*a*c)/(b^2)) ) - 1.0 ) / (2.0*a);
return

%!test
%!	setprngstates();
%!	a = randn()
%!	b = randn()
%!	c = randn()
%!	funchY = @(x)( (a*(x.^2)) + (b*x) + c );
%!	xGP = getQuadGoodPt( a, b, c )
%!	yGP = funchY(xGP)
%!	xVals = linspace( -1.5*abs(xGP)-1.0, 1.5*abs(xGP)+1.0, 1000 );
%!	yVals = funchY(xVals);
%!	y0 = funchY(0.0)
%!	zVals = b*xVals+c;
%!	plot( ...
%!	  xVals, yVals, 'o-', ...
%!	  xVals, zVals, '-', ...
%!	  0.0, y0, 's', 'markersize', 25, 'linewidth', 4, ...
%!	  xGP, yGP, '*', 'markersize', 25, 'linewidth', 4, ...
%!	  xVals, 0.0*yVals, 'k-', 'linewidth', 2 );
%!	grid on;
%!	isAZero = ( abs(yGP) <= sqrt(eps) )
%!	isAnExt = ( abs(b+(2.0*a*xGP)) <= sqrt(eps) )
%!	assert( abs(yGP) <= abs(y0) );
%!	assert( isAZero || isAnExt )
