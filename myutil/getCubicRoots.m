function xVals = getCubicRoots( c0, c1, c2, c3, prm=[] )
	assert( c3 != 0.0 );
	%
	function f = funcF( dummyX, c0, c1, c2, c3 )
		f = c0 + (dummyX.*( c1 + (dummyX.*( c2 + (dummyX*c3) )) ));
	return;
	end
	funchF = @(dummyX) funcF( dummyX, c0, c1, c2, c3 );
	%
	discrim = ((2.0*c2)^2) - (12.0*c1*c3);
	if ( 0.0 >= discrim )
		% Cubic has only one real root.
		x0 = mygetfield( prm, "x0", mypow(-c0/c3,1.0/3.0) );
		xVals = fzero( funchF, x0 );
		return;	
	end
	%
	% Cubic may have up to three real roots.
	sqrtDiscrim = sqrt(discrim);
	if ( c3 > 0.0 )
		xm = ( (-2.0*c2) - sqrtDiscrim ) / (6.0*c3);
		xp = ( (-2.0*c2) + sqrtDiscrim ) / (6.0*c3);
	else
		xm = ( (-2.0*c2) + sqrtDiscrim ) / (6.0*c3);
		xp = ( (-2.0*c2) - sqrtDiscrim ) / (6.0*c3);
	end
	fm = funchF(xm);
	fp = funchF(xp);
	%
	if ( 0.0 < fm*fp )
		% Only one real root.
		if ( 0.0 < fm*c3 )
			x0 = mygetfield( prm, "x0", xm-(xp-xm) )
		else
			x0 = mygetfield( prm, "x0", xp+(xp-xm) )
		end
		xVals = fzero( funchF, x0 );
		return;
	end
	%
	% Consdier "marginal" cases, which have two roots.
	if ( 0.0 == fm )
		msg( __FILE__, __LINE__, "NOT IMPLEMENTED!" );
		error( "NOT IMPLEMENTED!" );
	elseif ( 0.0 == fp )
		msg( __FILE__, __LINE__, "NOT IMPLEMENTED!" );
		error( "NOT IMPLEMENTED!" );
	end
	%
	% One of the roots should be to the left of xm, one between xm and xp, one to the right of xp.
	msg( __FILE__, __LINE__, "THIS CODES IS A FAIR ATTEMPT, BUT OFTEN FAILS;" );
	msg( __FILE__, __LINE__, "TODO: RIGIROUSLY RE-IMPLEMENT UN-BOUND CASES." );
	x1 = xm-(eps^0.50)*(xp-xm)
	xr = [ xm, xp ]
	x3 = xp+(eps^0.50)*(xp-xm)
	xVals(1) = fzero( funchF, x1 );
	xVals(2) = fzero( funchF, [ xm, xp ] );
	xVals(3) = fzero( funchF, x3 )
	%
	if ( xVals(1) >= xm || xm >= xVals(2) || xVals(2) >= xp || xp >= xVals(3) )
		msg( __FILE__, __LINE__, "Calculation failed, but a more rigirous version might work." );
		error( "Calculation failed, but a more rigirous version might work." );
	end
return
end

%!test
%!	%setprngstates(72583808);
%!	%setprngstates(55260128);
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	a = randn();
%!	b = randn();
%!	c = randn();
%!	d = randn();
%!	funchF = @(x) a + (x.*( b + (x.*( c + (x*d) )) ));
%!	%
%!	xVals = getCubicRoots( a, b, c, d );
%!	numVals = size(xVals,2);
%!	assert( 3 >= numVals );
%!	assert( isrealarray(xVals,[1,numVals]) );
%!	fVals = funchF(xVals);
%!	%
%!	xLo = min(xVals) - 1.0 - 1.5*(max(xVals)-min(xVals));
%!	xHi = max(xVals) + 1.0 + 1.5*(max(xVals)-min(xVals));
%!	%xLo = min([-10.0, xLo ] );
%!	%xHi = max([-10.0, xHi ] );
%!	numPts = 1001;
%!	xPts = linspace( xLo, xHi, numPts );
%!	fPts = funchF(xPts);
%!	%
%!	numFigs++; figure(numFigs);
%!	plot( xPts, fPts, 'o-' );
%!	hold on;
%!	switch(numVals)
%!	case 1
%!		plot( ...
%!		  xVals, fVals, 'rs', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals, fVals, 'rx', 'markersize', 20, 'linewidth', 3 );
%!	case 2
%!		plot( ...
%!		  xVals(1), fVals(1), 'gs', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(1), fVals(1), 'gx', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(2), fVals(2), 'gs', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(2), fVals(2), 'gx', 'markersize', 20, 'linewidth', 3 );
%!	case 3
%!		plot( ...
%!		  xVals(1), fVals(1), 'rs', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(1), fVals(1), 'rx', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(2), fVals(2), 'cs', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(2), fVals(2), 'cx', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(3), fVals(3), 'rs', 'markersize', 20, 'linewidth', 3, ...
%!		  xVals(3), fVals(3), 'rx', 'markersize', 20, 'linewidth', 3 );
%!	otherwise
%!		error( "Invalid numVals." );
%!	end
%!	hold off;
%!	grid on;
%!	%
%!	if ( 2 == numVals )
%!		assert( xVals(1) <= xVals(2) );
%!		if ( 3 == numVals )
%!			assert( xVals(2) <= xVals(3) );
%!		end
%!	end
%!	%
%!	for n=1:numVals
%!		x = xVals(n);
%!		f = a + (x.*( b + (x.*( c + (x*d) )) ));
%!		epsF = (eps^0.50)*( abs(a) + abs(b*x) + abs(c*(x^2)) + abs(d*(x^3)) );
%!		assert( abs(f) < epsF );
%!	end
