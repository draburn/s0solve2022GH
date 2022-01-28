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
	if ( 0.0 <= fm*fp )
		% Only one real root, or, maginally two roots;
		%  but we can fairly handle this as one root.
		x0 = mygetfield( prm, "x0", mypow(-c0/c3,1.0/3.0) );
		xVals = fzero( funchF, x0 );
		return;
	end
	%
	% One of the roots should be to the left of xm, one between xm and xp, one to the right of xp.
	dx = xp-xm;
	iterMax = mygetfield( prm, "iterMax", 30 );
	iterNum = 1;
	while (1)
		xlo = xm-((2.0^iterNum)*dx);
		flo = funchF(xlo);
		if ( flo*fm < 0.0 )
			break;
		end
		iterNum++;
		if ( iterNum > iterMax )
			msg( __FILE__, __LINE__, "Could not find bounding point." );
			error( "Could not find bounding point." );
		end
	end
	while (1)
		xhi = xp+((2.0^iterNum)*dx);
		fhi = funchF(xhi);
		if ( fhi*fp < 0.0 )
			break;
		end
		iterNum++;
		if ( iterNum > iterMax )
			msg( __FILE__, __LINE__, "Could not find bounding point." );
			error( "Could not find bounding point." );
		end
	end
	xVals(1) = fzero( funchF, [ xlo, xm ] );
	xVals(2) = fzero( funchF, [ xm, xp ] );
	xVals(3) = fzero( funchF, [ xp, xhi ] );
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
