% function xMCRoot = getMCRootOfFit( xL, fL, dfdxL, xR, fR, dfdxR ).
% Gets root of monotonic cubic fit to data, if it exists.
% DRaburn 2022-12-18.

function xMCRoot = getMCRootOfFit( xL, fL, dfdxL, xR, fR, dfdxR, prm=[] )
	assert( isrealscalar(xL) );
	assert( isrealscalar(fL) );
	assert( isrealscalar(dfdxL) );
	assert( isrealscalar(xR) );
	assert( isrealscalar(fR) );
	assert( isrealscalar(dfdxR) );
	if ( xL >= xR )
		error( "Initial points are not well-ordered." );
	elseif ( fL*fR >= 0.0 )
		error( "Initial points do not bracket a root." );
	endif
	xSpan = xR - xL;
	foo = [ 3.0, -1.0; -2.0, 1.0 ] * [ fR - fL - dfdxL * xSpan; (dfdxR - dfdxL) * xSpan ];
	c0 = fL;
	c1 = dfdxL * xSpan;
	c2 = foo(1);
	c3 = foo(2);
	%
	sExts = roots([ 3.0*c3, 2.0*c2, c1 ]);
	for s = reshape( sExts, 1, [] )
	if ( isreal(s) && 0.0 < s && s < 1.0 )
		if ( mygetfield( prm, "verbLev", 0 ) > 0 )
			msg( __FILE__, __LINE__, "Cubic is not monotonic over interval." );
			xExts = sExts*xSpan
		endif
		xMCRoot = [];
		return;
	endif
	endfor
	%
	foundMultipleRoots = false;
	xMCRoot = [];
	sRoots = roots([ c3, c2, c1, c0 ]);
	for s = reshape( sRoots, 1, [] )
	if ( isreal(s) && 0.0 <= s && s <= 1.0 )
		if ( ~isempty(xMCRoot) )
			msg( __FILE__, __LINE__, "ERROR: Found multiple roots, despite checks so far. This should be impossible." );
			xMCRoot_prev = xMCRoot
			foundMultipleRoots = true;
		endif
		xMCRoot = xL + s * xSpan;
	endif
	endfor
	assert( ~foundMultipleRoots );
	assert( ~isempty(xMCRoot) );
return;
endfunction

%!test
%!	clear;
%!	setprngstates();
%!	%
%!	xL = randn()*exp(3.0*randn());
%!	xR = xL + abs(randn()*exp(3.0*randn()));
%!	fL = randn()*exp(3.0*randn());
%!	fR = -sign(fL)*abs(randn()*exp(3.0*randn()));
%!	dfdxL = randn()*exp(3.0*randn());
%!	dfdxR = sign(dfdxL)*abs(randn()*exp(3.0*randn()));
%!	%
%!	prm = [];
%!	prm.verbLev = 100;
%!	xMCRoot = getMCRootOfFit( xL, fL, dfdxL, xR, fR, dfdxR, prm );
%!	%
%!	sViz = linspace(0.0,1.0,101);
%!	xSpan = xR - xL;
%!	xViz = xL + xSpan*sViz;
%!	foo = [ 3.0, -1.0; -2.0, 1.0 ] * [ fR - fL - dfdxL * xSpan; (dfdxR - dfdxL) * xSpan ];
%!	c0 = fL;
%!	c1 = dfdxL * xSpan;
%!	c2 = foo(1);
%!	c3 = foo(2);
%!	fViz = c0 + c1*sViz + c2*(sViz.^2) + c3*(sViz.^3);
%!	fLViz = fL + dfdxL*xSpan*sViz;
%!	fRViz = fR + dfdxR*xSpan*(sViz-1.0);
%!	%
%!	plot( ...
%!	  xViz, fViz, 'o-', ...
%!	  xViz, 0.0*xViz, 'k-', ...
%!	  xL, fL, 'v', 'markerSize', 20, ...
%!	  xR, fR, '^', 'markerSize', 20, ...
%!	  xViz, fLViz, '-', ...
%!	  xViz, fRViz, '-' );
%!	grid on;
%!	if ( isempty(xMCRoot) )
%!		msg( __FILE__, __LINE__, "MCRoot does not exist." );
%!	else
%!		msg( __FILE__, __LINE__, sprintf( "MCRoot is at %f.", xMCRoot ) );
%!		hold on;
%!		plot( xMCRoot, 0.0, '*', 'markerSize', 20 );
%!		hold off;
%!	endif
%!	msg( __FILE__, __LINE__, "Please inspect the graph." );
