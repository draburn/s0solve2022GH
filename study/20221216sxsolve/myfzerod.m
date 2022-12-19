function [ xZero, datOut ] = myfzerod( funchF, xL, xR, prm=[] )
	xZero = [];
	datOut = [];
	datOut.fevalCount = 0;
	assert( isrealscalar(xL) );
	assert( isrealscalar(xR) );
	if (isempty(mygetfield(prm,"fL",[])))
		[ fL, dfdxL ] = funchF(xL);
		datOut.fevalCount++;
	else
		fL = prm.fL;
		dfdxL = prm.dfdxL;
	endif
	if (isempty(mygetfield(prm,"fR",[])))
		[ fR, dfdxR ] = funchF(xR);
		datOut.fevalCount++;
	else
		fR = prm.fR;
		dfdxR = prm.dfdxR;
	endif
	assert( isrealscalar(fL) );
	assert( isrealscalar(dfdxL) );
	assert( isrealscalar(fR) );
	assert( isrealscalar(dfdxR) );
	if ( fL * fR >= 0.0 )
		error( "Initial points do not bracket a zero." );
	endif
	%
	fTol = mygetfield( prm, "fTol", sqrt(eps)*sqrt( fL^2 + fR^2 ) );
	xTol = mygetfield( prm, "xTol", sqrt(eps)*sqrt( eps*(xL^2+xR^2) + (xR-xL)^2 ) );
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	%
	fAbsBestPrev = [];
	xAbsSpanPrev = [];
	for iterCount = 1 : iterLimit+1
		%[ xL, xR, fL, fR, dfdxL, dfdxR ]
		if ( min(abs([fR,fL])) < fTol )
			%msg( __FILE__, __LINE__, "SUCCESS: fAbsBest < fTol." );
			if ( abs(fL) < abs(fR) )
				xZero = xL;
			else
				xZero = xR;
			endif
			break;
		elseif ( abs(xR-xL) < xTol )
			%msg( __FILE__, __LINE__, "QUALIFIED SUCCESS: xAbsSpan < xTol." );
			if ( abs(fL) < abs(fR) )
				xZero = xL;
			else
				xZero = xR;
			endif
			break;
		elseif ( iterCount >= iterLimit )
			%msg( __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterLimit." );
			xZero = [];
			break;
		endif
		xCandL = __getCand_quad( xL, fL, dfdxL, xR, fR );
		xCandR = __getCand_quad( xR, fR, dfdxR, xL, fL );
		%
		haveHeadOnCollision = ( ~isempty(xCandL) && ~isempty(xCandR) && abs(xCandL-xCandR) < 0.01*abs(xR-xL) );
		haveStall = ( ~isempty(fAbsBestPrev) && ~isempty(xAbsSpanPrev) && min(abs([fR,fL])) > 0.5 * fAbsBestPrev && abs(xR-xL) > 0.6 * xAbsSpanPrev );
		if (haveHeadOnCollision)
			% We'll take whichever candidate is closer to the midpoint.
			if ( abs( xCandL - (xR+xL)/2.0 ) < abs( xCandR - (xR+xL)/2.0 ) )
				%msg( __FILE__, __LINE__, "Head-on L" );
				xNext = xCandL;
			else
				%msg( __FILE__, __LINE__, "Head-on R" );
				xNext = xCandR;
			endif
		elseif (haveStall)
			%msg( __FILE__, __LINE__, "Stall" );
			xNext = ( xR + xL )/2.0; % Bisection.
		elseif ( isempty(xCandL) && isempty(xCandR) )
			%msg( __FILE__, __LINE__, "Secant" );
			xNext = ( fR*xL - fL*xR ) / ( fR - fL ); % Secant.
		elseif ( isempty(xCandL) )
			%msg( __FILE__, __LINE__, "Only: R" );
			xNext = xCandR;
		elseif ( isempty(xCandR) )
			%msg( __FILE__, __LINE__, "Only: L" );
			xNext = xCandL;
		% Below here, both exist, so take whicheve has a smaller |f|.
		elseif ( abs(fL) < abs(fR) )
			%msg( __FILE__, __LINE__, "Better: L" );
			xNext = xCandL;
		else
			%msg( __FILE__, __LINE__, "Better: R" );
			xNext = xCandR;
		endif
		%
		[ fNext, dfdxNext ] = funchF(xNext);
		datOut.fevalCount++;
		%
		xAbsSpanPrev = abs( xR - xL );
		fAbsBestPrev = min(abs([ fR, fL ]));
		if ( fNext * fL > 0.0 )
			xL = xNext;
			fL = fNext;
			dfdxL = dfdxNext;
		else
			xR = xNext;
			fR = fNext;
			dfdxR = dfdxNext;
		endif
	endfor
	datOut.xL = xL;
	datOut.xR = xR;
	datOut.fL = fL;
	datOut.fR = fR;
	datOut.dfdxL = dfdxL;
	datOut.dfdxR = dfdxR;
return;
endfunction;


function xCand = __getCand_quad( x0, f0, dfdx0, x1, f1 )
	xSpan = x1 - x0; % Could be negative.
	if ( f0*dfdx0*xSpan > 0.0 )
		% Do not attempt a step based on a point where the slope is in the wrong direction.
		xCand = [];
		return;
	endif
	sRoots = roots([ f1 - ( f0 + dfdx0*xSpan ), dfdx0*xSpan, f0 ]);
	% Note: the descriminant must be positive, assuming (f0*f1) < 0.0.
	xCand = []; % For now.
	haveMultipleValidSolutions = false; % Unless...
	for s = reshape( sRoots, 1, [] )
	if ( isreal(s) && 0.0 <= s && s <= 1.0 )
		if (~isempty(xCand))
			haveMultipleValidSolutions = true;
		endif
		xCand = x0 + xSpan * s;
	endif
	endfor
	if ( ~isempty(xCand) && ~haveMultipleValidSolutions )
		return;
	endif
	msg( __FILE__, __LINE__, "Infodump..." );
	x0
	f0
	dfdx0
	x1
	f1
	xSpan
	sRoots
	if (haveMultipleValidSolutions)
		error( "SHOULD-BE-IMPOSSIBLE-ERROR: There was more than one root in the interval." );
	else
		error( "SHOULD-BE-IMPOSSIBLE-ERROR: There was no root in the interval." );
	endif
return;
endfunction

%!function [ f, dfdx ] = foo( x )
%!	f = x.^4 - 2.0;
%!	dfdx = 4.0*(x.^3);
%!endfunction

%!test
%!	clear;
%!	[ xZero, datOut ] = myfzerod( @(x)( foo(x) ), 0.0, 2.0 )
%!	abs(foo(xZero))
