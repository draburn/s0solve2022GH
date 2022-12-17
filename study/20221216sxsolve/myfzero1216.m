function [ x, datOut ] = myfzero1216( funchFC, x0, prm=[] )
	msg( __FILE__, __LINE__, "WARNING: THIS CODE HAS BEEN ABANDONED. (However, it may still work.)" );
	x = [];
	datOut = [];
	%
	xL = x0(1);
	xR = x0(2);
	if ( xL >= xR )
		error( "Initial points are not properly ordered." );
	endif
	if ( isempty(mygetfield( prm, "fL", []) ) )
		[ fL, cL, fpL, cpL ] = funchFC(xL);
	else
		fL = prm.fL;
		cL = prm.cL;
		fpL = prm.fpL;
		cpL = prm.cpL;
	endif
	if ( isempty(mygetfield( prm, "fR", []) ) )
		[ fR, cR, fpR, cpR ] = funchFC(xR);
	else
		fR = prm.fR;
		cR = prm.cR;
		fpR = prm.fpR;
		cpR = prm.cpR;
	endif
	if ( fL * fR > 0.0 )
		error( "Initial points do not bracket a zero." );
	endif
	charSide = '';
	if ( cL>=0.0 && cR>=0.0 )
		if ( abs(fL) < abs(fR) )
			charSide = 'L';
		else
			charSide = 'R';
		endif
	elseif ( cL>=0.0 )
		charSide = 'L';
	elseif ( cR>=0.0 )
		charSide = 'R';
	else
		error( "Initial points both fail constraint." );
	endif
	if ( 'L' == charSide )
		xCBest = xL;
		fCBest = fL;
		cCBest = cL;
		fpCBest = fpL;
		cpCBest = cpL;
	elseif ( 'R' == charSide )
		xCBest = xR;
		fCBest = fR;
		cCBest = cR;
		fpCBest = fpR;
		cpCBest = cpR;
	endif
	fTol = mygetfield( prm, "fTol", sqrt(eps)*sqrt(fL^2+fR^2) );
	xTol = mygetfield( prm, "xTol", sqrt(eps)*sqrt( eps*(xL^2+xR^2) + (xR-xL)^2 ) );
	cTol = mygetfield( prm, "cTol", sqrt(eps) );
	stepLimitCoeff = mygetfield( prm, "stepLimitCoeff", 0.7 );
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	%
	% First, we find the zero of f without worrying about the constraint.
	for n = 1 : iterLimit
		%[ xL, fL, cL, xR, fR, cR]
		if ( abs(xR-xL) <= xTol )
			break;
		elseif ( abs(fL) <= fTol || abs(fR) <= fTol )
			break;
		endif
		%
		dL = []; % Unless...
		if ( fL * fpL < 0.0 )
			dL = - fL / fpL; % But, may be capped.
			if ( abs(dL) > stepLimitCoeff*abs(xR-xL) )
				dL = stepLimitCoeff*abs(xR-xL);
			endif
		endif
		% The cap is based on the idea that, if the solution actually is closer to the other side,
		% then that should be reflected in the step from the other side.
		%
		dR = []; % Unless...
		if ( fR * fpR > 0.0 )
			dR = - fR / fpR; % But, may be capped.
			if ( abs(dR) > stepLimitCoeff*abs(xR-xL) )
				dR = -stepLimitCoeff*abs(xR-xL);
			endif
		endif
		%
		if ( isempty(dL) && isempty(dR) )
			xC = (xR+xL)/2.0;
		elseif ( isempty(dL) )
			xC = xR + dR;
		elseif ( isempty(dR) )
			xC = xL + dL;
		elseif ( abs(dL) < abs(dR) )
			xC = xL + dL;
		else
			xC = xR + dR;
		endif
		%
		[ fC, cC, fpC, cpC ] = funchFC(xC);
		if ( cC >= 0.0 && abs(fC) < abs(fCBest) )
			xCBest = xC;
			fCBest = fC;
			cCBest = cC;
			fpCBest = fpC;
			cpCBest = cpC;
		endif
		%
		charSide = '';
		if ( fC*fL >= 0.0 )
			charSide = 'L';
		elseif ( fC*fR >= 0.0 )
			charSide = 'R';
		else
			error( "We appear to have lost bracketting of the zero. This should be impossible." );
		endif
		if ( 'L' == charSide )
			xL = xC;
			fL = fC;
			cL = cC;
			fpL = fpC;
			cpL = cpC;
		elseif ( 'R' == charSide )
			xR = xC;
			fR = fC;
			cR = cC;
			fpR = fpC;
			cpR = cpC;
		endif
	endfor
	if ( abs(fL) < abs(fR) )
		xFBest = xL;
		fFBest = fL;
		cFBest = cL;
		fpFBest = fpL;
		cpFBest = cpL;
	else
		xFBest = xR;
		fFBest = fR;
		cFBest = cR;
		fpFBest = fpR;
		cpFBest = cpR;
	endif
	%
	if ( cFBest >= 0.0 )
		x = xFBest;
		return;
	endif
	%
	% Repeat process for "f" with "c".
	if ( xCBest < xFBest )
		xL = xCBest;
		fL = fCBest;
		cL = cCBest;
		fpL = fpCBest;
		cpL = cpCBest;
		xR = xFBest;
		fR = fFBest;
		cR = cFBest;
		fpR = fpFBest;
		cpR = cpFBest;
	else
		xR = xCBest;
		fR = fCBest;
		cR = cCBest;
		fpR = fpCBest;
		cpR = cpCBest;
		xL = xFBest;
		fL = fFBest;
		cL = cFBest;
		fpL = fpFBest;
		cpL = cpFBest;
	endif
	for n = 1 : iterLimit
		%[ xL, fL, cL, xR, fR, cR]
		if ( abs(xR-xL) <= xTol )
			break;
		elseif ( abs(cL) <= cTol || abs(cR) <= cTol )
			break;
		endif
		%
		dL = []; % Unless...
		if ( cL * cpL < 0.0 )
			dL = - cL / cpL; % But, may be capped.
			if ( abs(dL) > stepLimitCoeff*abs(xR-xL) )
				dL = stepLimitCoeff*abs(xR-xL);
			endif
		endif
		% The cap is based on the idea that, if the solution actually is closer to the other side,
		% then that should be reflected in the step from the other side.
		%
		dR = []; % Unless...
		if ( cR * cpR > 0.0 )
			dR = - cR / cpR; % But, may be capped.
			if ( abs(dR) > stepLimitCoeff*abs(xR-xL) )
				dR = -stepLimitCoeff*abs(xR-xL);
			endif
		endif
		%
		if ( isempty(dL) && isempty(dR) )
			xC = (xR+xL)/2.0;
		elseif ( isempty(dL) )
			xC = xR + dR;
		elseif ( isempty(dR) )
			xC = xL + dL;
		elseif ( abs(dL) < abs(dR) )
			xC = xL + dL;
		else
			xC = xR + dR;
		endif
		%
		[ fC, cC, fpC, cpC ] = funchFC(xC);
		%
		charSide = '';
		if ( cC*cL >= 0.0 )
			charSide = 'L';
		elseif ( cC*cR >= 0.0 )
			charSide = 'R';
		else
			error( "We appear to have lost bracketting of the constraint. This should be impossible." );
		endif
		if ( 'L' == charSide )
			xL = xC;
			fL = fC;
			cL = cC;
			fpL = fpC;
			cpL = cpC;
		elseif ( 'R' == charSide )
			xR = xC;
			fR = fC;
			cR = cC;
			fpR = fpC;
			cpR = cpC;
		endif
	endfor
	if ( abs(cL) < abs(cR) )
		xCBest = xL;
		fCBest = fL;
		cCBest = cL;
		fpCBest = fpL;
		cpCBest = cpL;
	else
		xCBest = xR;
		fCBest = fR;
		cCBest = cR;
		fpCBest = fpR;
		cpCBest = cpR;
	endif
	x = xCBest;
	return;
return;
endfunction
