% Function...

function [ vecXF, datOut ] = findZero_fullH( vecX0, funchF, prm=[] )
	% Look at basics.
	fevalCount = 0;
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	cfdjk_prm = mygetfield( prm, "cfdjk_prm", [] );
	[ ary3Kappa0, cfdjk_datOut ] = calcFDJKappa( vecX0, funchF, cfdjk_prm );
	vecF0 = cfdjk_datOut.vecF0;
	matJ0 = cfdjk_datOut.matJ0;
	omega0 = sumsq(vecF0,1)/2.0;
	vecG0 = cfdjk_datOut.vecG0;
	matH0 = cfdjk_datOut.matH0;
	%
	sizeF = size( vecF0, 1 );
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( isrealarray(matJ0,[sizeF,sizeX]) );
	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH0,[sizeX,sizeX]) );
	matIX = eye(sizeX,sizeX);
	%
	%
	% Do other prep.
	% ?
	%
	%
	% Set initial iterates here.
	bestIsBetterThan0 = false;
	vecX_best = vecX0;
	vecF_best = vecF0;
	omega_best = omega0;
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	ary3Kappa = ary3Kappa0;
	omega = omega0;
	vecG = vecG0;
	matH = matH0;
	iterCount = 0;
	while (1)
		%
		% Check pre-iter stop crit here.
		omegaMin = mygetfield( prm, "omegaMin", eps*omega0 );
		gNormMin = mygetfield( prm, "gNormMin", eps*norm(vecG0) );
		iterMax = mygetfield( prm, "iterMax", 1000 );
		if ( omega <= omegaMin )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: omega <= omegaMin." );
			break;
		elseif ( norm(vecG) <= gNormMin )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "WEAK SUCCESS: gNorm <= gNormMin." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		%
		%
		% Start iteration.
		iterCount++;
		%
		%
		%
		hNorm = sqrt( sum(sumsq(matH))/sizeX );
		assert( 0.0 < hNorm );
		findZero_fullH__regularize;
		% matHRegu is positive-definite and matR = chol(matHRegu).
		%
		%
		% Get scaling.
		matS_curve = mygetfield( prm, "matS_curve", [] );
		if ( isempty(matS_curve) )
			stepScalingType = mygetfield( prm, "stepScalingType", "" );
			switch ( tolower(stepScalingType) )
			case { "", "none", "identity" }
				matS_curve = matIX;
			case { "marquard" }
				matS_curve = diag(abs(matHRegu));
			otherwise
				error( "Invalid value of stepScalingType." );
			endswitch
		endif
		assert( isrealarray(matS_curve,[sizeX,sizeX]) );
		%
		% Set funchDeltaOfS here.
		stepCurveType = mygetfield( prm, "stepCurveType", "" );
		switch ( tolower(stepCurveType) )
		case { "newton" }
			error( "Not implemented." );
		case { "", "levenberg" }
			funchDeltaOfP = @(p) ( p*matHRegu + (1.0-p)*matS_curve ) \ (-p*vecG);
		case { "powell", "dog leg" }
			error( "Not implemented." );
		case { "gradesc", "gradient descent curve" }
			error( "Not implemented." );
		case { "cauchy", "gradient descent segment" }
			error( "Not implemented." );
		otherwise
			error( "Invalid case." );
		endswitch
		funchOmegaOfDelta = @(delta)( omega + vecG'*delta + 0.5*delta'*matH*delta );
		funchOmegaOfP = @(p) funchOmegaOfDelta( funchDeltaOfP(p) );
		%
		%
		%
		stepLengthType = mygetfield( prm, "stepLengthType", "" );
		switch ( tolower(stepLengthType) )
		case { "", "fminbnd", "min", "min scan" }
			pOfMin = fminbnd( funchOmegaOfP, 0.0, 1.0 );
			vecDelta = funchDeltaOfP( pOfMin );
		otherwise
			error( "Invalid case." );
		endswitch
		assert( isrealarray(vecDelta,[sizeX,1]) );
		vecX_trial = vecX + vecDelta;
		vecDelta = vecX_trial - vecX;
		assert( 0.0 ~= norm(vecDelta) );
		%
		error( "END OF VALID CODE, AS EXPECTED." );
	endwhile
	error( "END OF VALID CODE, SURPRISINGLY." );
	%
	vecXF = vecXBest;
return;
endfunction



%!test
%!	tic();
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 2;
%!	sizeF = 2;
%!	vecXE = randn(sizeX,1);
%!	matJE = randn(sizeF,sizeX);
%!	matA0 = 0.1*randn(sizeF,sizeX);
%!	matA1 = randn(sizeX,sizeX);
%!	matA2 = randn(sizeX,sizeX);
%!	matB0 = 0.01*randn(sizeF,sizeX);
%!	matB1 = randn(sizeX,sizeX);
%!	matB2 = randn(sizeX,sizeX);
%!	matB3 = randn(sizeX,sizeX);
%!	%
%!	y = @(x)( x - vecXE );
%!	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
%!	%funchF = @(x)( matJE * ( x - vecXE ) );
%!	%
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	[ vecXF, datOut ] = findZero_fullH( vecX0, funchF );
%!	%
%!	toc();

