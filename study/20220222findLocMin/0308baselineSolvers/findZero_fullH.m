% Function...

function [ vecXF, datOut ] = findZero_fullH( vecX0, funchF, prm=[] )
	% Look at basics.
	setCommon;
	fevalCount = 0;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__UNLIMITED );
	%
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
		norm(vecG)
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
			pOfMin = fminbnd( funchOmegaOfP, 0.0, 1.0 )
			vecDelta = funchDeltaOfP( pOfMin );
		otherwise
			error( "Invalid case." );
		endswitch
		assert( isrealarray(vecDelta,[sizeX,1]) );
		vecX_trial = vecX + vecDelta;
		vecDelta = vecX_trial - vecX;
		assert( 0.0 ~= norm(vecDelta) );
		omegaModel_trial = funchOmegaOfP( pOfMin )
		%
		vecF_trial = funchF( vecX_trial ); fevalCount++;
		omega_trial = sumsq(vecF_trial,1)/2.0;
		if ( omega_trial >= omega )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Failed to decrease omega." );
			break;
		endif
		%
		%
		vecX = vecX_trial;
		vecF = vecF_trial;
		[ ary3Kappa, cfdjk_datOut ] = calcFDJKappa( vecX, funchF, cfdjk_prm );
		omega = omega_trial;
		cfdjk_prm = mygetfield( prm, "cfdjk_prm", [] );
		cfdjk_prm.vecF0 = vecF;
		matJ = cfdjk_datOut.matJ0;
		vecG = cfdjk_datOut.vecG0;
		matH = cfdjk_datOut.matH0;
	endwhile
	%
	vecXF = vecX;
	if ( nargout >= 2 )
		datOut = [];
	endif
return;
endfunction


%!function vecF = funcFQuad( vecX, vecX0, vecF0, matJ0, ary3Kappa0 )
%!	sizeX = size(vecX0,1);
%!	assert( isrealarray(vecX0,[sizeX,1]) );
%!	assert( isrealarray(vecX,[sizeX,1]) );
%!	sizeF = size(vecF0,1);
%!	assert( isrealarray(vecF0,[sizeF,1]) );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	%
%!	vecD = vecX - vecX0;
%!	vecF = vecF0 + matJ0*vecD;
%!	for n=1:sizeF
%!		vecF(n) += 0.5*( vecD' * reshape( ary3Kappa0(n,:,:), [ sizeX, sizeX ] ) * vecD );
%!	endfor
%!endfunction


%!test
%!	setprngstates(0);
%!	sizeX = 2;
%!	sizeF = 2;
%!	vecX_secret = (1:sizeX)';%randn(sizeX,1);
%!	vecF_secret = zeros(sizeF,1);
%!	matJ_secret = eye(sizeF,sizeX);%randn(sizeF,sizeX);
%!	%
%!	sizeA = 2;
%!	for nf=1:sizeF
%!		matA = 0*randn(sizeA,sizeX);
%!		matKappa = matA'*matA;
%!		ary3Kappa_secret(nf,:,:) = matKappa;
%!	endfor
%!	%
%!	funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret );
%!	vecX0 = randn(sizeX,1);
%!	%
%!	[ vecXF, datOut ] = findZero_fullH( vecX0, funchF, prm=[] )
%!	rd = reldiff( vecXF, vecX_secret )
