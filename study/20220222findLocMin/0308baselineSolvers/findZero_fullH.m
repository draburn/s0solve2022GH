% Function...

function [ vecXF, datOut ] = findZero_fullH( vecX0, funchF, prm=[] )
	% Look at basics.
	setCommon;
	fevalCount = 0;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
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
	% Any other prep?
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	%
	% Set initial iterates here.
	bestIsNot0 = false;
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
	%echo__normG0 = norm(vecG0)
	while (1)
		if ( nargout >= 2 )
			datOut.fevalCountVals(iterCount+1) = fevalCount;
			datOut.omegaVals(iterCount+1) = omega;
			datOut.gNormVals(iterCount+1) = norm(vecG);
		endif
		%
		% Check pre-iter stop crit.
		omegaMin = mygetfield( prm, "omegaMin", eps^2*omega0 );
		gNormMin = mygetfield( prm, "gNormMin", eps^2*norm(vecG0) );
		iterMax = mygetfield( prm, "iterMax", 20 );
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
		% matHRegu is positive-definite and matR_hRegu = chol(matHRegu).
		%
		%
		% Get scaling.
		matS_curve = mygetfield( prm, "matS_curve", [] );
		if ( isempty(matS_curve) )
			stepScalingType = mygetfield( prm, "stepScalingType", "" );
			switch ( tolower(stepScalingType) )
			case { "", "none", "identity" }
				matS_curve = matIX;
			case { "marquardt" }
				matS_curve = diag(abs(matHRegu));
			otherwise
				error( "Invalid value of stepScalingType." );
			endswitch
		endif
		assert( isrealarray(matS_curve,[sizeX,sizeX]) );
		%
		% Set funchDeltaOfS here.
		vecDeltaNewton = matR_hRegu \ ( matR_hRegu' \ (-vecG) );
		pCauchy = calcLinishRootOfQuad( omega, -sumsq(vecG), 0.5*(vecG'*matH*vecG) );
		vecDeltaCauchy = pCauchy*(-vecG);
		stepCurveType = mygetfield( prm, "stepCurveType", "" );
		switch ( tolower(stepCurveType) )
		case { "newton" }
			funchDeltaOfP = @(p) ( p * vecDeltaNewton );
		case { "", "levenberg" }
			funchDeltaOfP = @(p) ( p*matHRegu + (1.0-p)*matS_curve ) \ (-p*vecG);
		case { "powell", "dog leg" }
			funchDeltaOfP = @(p) ( 2.0*p*vecDeltaCauchy + ...
			  (p>0.5) * ( (2.0*p-1.0)*vecDeltaNewton + 4.0*(0.5-p)*vecDeltaCauchy ) );
		case { "gradesc", "gradient descent curve" }
			error( "Not implemented." );
		case { "cauchy", "gradient descent segment" }
			funchDeltaOfP = @(p) ( p * vecDeltaCauchy );
		otherwise
			error( "Invalid value of stepCurveType." );
		endswitch
		funchOmegaModelOfDelta = @(delta)( omega + vecG'*delta + 0.5*delta'*matH*delta );
		funchOmegaOfDelta = @(delta)( sumsq(funchF(vecX+delta),1)/2.0 );
		funchOmegaOfP = @(p) funchOmegaOfDelta( funchDeltaOfP(p) );
		%
		%
		%
		stepLengthType = mygetfield( prm, "stepLengthType", "" );
		switch ( tolower(stepLengthType) )
		case { "tr", "trust region" }
			error( "Not implemented." );
		case { "", "fminbnd", "min", "min scan" }
			fminbnd_prm = optimset( "TolX", eps^2, "TolFun", eps^2 );
			[ pOfMin, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchOmegaOfP, 0.0, 1.0, fminbnd_prm );
			fevalCount += fminbnd_output.funcCount;
			vecDelta = funchDeltaOfP( pOfMin );
		otherwise
			error( "Invalid value of stepLengthType." );
		endswitch
		assert( isrealarray(vecDelta,[sizeX,1]) );
		vecX_trial = vecX + vecDelta;
		vecDelta = vecX_trial - vecX;
		if ( 0.0 == norm(vecDelta) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: 0.0 == norm(vecDelta)." );
			break;
		endif
		omegaModel_trial = funchOmegaModelOfDelta( vecDelta );
		%
		vecF_trial = funchF( vecX_trial ); fevalCount++;
		omega_trial = sumsq(vecF_trial,1)/2.0;
		if ( omega_trial >= omega )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: omega_trial >= omega." );
			break;
		endif
		if ( omega_trial < omega_best )
			assert( 0.0 < norm(vecX_trial-vecX_best) );
			omega_best = omega_trial;
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			bestIsNot0 = true;
		endif
		%
		%
		vecX = vecX_trial;
		vecF = vecF_trial;
		cfdjk_prm = mygetfield( prm, "cfdjk_prm", [] );
		cfdjk_prm.vecF0 = vecF;
		[ ary3Kappa, cfdjk_datOut ] = calcFDJKappa( vecX, funchF, cfdjk_prm );
		omega = omega_trial;
		matJ = cfdjk_datOut.matJ0;
		vecG = cfdjk_datOut.vecG0;
		matH = cfdjk_datOut.matH0;
		%echo__normG = norm(vecG)
	endwhile
	%
	vecXF = vecX_best;
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.iterCount = iterCount;
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


%!test
%!	tic();
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 20;
%!	sizeF = 20;
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
%!	%vecX0 = vecXE + 0.001*randn(sizeX,1)
%!	%
%!	[ vecXF, datOut ] = findZero_fullH( vecX0, funchF );
%!	%
%!	msg( __FILE__, __LINE__, "*** PLEASE USE findZero_fullH__test INSTEAD! ***" );
%!	toc();
