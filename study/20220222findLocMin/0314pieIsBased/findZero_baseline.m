% Function...
%  findZero_baseline
% Focus is on "fully recalculated Jacobian",
% with varying degrees of support for other stuff.

function [ vecXF, datOut ] = findZero_baseline( vecX0, funchF, prm=[] )
	setCommon;
	fevalCount = 0;
	time0 = time();
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = funchF( vecX0 ); fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	omega0 = sumsq(vecF0)/2.0;
	%
	cfdj_prm = mygetfield( prm, "cdfj_prm", [] );
	cfdj_prm.vecF0 = vecF0;
	[ matJ0, cfdj_datOut ] = calcFDJ( vecX0, funchF, cfdj_prm );
	fevalCount += cfdj_datOut.fevalCount;
	assert( isrealarray(matJ0,[sizeF,sizeX]) );
	vecG0 = matJ0'*vecF0;
	assert( 0.0 ~= norm(vecG0) );
	%
	%
	bestIsNot0 = false;
	vecX_best = vecX0;
	vecF_best = vecF0;
	omega_best = omega0;
	%
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0;
	vecG = vecG0;
	matJTJ = matJ0'*matJ0;
	trustRegionSize = [];
	iterCount = 0;
	while (1)
		calcH_prm = mygetfield( prm, "calcH_prm", [] );
		[ matH, ch_datOut ] = findZero_baseline__calcH( vecX, vecF, matJ, funchF, calcH_prm );
		fevalCount += ch_datOut.fevalCount;
		funchFModel = @(x)( funcVecQuad( x, vecX, vecF, matJ, ch_datOut.ary3KappaA ) );
		funchFModelOfDelta = @(delta)( funcVecQuad( vecX+delta, vecX, vecF, matJ, ch_datOut.ary3KappaA ) );
		%
		%
		%
		scale_prm = mygetfield( prm, "scale_prm", [] );
		[ matSC, vecGSC, matHSC ] = findZero_baseline__scale( vecG, matH, scale_prm );
		%
		%
		%
		regu_prm = mygetfield( prm, "regu_prm", [] );
		matHSCRegu = calcHRegu( matHSC, regu_prm );
		%
		%
		%
		curve_prm = mygetfield( prm, "curve_prm", [] );
		[ funchDeltaOfP, curve_datOut ] = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, curve_prm );
		funchOmegaHessModelOfDelta = @(delta)( omega + vecG'*delta + 0.5*delta'*matH*delta );
		%
		%
		%
		doCurveViz = false;
		if (doCurveViz)
			findZero_baseline__curveViz;
			error( "Finizhed viz." );
		endif
		%
		%
		%
		%
		step_prm = mygetfield( prm, "step_prm" );
		step_prm.trustRegionSize = trustRegionSize;
		step_prm.funchFModel = funchFModel;
		[ vecX_step, vecF_step, step_datOut ] = findZero_baseline__step( vecX, funchDeltaOfP, funchF, step_prm )
		sumsq(vecF)/2.0
		sumsq(vecF_step)/2.0
		%
		msg( __FILE__, __LINE__, "HALT." ); error("HALT.");
		%
		%
		%
		calcNext_prm = mygetfield( prm, "calcNext_prm", [] );
		prm.trustRegionSize = trustRegionSize; %???
		[ p_trial, vecDelta_trial, vecF_trial, calcNext_datOut ] = findZero_baseline__calcNext( funchDeltaOfP, funchOmegaModelOfDelta, calcNext_prm );
		fevalCount += calcNext_datOut.fevalCount;
		trustRegionSize = calcNext_datOut.trustRegionSize; %???
		%
		%
		%
		error( "Pick step using TR or minscan or whatever." );
		error( "  Check continually: WEAK SUCCESS: omegaModelMin > omega - fTol*sqrt(omega)???." );
		error( "Update best." );
		iterCount++;
		error( "STRONG SUCCESS: omega < fTol^2." );
		error( "WEAK SUCCESS: fall < fTol*sqrt(omega)???" );
		error( "IMPOSED STOP: iterCount > iterLimit." );
		error( "Log progress." );
		error( "Update matJ." );
	endwhile
	msg( __FILE__, __LINE__, "HALT." ); error("HALT.");
	error( "Log final." );
	
	
	error( "END OF DEVELOPMENT CODE." );
	
	
	
	%
	cfdjk_prm = mygetfield( prm, "cfdjk_prm", [] );
	[ ary3Kappa0, cfdjk_datOut ] = calcFDJKappa( vecX0, funchF, cfdjk_prm );
	fevalCount += cfdjk_datOut.fevalCount;
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
	trustRegionSize = [];
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
			% The fact that this does not generate a point with fevalCount = 1 is irksome.
			% The real problem is how I constructed my loop.
			datOut.fevalCountVals(iterCount+1) = fevalCount;
			datOut.omegaVals(iterCount+1) = omega;
			datOut.gNormVals(iterCount+1) = norm(vecG);
		endif
		%
		%
		% Log progress.
		if ( verbLev >= VERBLEV__PROGRESS )
		%if ( abs( iterCount - round(sqrt(iterCount))^2 ) < 0.001 )
		if ( 0 == iterCount )
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d;  %10.3e, %10.3e;  %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, ...
			  sumsq(vecF)/2.0, -1.0, ...
			  norm(vecG) ) );
		else
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d;  %10.3e, %10.3e;  %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, ...
			  sumsq(vecF)/2.0, (sumsq(vecF_prev)-sumsq(vecF))/2.0, ...
			  norm(vecG) ) );
		endif
		%endif
		endif
		%
		% Check pre-iter stop crit.
		omegaTol = mygetfield( prm, "omegaTol", eps^2*omega0 );
		gNormMin = mygetfield( prm, "gNormMin", eps^2*norm(vecG0) );
		iterMax = mygetfield( prm, "iterMax", 50 );
		if ( omega <= omegaTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: omega <= omegaTol." );
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
		hNorm = sqrt( sum(sumsq(matH))/sizeX );
		assert( 0.0 < hNorm );
		%
		% Scale.
		matS_curve = mygetfield( prm, "matS_curve", [] );
		if ( isempty(matS_curve) )
			stepScalingType = mygetfield( prm, "stepScalingType", "" );
			switch ( tolower(stepScalingType) )
			case { "", "none", "identity" }
				matS_curve = matIX;
			case { "marquardt" }
				matS_curve = diag( 1.0./( abs(diag(matH)) + eps*hNorm ) );
			otherwise
				error( "Invalid value of stepScalingType." );
			endswitch
		endif
		assert( isrealarray(matS_curve,[sizeX,sizeX]) );
		vecGScaled = matS_curve'*vecG;
		matHScaled = matS_curve'*matH*matS_curve;
		%
		%
		%
		findZero_fullH__regularize;
		% matHRegu is positive-definite and matR_hRegu = chol(matHRegu).
		%
		% Set funchDeltaOfS here.
		vecDeltaNewton = matS_curve * ( matR_hRegu \ ( matR_hRegu'\(-vecGScaled) ) );
		pCauchy = calcLinishRootOfQuad( 0.5*(vecGScaled'*matHScaled*vecGScaled), -sumsq(vecGScaled), omega );
		assert( pCauchy > 0.0 );
		vecDeltaCauchy = pCauchy*matS_curve*(-vecGScaled);
		stepCurveType = mygetfield( prm, "stepCurveType", "" );
		switch ( tolower(stepCurveType) )
		case { "newton" }
			funchDeltaOfP = @(p) ( p * vecDeltaNewton );
		case { "", "levenberg" }
			funchDeltaOfP = @(p) matS_curve * (( p * matHRegu + (1.0-p)*eye(sizeX,sizeX) ) \ (-p*vecGScaled));
		case { "powell", "dog leg" }
			funchDeltaOfP = @(p) ( 2.0*p*vecDeltaCauchy + ...
			  (p>0.5) * ( (2.0*p-1.0)*vecDeltaNewton + 4.0*(0.5-p)*vecDeltaCauchy ) );
		case { "gradesc", "gradient descent curve" }
			% I suspect we could get a faster run time by doing an ODE solve then interpolating,
			% but, this is easier to code...
			[ matPsi_hRegu, matLambda_hRegu ] = eig( matHRegu );
			vecLambda_hRegu = diag(matLambda_hRegu);
			vecLIPNG = matLambda_hRegu \ ( matPsi_hRegu' * (-vecGScaled) );
			matSP = matS_curve * matPsi_hRegu;
			funchDeltaOfP = @(p) ( matSP * (diag( 1.0 - (1.0-p).^vecLambda_hRegu ) * vecLIPNG) );
		case { "cauchy", "gradient descent segment" }
			funchDeltaOfP = @(p) ( p * vecDeltaCauchy );
		otherwise
			error( "Invalid value of stepCurveType." );
		endswitch
		funchOmegaModelOfDelta = @(delta)( omega + vecG'*delta + 0.5*delta'*matH*delta );
		funchOmegaModelOfP = @(p)( funchOmegaModelOfDelta(funchDeltaOfP(p)) );
		%
		testPCauchy = true;
		if ( testPCauchy )
			omLo = funchOmegaModelOfDelta( 0.999 * vecDeltaCauchy );
			omAt = funchOmegaModelOfDelta( 1.000 * vecDeltaCauchy );
			omHi = funchOmegaModelOfDelta( 1.001 * vecDeltaCauchy );
			if ( abs(omAt) > eps*omega0  && omLo >= 0.0 && omAt >= 0.0  && omHi >= 0.0 )
				%[ omLo, omAt, omHi ]
				assert( omAt <= omLo );
				assert( omAt <= omHi );
			endif
		endif
		doCurveViz = false;
		if (doCurveViz)
			msg( __FILE__, __LINE__, "Doing curve viz hack." );
			numPts = 101;
			pPts = linspace(0.0,1.0,numPts);
			for n=1:numPts
				vecDeltaPts(:,n) = funchDeltaOfP( pPts(n) );
			endfor
			vecXPts = vecX + vecDeltaPts;
			figure(10);
			plot( vecXPts(1,:), vecXPts(2,:), 'o-' );
			grid on;
			hold on;
		break;
		endif
		%
		%
		if ( omega - funchOmegaModelOfP(1.0) < omegaTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: omega - funchOmegaModelOfP(1.0) < omegaTol." );
			break;
		endif
		%
		%
		%
		stepLengthType = mygetfield( prm, "stepLengthType", "" );
		switch ( tolower(stepLengthType) )
		case { "blind", "full" }
			vecDelta = funchDeltaOfP(1.0);
			vecX_trial = vecX + vecDelta;
			vecDelta = vecX_trial - vecX;
			omegaModel_trial = funchOmegaModelOfDelta(vecDelta);
			vecF_trial = funchF( vecX_trial ); fevalCount++;
			omega_trial = sumsq(vecF_trial,1)/2.0;
			%
			%
		case { "", "tr", "trust region", "bt", "backtracking" }
			switch (tolower(stepLengthType) )
			case { "bt", "backtracking" }
				% This is the only difference from trust region.
				trustRegionSize = [];
			endswitch
			matS_trustRegion = mygetfield( prm, "matS_trustRegion", eye(sizeX,sizeX) );
			funchStepSizeOfP = @(p)( norm(matS_trustRegion*funchDeltaOfP(p)) );
			%
			p = 1.0;
			btCount = 0;
			btLimit = mygetfield( prm, "btLimit", 20 );
			stepSizeMin = mygetfield( prm, "stepSizeMin", eps*(1.0+norm(vecX0)) );
			trThresh = mygetfield( prm, "trThresh", 1.1 );
			omThresh = mygetfield( prm, "omThresh", 1.1 );
			%
			failBTCoeff = mygetfield( prm, "failBTCoeff", 0.1 );
			uphillBTCoeff = mygetfield( prm, "uphillBTCoeff", 0.2 );
			accelTol = mygetfield( prm, "accelTol", 0.1 );
			accelCoeff = mygetfield( prm, "accelCoeff", 2.0 );
			while (1)
				%
				if ( ~isempty(trustRegionSize) )
				if ( trustRegionSize < stepSizeMin )
					msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: trustRegionSize < stepSizeMin." );
					break;
				endif
				endif
				%
				vecDelta = funchDeltaOfP(p);
				stepSize = funchStepSizeOfP(p);
				%
				if ( ~isempty(trustRegionSize) )
				if ( stepSize > trThresh*trustRegionSize )
					funchOverstepOfP = @(p)( funchStepSizeOfP(p) - trustRegionSize );
					p_next = fzero( funchOverstepOfP, [ 0.0, p ] );
					assert( p_next < p );
					p = p_next;
					continue;
				endif
				endif
				%
				if ( funchOmegaModelOfP(p) < (1.0-omThresh)*omega );
					msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "funchOmegaModelOfP is negative." );
					p_next = fzero( funchOmegaModelOfP, [ 0.0, p ] );
					assert( p_next < p );
					p = p_next;
					if (isempty(trustRegionSize))
						trustRegionSize = funchStepSizeOfP(p);
					elseif ( funchStepSizeOfP(p) < trustRegionSize )
						trustRegionSize = funchStepSizeOfP(p);
					endif
					continue;
				endif
				%
				vecX_trial = vecX + vecDelta;
				vecDelta = vecX_trial - vecX;
				omegaModel_trial = funchOmegaModelOfDelta(vecDelta);
				%
				%
				if ( omega - omegaModel_trial < omegaTol )
					msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: omega - omegaModel_trial < omegaTol." );
					break;
				endif
				%
				vecF_trial = funchF( vecX_trial ); fevalCount++;
				if ( ~isrealarray(vecF_trial,[sizeF,1]) )
					msg( __FILE__, __LINE__, "feval failed." );
					trustRegionSize = stepSize * failBTCoeff;
					continue;
				endif
				%
				omega_trial = sumsq(vecF_trial,1)/2.0;
				if ( omega_trial >= omega )
					trustRegionSize = stepSize * uphillBTCoeff;
					continue;
				endif
				%
				if ( ~isempty(trustRegionSize) && reldiff(omega_trial,omegaModel_trial) < accelTol && 1 >= btCount )
					msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Increasing trust region size." );
					trustRegionSize = stepSize * accelCoeff;
				endif
				break;
			endwhile
			%
			%
		case { "fminbnd", "min", "minscan", "min scan" }
			% Specifically checking p=1.0 seems important for gradesc. Not sure fhat fminbnd() is doing.
			vecDelta_full = funchDeltaOfP(1.0);
			vecX_full = vecX + vecDelta_full;
			vecDelta_full = vecX_full - vecX;
			omegaModel_full = funchOmegaModelOfDelta(vecDelta_full);
			vecF_full = funchF( vecX_full ); fevalCount++;
			omega_full = sumsq(vecF_full,1)/2.0;
			%
			fminbnd_prm = optimset( "TolX", eps^2, "TolFun", eps^2 );
			funchOmegaOfP = @(p)( sumsq(funchF(vecX+funchDeltaOfP(p)),1)/2.0 );
			[ p_scan, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchOmegaOfP, 0.0, 1.0, fminbnd_prm );
			fevalCount += fminbnd_output.funcCount;
			vecDelta_scan = funchDeltaOfP( p_scan );
			vecX_scan = vecX + vecDelta_scan;
			vecDelta_scan = vecX_scan - vecX;
			vecF_scan = funchF( vecX_scan ); fevalCount++;
			omega_scan = sumsq(vecF_scan,1)/2.0;
			%
			%
			if ( omega_scan < omega_full )
				vecDelta = vecDelta_scan;
				vecX_trial = vecX_scan;
				vecF_trial = vecF_scan;
			else
				vecDelta = vecDelta_full;
				vecX_trial = vecX_full;
				vecF_trial = vecF_full;
			endif
			omegaModel_trial = funchOmegaModelOfDelta( vecDelta );
			omega_trial = sumsq(vecF_trial)/2.0;
			%
			%
		otherwise
			error( "Invalid value of stepLengthType." );
		endswitch
		if ( 0.0 == norm(vecDelta) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: 0.0 == norm(vecDelta)." );
			break;
		endif
		omegaModel_trial = funchOmegaModelOfDelta( vecDelta );
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
		vecX_prev = vecX;
		vecF_prev = vecF;
		%
		vecX = vecX_trial;
		vecF = vecF_trial;
		cfdjk_prm = mygetfield( prm, "cfdjk_prm", [] );
		cfdjk_prm.vecF0 = vecF;
		[ ary3Kappa, cfdjk_datOut ] = calcFDJKappa( vecX, funchF, cfdjk_prm );
		fevalCount += cfdjk_datOut.fevalCount;
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


%!test
%!	tic();
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 5;
%!	sizeF = 5;
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
%!	[ vecXF, datOut ] = findZero_baseline( vecX0, funchF );
