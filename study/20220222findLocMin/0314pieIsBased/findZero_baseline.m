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
	% Set placeholder values for final log if don't complete an iter.
	vecX_prev = vecX0;
	vecF_prev = vecF0;
	matJ_prev = matJ0;
	omega_prev = omega0;
	vecG_prev = vecG0;
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0;
	vecG = vecG0;
	matJTJ = matJ0'*matJ0;
	trustRegionSize = [];
	iterCount = 0;
	%
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
		step_prm = mygetfield( prm, "step_prm" );
		step_prm.trustRegionSize = trustRegionSize;
		step_prm.funchFModel = funchFModel;
		[ vecX_step, vecF_step, step_datOut ] = findZero_baseline__step( vecX, funchDeltaOfP, funchF, step_prm );
		omega_step = sumsq(vecF_step)/2.0;
		%
		%
		omegaTol = mygetfield( prm, "omegaTol", eps^2*omega0 );
		if ( omega_step >= omega + omegaTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: omega_step >= omega + omegaTol." );
			break;
		elseif ( omega_step < omega_best )
			bestIsNot0 = true;
			vecX_best = vecX_step;
			vecF_best = vecF_step;
			omega_best = omega_step;
		endif
		%
		%
		iterCount++;
		%
		%
		iterMax = mygetfield( prm, "iterMax", 50 );
		if ( omega_step <= omegaTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: omega <= omegaTol." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		%
		% Move to next point.
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		omega_prev = omega;
		vecG_prev = vecG;
		matH_prev = matH;
		%
		vecX = vecX_step;
		vecF = vecF_step;
		omega = omega_step;
		%
		%
		%
		jupdate_prm = mygetfield( prm, "jupdate_prm", [] );
		[ matJ, jupdate_datOut ] = findZero_baseline__jupdate( vecX, vecF, vecX_prev, vecF_prev, matJ_prev, funchF, jupdate_prm );
		fevalCount += jupdate_datOut.fevalCount;
		vecG = matJ'*vecF;
		%
		%
		%
		msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e, %10.3e;  %10.3e, %10.3e;  %10.3e.", ...
		  time()-time0, iterCount, fevalCount, ...
		  norm(vecX-vecX0), norm(vecX-vecX_prev), ...
		  omega, omega_prev-omega, ...
		  norm(vecG) ) );
	endwhile
	msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e, %10.3e;  %10.3e, %10.3e;  %10.3e.", ...
	  time()-time0, iterCount, fevalCount, ...
	  norm(vecX-vecX0), norm(vecX-vecX_prev), ...
	  omega, omega_prev-omega, ...
	  norm(vecG) ) );
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
