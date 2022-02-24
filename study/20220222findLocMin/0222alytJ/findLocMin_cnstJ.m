% Function...
%  [ vecX, datOut ] = findLocMin_cnstJ( vecX0, vecF0, matJ, funchF, prm=[] )
% Searches for the local min of ||F|| starting assuming a constant Jacobian
%  but not assuming a constant Hessian.

function [ vecX, datOut ] = findLocMin_cnstJ( vecX0, vecF0, matJ, funchF, prm=[] )
	%
	%
	% Parse input.
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if ( debugMode )
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ,[sizeF,sizeX]) );
		%
		vecF0_test = funchF(vecX0);
		assert( isrealarray(vecF0_test,[sizeF,1]) );
		assert( reldiff(vecF0,vecF0_test) < sqrt(eps) );
		clear vecF0_test;
		%
		epsX_test = 0.001; % Arbitrary.
		vecDelta_test = epsX_test*ones(sizeX,1);
		vecJF_testA = ( funchF( vecX0 + vecDelta_test ) - funchF( vecX0 - vecDelta_test ) ) / 2.0;
		vecJF_testB = matJ*vecDelta_test;
		if ( sumsq(vecJF_testA-vecJF_testB) > sqrt(eps)*(sumsq(vecJF_testA)+sumsq(vecJF_testA)) )
			msg( __FILE__, __LINE__, "*** WARNING: Jacobian calculated by funchFJ appears to be incorrect. ***" );
		endif
		clear vecDelta_test;
		clear vecJF_testA;
		clear vecJF_testB;
	endif
	%
	matI = eye(sizeX,sizeX);
	omega0 = sumsq(vecF0,1)/2.0;
	vecG = matJ'*vecF0;
	matJTJ = matJ'*matJ;
	fevalCount = 0;
	%
	matK0 = zeros(sizeX,sizeX);
	iterLimit = 100;
	deltaNormMax = [];
	deltaNormMaxRelTol = 0.4;
	fallThresh_success = 0.9;
	fallThresh_giveup = 1.0e-5;
	fallThresh_okay = 1.0e-2;
	coeff_reduceTrustRegionOnFevalFail = 0.1;
	coeff_declareModelIsRadicallyWrong = 10.0;
	coeff_reduceTrustRegionOnRadicallyWrong = 0.2;
	doKUpdating = true;
	if ( ~isempty(prm) )
		matK0 = mygetfield( prm, "matK0", matK0 );
		deltaNormMax = mygetfield( prm, "deltaNormMax", deltaNormMax );
		doKUpdating = mygetfield( prm, "doKUpdating", doKUpdating );
	endif
	useDeltaNormMax = ~isempty(deltaNormMax);
	if ( debugMode )
		assert( isrealarray(matK0,[sizeX,sizeX]) );
		assert( issymmetric(matK0) );
		if ( useDeltaNormMax )
			assert( isrealscalar(deltaNormMax) );
			assert( 0.0 < deltaNormMax );
		endif
		assert( isscalar(doKUpdating) );
		assert( isbool(doKUpdating) );
	endif
	%
	%
	% Set some stuff...
	gNorm = norm(vecG);
	matI = eye(sizeX,sizeX);
	%
	% Set default return values.
	vecDelta = zeros(sizeX,1);
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	% Bail on "zero corner" cases.
	if ( 0.0 == omega0 )
		msg( __FILE__, __LINE__, "WARNING: vecF0 is zero." );
		return;
	elseif ( 0.0 == gNorm )
		msg( __FILE__, __LINE__, "WARNING: vecG is zero." );
		return;
	endif
	%
	%
	%
	vecX_best = vecX0;
	vecF_best = vecF0;
	omega_best = omega0;
	havePrev = false;
	vecX_prev = [];
	vecF_prev = [];
	omega_prev = [];
	havePrevPrev = false;
	vecX_prevPrev = [];
	vecF_prevPrev = [];
	omega_prevPrev = [];
	%
	vecX = vecX0;
	matK = matK0;
	iterCount = 0;
	haveNotedSepx = false;
	while ( 1 )
		%
		iterCount++;
		if ( iterCount > iterLimit )
			msg( __FILE__, __LINE__, "Imposed Stop: Reached iterLimit." );
			break;
		end
		%
		%
		% Check short-circuit success criteria.
		if ( omega0 - omega_best > (1.0-fallThresh_success)*(omega0-0.0) )
			% Since omegaModel can't be less than 0, we can check this before calculating the model.
			msgif( debugMode, __FILE__, __LINE__, "Success: omega_best is very good." );
			break;
		endif
		%
		% Calculate trial step.
		matH = matJTJ + matK;
		cdlPrm = [];
		cdlPrm.deltaNormMax = deltaNormMax; % May be empty.
		cdlPrm.deltaNormMaxRelTol = deltaNormMaxRelTol;
		[ vecDelta, cdlDatOut ] = calcDeltaLev( omega0, vecG, matH, cdlPrm )
		vecFModel = vecF0 + (matJ*vecDelta);
		omegaModel = cdlDatOut.omegaModel;
		%
		%
		% Decide whether or not to do a feval.
		if ( omega0 - omega_best > (1.0-fallThresh_success)*(omega0-omegaModel) )
			% omega_best offers a sufficient percentage of the potential decrease.
			msgif( debugMode, __FILE__, __LINE__, "Success: omega_best is sufficiently good compared to model." );
			break;
		elseif ( omegaModel > (1.0-fallThresh_giveup)*omega0 )
			msgif ( debugMode, __FILE__, __LINE__, "Imposed stop: omegaModel suggests very little decrease is possible." );
			break;
		elseif ( havePrevPrev )
			assert( havePrev );
			if ( omega_prev > omega_prevPrev )
				msgif( debugMode, __FILE__, __LINE__, "Objective previously got worse with backtracking." );
				if ( omega_prev > omega0 )
					msgif( ~haveNotedSepx, __FILE__, __LINE__, "  omega_prev > omega0; this suggests a sepratrix." );
					haveNotedSepx = true;
				elseif ( omega0 - omega_best > (1.0-fallThresh_okay)*(omega0-omegaModel) )
					msgif( debugMode, __FILE__, __LINE__, "Imposed Stop: omega_best is acceptable and further backtracking may hurt." );
					break;
				endif
			endif
		endif
		%
		%
		% Do the feval.
		vecX = vecX0 + vecDelta;
		vecF = funch( vecX );
		fevalCount++;
		%
		% If feval failed, cut trust region size.
		if ( ~isrealarray(vecF,[sizeF,1]) )
			msgif( debugMode, __FILE__, __LINE__, "Function evaluation failed." );
			if ( havePrev )
				msgif( ~haveNotedSepx, __FILE__, __LINE__, "  Function evaluation failed even though it succeeded for an earlier trial step!" );
				msgif( ~haveNotedSepx, __FILE__, __LINE__, "  This suggests the invalid space is not simply at distant values of x." );
				haveNotedSepx = true;
			endif
			deltaNormMax = coeff_reduceTrustRegionOnFevalFail * norm(vecDelta);
			msgif( debugMode, __FILE__, __LINE__, sprintf( "Set deltaNormMax to %g.", deltaNormMax ) );
			continue;
		endif
		%
		%
		% If vecF is radically different from vecFModel, treat like feval fail.
		% We could check to see how close vecF is to (1-b)*vecF0 + b*vecFModel,
		%  for the value of b in [0,1] that minimizes the difference, but,
		%  assuming norm( vecFModel - vecF0 ) <= norm( vecF0 ),
		%  and coeff_declareModelIsRadicallyWrong is > 2 ish(?), this isn't necessary.
		if ( norm( vecF - vecFModel ) >= coeff_declareModelIsRadicallyWrong * norm(vecF0) )
			msgif( debugMode, __FILE__, __LINE__, "Function was radically different at trial point." );
			if ( havePrev )
				msgif( ~haveNotedSepx, __FILE__, __LINE__, "  Function was radically different even though it was reasonable for an earlier trial step!" );
				msgif( ~haveNotedSepx, __FILE__, __LINE__, "  This suggests the presence of a sepratrix." );
				haveNotedSepx = true;
			endif
			deltaNormMax = coeff_reduceTrustRegionOnRadicallyWrong * norm(vecDelta);
			msgif( debugMode, __FILE__, __LINE__, sprintf( "Set deltaNormMax to %g.", deltaNormMax ) );
			continue;
		endif
		%
		%
		%
		omega = sumsq(vecF,1)/2.0
		if ( havePrev )
		if ( omega > omega_prev )
			msgif( ~haveNotedSepx, __FILE__, __LINE__, "omega > omega_prev; this suggests a sepratrix." );
			haveNotedSepx = true;
		endif
		endif
		%
		if ( doKUpdating )
			msgif( debugMode, __FILE__, __LINE__, "Updating K!" );
			matK += (2.0*(omega-omegaModel)/sumsq(vecDelta)^2)*(vecDelta*(vecDelta'));
		endif
		%
		%
		%
		if ( omega < omega_best )
			vecX_best = vecX;
			vecF_best = vecF;
			omega_best = omega;
		endif
		if ( havePrev )
			havePrevPrev = true;
			vecX_prevPrev = vecX_prev;
			vecF_prevPrev = vecF_prev;
			omega_prevPrev = omega_prev;
		endif
		havePrev = true;
		vecX_prev = vecX;
		vecF_prev = vecF;
		omega_prev = omega;
		%
		btForceFactor = 0.9; % Arbitrary
		deltaNormMax = btForceFactor*norm(vecDelta);
		msgif( debugMode, __FILE__, __LINE__, sprintf( "Set deltaNormMax to %g.", deltaNormMax ) );
	endwhile
	%
	%
	% Prep proper output.
	vecX = vecX_best;
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.vecF = vecF_best;
	endif
return;
endfunction
