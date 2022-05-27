function [ vecX_best, vecF_best, retCode, fevalCount, stepsCount, datOut ] = zlinsolf200( funchF, vecX_initial, vecF_initial=[], prmIn=[] )
	startTime = time();
	if ( stopsignalpresent() )
		msg(__FILE__, __LINE__, "ERROR: Stop signal already present." );
		retCode = RETCODE__IMPOSED_STOP;
		break;
	endif
	% INIT
	mydefs;
	vecX_best = [];
	vecF_best = [];
	retCode = RETCODE__NOT_SET;
	fevalCount = 0;
	stepsCount = 0;
	datOut = [];
	%
	[ retCode, fevalIncr, vecF_initial, fModelDat, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn );
	fevalCount += fevalIncr; clear fevalIncr;
	if ( 0~= retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	vecF_best = vecF_initial;
	%
	%
	% MAIN LOOP
	vecX = vecX_initial;
	vecF = vecF_initial;
	vecX_next = [];
	vecF_next = [];
	iterCount = 0;
	while (1)
		% Tier 0a actions: simple stoping criteria.
		if ( norm(vecF_best) <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: norm(vecF_best) <= prm.fTol." );
			retCode = RETCODE__SUCCESS;
			break;
		elseif ( prm.timeMax >= 0.0 && time()-startTime >= prm.timeMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: time()-startTime >= prm.timeMax." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		elseif ( prm.iterMax >= 0 & iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterMax." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		elseif ( prm.fevalMax >= 0 && fevalCount >= prm.fevalMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: fevalCount >= prm.fevalMax." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		elseif ( prm.stepsMax >= 0 && stepsCount >= prm.stepsMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: stepsCount >= prm.stepsMax." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		elseif ( stopsignalpresent() )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		endif
		%
		%
		iterCount++;
		%
		%
		% Tier 0b actions: mandatory actions.
		if ( isempty(fModelDat) )
			[ retCode, fevalIncr, fModelDat ] = __initModel( funchF, vecX, vecF, prm );
			fevalCount += fevalIncr; clear fevalIncr;
			if ( 0~=retCode )
				msg( __FILE__, __LINE__, "Passing an error." );
				break;
			endif
			continue;
		endif
		if ( prm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, "Progress..." );
			msg( __FILE__, __LINE__, sprintf( ...
			  "  elapsed time: %10.3e/%0.3e;  iter: %5d/%d;  feval: %5d/%d;  steps: %5d/%d;  size: %5d/%5d/%dx%d", ...
			  time()-startTime, prm.timeMax, ...
			  iterCount, prm.iterMax, ...
			  fevalCount, prm.fevalMax, ...
			  stepsCount, prm.stepsMax, ...
			  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX,1), size(vecF,1) ) );
		endif
		[ retCode, fevalIncr, fModelDat ] = __analyzeModel( funchF, fModelDat, prm );
		fevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		%
		%
		% Tier 1 actions: rather clear-cut.
		
		% If we have a "next" candidate and it looks sufficiently good, move to it.
		
		% If high variation omega of zero variation (bounded) step is below cnvg tol, try it.
		
		% If zero variation omega of (zero variation) unbounded step is large, expand subspace.
		% We'll be particularly inclined to do this if we have yet to take a step.
		% But we cannot do this is we've already explored the full X space.
		if ( 0 == stepsCount )
		if ( fModelDat.etaZeroV_ideal > sqrt( prm.omegaTol*sumsq(vecF)/2.0 ) )
			[ retCode, fevalIncr, fModelDat ] = __expandModel( funchF, fModelDat, prm );
			fevalCount += fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				break
			endif
			continue;
		endif
		endif
		
		% If zero variation omega of zero variation (bounded) step offers only a too small decrease,
		%  declare that the trust region size has gotten too small.
		
		% If hivar omega of hivar (b) step offer suffic reduction, try it.
		
		msgif( prm.verbLev >= VERBLEV__ERROR, __FILE__, __LINE__, "ERROR: No acceptable action." );
		break;
	endwhile
	if ( prm.verbLev >= VERBLEV__PROGRESS )
		msg( __FILE__, __LINE__, "Final..." );
		msg( __FILE__, __LINE__, sprintf( ...
		  "  elapsed time: %10.3e/%10.3e;  iter: %5d/%5d;  feval: %5d/%5d;  steps: %5d/%5d.", ...
		  time()-startTime, prm.timeMax, iterCount, prm.iterMax, fevalCount, prm.fevalMax, stepsCount, prm.stepsMax ) );
	endif
return;
endfunction


function [ retCode, fevalIncr, vecF_initial, fModelDat, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn )
	% If an error is thrown during init, let it be thrown for useful reporting to the command window.
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	vecF_initial = [];
	prm = [];
	%
	if ( isempty(vecF_initial) )
		fevalIncr++;
		vecF_initial = funchF(vecX_initial);
	endif
	sizeX = size(vecX_initial,1);
	sizeF = size(vecF_initial,1);
	%
	prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__ZERO; % Production.
	prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration.
	prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Performance testing.
	prm.verbLev = VERBLEV__COPIOUS; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Debug.
	%
	%%%prm.timeMax = -1.0; %%%
	prm.timeMax = 3.0; %%%
	prm.iterMax = ceil( 100 + 10*sqrt(sizeX+sizeF) + sizeX );
	prm.fevalMax = prm.iterMax;
	prm.stepsMax = 100;
	prm.fTol = sizeF*100.0*eps;
	%
	prm.fModelDat_initial = [];
	%
	prm.precon_funchPrecon = [];
	prm.precon_matL = [];
	prm.precon_matU = [];
	prm.precon_matJA = [];
	%
	prm.epsFD = 1.0e-3;
	prm.orthoTol = 1.0e-10;
	%
	prm.epsB = sqrt(eps);
	prm.curveType = "lev"; %%%
	%%%prm.curveType = "powell"; %%%
	prm.curveScaling = "b";
	prm.matC = [];
	prm.cholRelTol = sqrt(eps);
	prm.epsRegu = sqrt(eps);
	%%%prm.candStepRelTol = 0.2;
	prm.candStepRelTol = 1000.0*eps;
	prm.levIterMax = 100;
	%
	prm = overwritefields( prm, prmIn );
	%
	if ( ~isempty(prm.precon_matJA) )
		if ( isempty(prm.precon_matL) && isempty(prm.precon_matU) )
			[ prm.precon_matL, prm.precon_matU ] = lu( prm.precon_matJA );
		else
			msgif( prm.verbLev >= VERBLEV__ERROR, "ERROR: prm.precon_matJA is non-empty but at least one of _matL and _matU is also non-empty." );
			retCode = RETCODE__BAD_INPUT;
			return;
		endif
	endif
	prm.omegaTol = (prm.fTol^2)/2.0;
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( isrealscalar(prm.verbLev) );
		assert( isrealscalar(prm.valdLev) );
		%
		assert( isrealarray(vecX_initial,[sizeX,1]) );
		assert( isrealarray(vecF_initial,[sizeF,1]) );
		%
		assert( isrealscalar(prm.timeMax) );
		assert( isrealscalar(prm.iterMax) );
		assert( isrealscalar(prm.fevalMax) );
		assert( isrealscalar(prm.stepsMax) );
		assert( abs(prm.iterMax-round(prm.iterMax)) < sqrt(eps) );
		assert( abs(prm.fevalMax-round(prm.fevalMax)) < sqrt(eps) );
		assert( abs(prm.stepsMax-round(prm.stepsMax)) < sqrt(eps) );
		%
		assert( isrealscalar(prm.fTol) );
		assert( isrealscalar(prm.epsFD) );
		assert( isrealscalar(prm.orthoTol) );
		assert( 0.0 < prm.fTol );
		assert( 0.0 < prm.epsFD );
		assert( 0.0 < prm.orthoTol );
		%
		if ( ~isempty(prm.precon_funchPrecon) )
			assert( isempty(prm.precon_matL) );
			assert( isempty(prm.precon_matU) );
			assert( isempty(prm.precon_matJA) );
		elseif ( ~isempty(prm.precon_matL) )
			assert( isrealarray(prm.precon_matL,[sizeF,min([sizeX,sizeF])]) );
			assert( ~isempty(prm.precon_matU) );
			assert( isrealarray(prm.precon_matL,[min([sizeX,sizeF]),sizeX]) );
			% Allow precon_matJA to be whatever.
		elseif ( isempty(prm.precon_matL) )
			assert( isempty(prm.precon_matU) );
		endif
	endif
	%
	fModelDat = prm.fModelDat_initial;
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __initModel( funchF, vecX, vecF, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	fModelDat = [];
	%
	vecRhoF = vecF;
	vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, [], prm );
	if ( norm(vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "__applyPrecon() failed to generate a linearly independent vector." );
		sizeX = size(vecX,1);
		for n=1:sizeX
			vecU = zeros(size(vecX));
			vecU(n) = 1.0;
			vecV = __calcOrthonorm( vecU, [], prm );
			if ( norm(vecV) > sqrt(eps) )
				break;
			endif
		endfor
	endif
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalIncr++;
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = 1;
	sizeVLocal = 1;
	sizeB = 1;
	%
	fModelDat.matV = [ vecV ];
	fModelDat.matW_frozen = [ vecW ]; % Only for dev/debug?
	fModelDat.matALo_frozen = [ 0.0 ]; % Only for dev/debug?
	fModelDat.matAHi_frozen = [ 0.0 ]; % Only for dev/debug?
	fModelDat.matB_frozen = [ 0.0 ]; % Only for dev/debug?
	%
	fModelDat.matW = [ vecW ];
	fModelDat.matALo = [ 0.0 ];
	fModelDat.matAHi = [ 0.0 ];
	fModelDat.matB = [ 0.0 ];
	%
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matVLocal = [ vecV ];
	fModelDat.matWLocal = [ vecW ]; % Only for dev/debug?
	%
	fModel.sizeX = size(vecX,1);
	fModel.sizeF = size(vecF,1);
	fModel.sizeV = 1;
	fModel.sizeVLocal = 1;
	fModel.sizeB = 1;
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function vecU = __applyPrecon( vecRhoF, prm, vecX, vecF )
	if ( ~isempty(prm.precon_funchPrecon) )
		vecU = prm.precon_funchPrecon( vecRhoF, vecX, vecF );
	elseif ( ~isempty(prm.precon_matU) )
		vecU = prm.precon_matU \ ( prm.precon_matL \ vecRhoF );
	elseif ( size(vecX,1) == size(vecF,1) )
		vecU = vecRhoF;
	else
		sizeX = size(vecX,1);
		sizeF = size(vecF,1);
		vecU = interp1( (0:sizeF-1), vecRhoF, linspace(0.0,sizeF-1.0,sizeX)' );
	endif
	return;
endfunction


function vecV = __calcOrthonorm( vecU, matV, prm )
	u0 = norm(vecU);
	if (0.0==u0)
		vecV = zeros(size(vecU));
		return;
	elseif (isempty(matV))
		vecV = vecU/u0;
		return;
	endif
	vecV = vecU;
	for n=1:2
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( v <= prm.orthoTol*u0 )
			vecV(:) = 0.0;
			return;
		else
			vecV /= v;
		endif
	endfor
	return;
endfunction


function vecW = __calcJV( vecV, funchF, vecX, vecF, prm )
	mydefs;
	v = norm(vecV);
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < v );
	endif
	vecFP = funchF( vecX + prm.epsFD*vecV );
	vecW = ( vecFP - vecF ) / (prm.epsFD * v);
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __analyzeModel( funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	matV = fModelDat.matV; % Subspace basis matrix.
	matW_frozen = fModelDat.matW_frozen; % Only for dev/debug?
	matALo_frozen = fModelDat.matALo_frozen; % Only for dev/debug?
	matAHi_frozen = fModelDat.matAHi_frozen; % Only for dev/debug?
	matB_frozen = fModelDat.matB_frozen; % Only for dev/debug?
	%
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Only for dev/debug?
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeB == sizeV ); % Not strictly required, but it's what we're doing.
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW_frozen,[sizeF,sizeV]) );
		assert( isrealarray(matALo_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matAHi_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matB_frozen,[sizeB,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( isrealarray(matWLocal,[sizeF,sizeVLocal]) );
		%
		assert( issymmetric(fModelDat.matALo) );
		assert( issymmetric(fModelDat.matAHi) );
		assert( issymmetric(fModelDat.matALo_frozen) );
		assert( issymmetric(fModelDat.matAHi_frozen) );
		% We're making matB sym too...
		assert( issymmetric(fModelDat.matB) );
		assert( issymmetric(fModelDat.matB_frozen) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(fModelDat.matV'*fModelDat.matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(fModelDat.matVLocal'*fModelDat.matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
	endif
	%
	matH = matW'*matW;
	vecG = matW'*vecF;
	matBTB = matB'*matB;
	matIV = eye(sizeV,sizeV);
	% "eta" is estimate for cost function;
	% "omega" is observed values.
	funchEta_zeroV = @(y)( sumsq(vecF)/2.0 + vecG'*y + (y'*matH*y)/2.0 );
	funchEta_loVar = @(y)( funchEta_zeroV(y) + (y'*matALo*y)/2.0 );
	funchEta_hiVar = @(y)( funchEta_zeroV(y) + (y'*matAHi*y)/2.0 );
	%
	if ( ~isempty(prm.matC) )
		matC = prm.matC;
	else
		switch ( tolower(prm.curveScaling) )
		case { "1", "eye", "i" }
			matC = matIV;
		case { "b", "btb", "boundary", "optimal" }
			% This is "optimal" in that, used with the Levenberg curve,
			%  it will produce the point on the boundary which minimizes
			%  the (estimated) objective function.
			matC = matBTB;
		case { "ddbtb" }
			matC = diag(diag(matBTB));
		case { "wtw", "newton" }
			matC = matH;
		case { "m", "marq", "marquardt", "ddwtw" }
			matC = diag(diag(matH));
		otherwise
			error( "Invalid value of curveScaling." );
		endswitch
		cScale = norm(diag(matC));
		if ( 0.0 == cScale )
			msgif( prm.verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, "Curve scaling matrix was zero; setting to I." );
			matC = matIV;
		else
			[ matRC, cholFlag ] = chol(matC);
			if ( 0 ~= cholFlag || min(diag(matRC)) < prm.cholRelTol * max(abs(diag(matRC))) )
				msgif( prm.verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Curve scaling matrix was non positive-definite; applying regularization." );
				matC += cScale * prm.epsRegu * matIV;
			endif
			cScale = norm(diag(matC));
			matRC = chol(matC);
		endif
	endif
	%
	vecY_ideal = __findCandStep( matH, vecG, matC, [], [], prm );
	vecY_zeroV = __findCandStep( matH, vecG, matC, matB, 1.0, prm );
	vecY_loVar = __findCandStep( matH + matALo, vecG, matC, matB, 1.0, prm );
	vecY_hiVar = __findCandStep( matH + matAHi, vecG, matC, matB, 1.0, prm );
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		assert( norm(matB*vecY_zeroV) <= 1.0 + prm.candStepRelTol );
		assert( norm(matB*vecY_loVar) <= 1.0 + prm.candStepRelTol );
		assert( norm(matB*vecY_hiVar) <= 1.0 + prm.candStepRelTol );
		if ( prm.candStepRelTol <= 1000.0*eps );
			assert( funchEta_zeroV(vecY_ideal) <= (1.0+sqrt(eps))*abs(funchEta_zeroV(vecY_zeroV)) + eps*sumsq(vecF) );
			assert( funchEta_zeroV(vecY_zeroV) <= (1.0+sqrt(eps))*abs(funchEta_zeroV(vecY_loVar)) + eps*sumsq(vecF) );
			assert( funchEta_zeroV(vecY_loVar) <= (1.0+sqrt(eps))*abs(funchEta_zeroV(vecY_hiVar)) + eps*sumsq(vecF) );
			assert( funchEta_zeroV(vecY_hiVar) <= (1.0+sqrt(eps))*sumsq(vecF)/2.0 );
			assert( funchEta_loVar(vecY_loVar) <= (1.0+sqrt(eps))*abs(funchEta_loVar(vecY_zeroV)) + eps*sumsq(vecF) );
			assert( funchEta_loVar(vecY_loVar) <= (1.0+sqrt(eps))*abs(funchEta_loVar(vecY_hiVar)) + eps*sumsq(vecF) );
			assert( funchEta_hiVar(vecY_hiVar) <= (1.0+sqrt(eps))*abs(funchEta_hiVar(vecY_zeroV)) + eps*sumsq(vecF) );
			assert( funchEta_hiVar(vecY_hiVar) <= (1.0+sqrt(eps))*abs(funchEta_hiVar(vecY_loVar)) + eps*sumsq(vecF) );
		endif
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
	switch ( tolower(prm.curveScaling) )
	case { "b", "btb", "boundary", "optimal" }
	%%%otherwise %%%
	switch ( prm.curveType )
	case { "l", "lev", "levenberg" }
	%%%otherwise %%%
		state0 = rand( "state" );
		rand( "state", 0 );
		for m = 1 : 4
			switch (m)
			case 1
				vecY = vecY_ideal;
				funchEta = funchEta_zeroV;
			case 2
				vecY = vecY_zeroV;
				funchEta = funchEta_zeroV;
			case 3
				vecY = vecY_hiVar;
				funchEta = funchEta_hiVar;
			case 4
				vecY = vecY_loVar;
				funchEta = funchEta_loVar;
			otherwise
				error( "Invalid case." );
			endswitch
			b = norm(matB*vecY);
			eta = funchEta( vecY );
			yNorm = norm(vecY);
			vecY_temp = vecY + yNorm * sqrt(eps)*(2.0*rand(sizeV,1)-1.0);
			b_temp = norm(matB*vecY_temp);
			eta_temp = funchEta(vecY_temp);
			assert( eta_temp >= eta-1000.0*eps*eta || b_temp >= b-1000.0*eps*b );
			eta_temp_etaTempMin = eta_temp;
			b_temp_etaTempMin = b_temp;
			eta_temp_bTempMin = eta_temp;
			b_temp_bTempMin = b_temp;
			for n=1:100
				vecY_temp = vecY + yNorm * sqrt(eps)*(2.0*rand(sizeV,1)-1.0);
				b_temp = norm(matB*vecY_temp);
				eta_temp = funchEta(vecY_temp);
				assert( eta_temp >= eta-1000.0*eps*eta || b_temp >= b-1000.0*eps*b );
			endfor
			clear eta_temp;
			clear b_temp;
			clear vecY_temp;
			clear yNorm;
			clear eta;
			clear b;
			clear funchEta;
			clear vecY;
			clear n;
		endfor
		clear m;
		rand( "state", state0 );
	endswitch
	endswitch
	endif
	%
	fModelDat.vecY_ideal = vecY_ideal;
	fModelDat.vecY_zeroV = vecY_zeroV;
	fModelDat.vecY_loVar = vecY_loVar;
	fModelDat.vecY_hiVar = vecY_hiVar;
	%
	fModelDat.etaZeroV_ideal = funchEta_zeroV( vecY_ideal );
	fModelDat.etaZeroV_zeroV = funchEta_zeroV( vecY_zeroV );
	fModelDat.etaZeroV_loVar = funchEta_zeroV( vecY_loVar );
	fModelDat.etaZeroV_hiVar = funchEta_zeroV( vecY_hiVar );
	%
	fModelDat.etaLoVar_ideal = funchEta_loVar( vecY_ideal );
	fModelDat.etaLoVar_zeroV = funchEta_loVar( vecY_zeroV );
	fModelDat.etaLoVar_loVar = funchEta_loVar( vecY_loVar );
	fModelDat.etaLoVar_hiVar = funchEta_loVar( vecY_hiVar );
	%
	fModelDat.etaHiVar_ideal = funchEta_hiVar( vecY_ideal );
	fModelDat.etaHiVar_zeroV = funchEta_hiVar( vecY_zeroV );
	fModelDat.etaHiVar_loVar = funchEta_hiVar( vecY_loVar );
	fModelDat.etaHiVar_hiVar = funchEta_hiVar( vecY_hiVar );
	%
	fModelDat.expand_vecU_suggested = vecF - matW*vecY_ideal;
	fModelDat.refresh_vecY_suggested = vecY_loVar;
	fModelDat.tryStep_vecY_suggested = vecY_hiVar;
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function vecY = __findCandStep( matH, vecG, matC, matB, bTrgt, prm )
	mydefs;
	if ( prm.valdLev >= VALDLEV__HIGH )
		sz = size(matH,1);
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH) );
		assert( isrealarray(vecG,[sz,1]) );
		assert( isrealarray(matC,[sz,sz]) );
		assert( issymmetric(matC) );
		if ( ~isempty(matB) )
			szb = size(matB,1);
			assert( isrealarray(matB,[szb,sz]) );
			clear szb;
			assert( ~isempty(bTrgt) );
		endif
		assert( min(diag(matH)) >= 0.0 );
		assert( min(diag(matC)) > 0.0 );
		if ( ~isempty(bTrgt) )
			assert( ~isempty(matB) );
			assert( bTrgt > 0.0 );
		endif
		clear sz;
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
		eigH = eig(matH);
		msgif( prm.verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, sprintf( "eig(matH): %0.3e ~ %0.3e", min(eigH), max(eigH) ) );
		eigC = eig(matC);
		msgif( prm.verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, sprintf( "eig(matC): %0.3e ~ %0.3e", min(eigC), max(eigC) ) );
		%
		if ( min(eigH) < -sqrt(eps)*max(abs(eigH)) )
			error( "Hessian matrix has a clearly negative eigenvalue." );
		elseif ( min(eigC) < -sqrt(eps)*max(abs(eigC)) )
			error( "Constraint matrix has a clearly negative eigenvalue." );
		endif
		if ( min(eigC) <= 0.0 )
			msgif ( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, ...
			  "WARNING: Constraint matrix appears to have a near-zero eigenvalue." );
		endif
		clear eigC;
		clear eigH;
	endif
	%
	switch ( prm.curveType )
	case { "l", "lev", "levenberg" }
		vecY = __findCandStep_lev( matH, vecG, matC, matB, bTrgt, prm );
	case { "p", "powell", "dog", "dog leg", "dog-leg", "dl" }
		vecY = __findCandStep_pow( matH, vecG, matC, matB, bTrgt, prm );
		%vecYLev = __findCandStep_lev( matH, vecG, matC, matB, bTrgt, prm );
		%norm(vecY-vecYLev)
		%deltaEtaPow = vecG'*vecY + (vecY'*matH*vecY)/2.0
		%deltaEtaLev = vecG'*vecYLev + (vecYLev'*matH*vecYLev)/2.0
	case { "g", "grad", "gradient", "gradient descent", "gradient-descent", "gradescent" }
		error( "Gradient curve is not (yet?) supported." );
	otherwise
		error( "Unsupported value of curveType." );
	endswitch
	%
	return;
endfunction


function vecY = __findCandStep_lev( matH, vecG, matC, matB, bTrgt, prm )
	hScale = norm(diag(matH));
	cScale = norm(diag(matC));
	[ matR, cholFlag ] = chol( matH );
	newtNeedsExtrap = ( 0 ~= cholFlag || min(diag(matR)) < prm.cholRelTol * max(abs(diag(matR))) );
	if ( newtNeedsExtrap )
		matR1 = chol( matH + matC * (1.0*prm.epsRegu*hScale/cScale) );
		matR2 = chol( matH + matC * (2.0*prm.epsRegu*hScale/cScale) );
		vecY1 = matR1 \ (matR1'\(-vecG));
		vecY2 = matR2 \ (matR2'\(-vecG));
		vecY_newt = (2.0*vecY1) - vecY2;
		vecYPrime1 = -( matR1 \ (matR1'\(matC*vecY1)) );
		vecYPrime2 = -( matR2 \ (matR2'\(matC*vecY1)) );
		vecYPime_newt = (2.0*vecYPrime1) - vecYPrime2;
	else
		vecY_newt = matR \ (matR'\(-vecG));
		vecYPrime_newt = -( matR \ (matR'\(matC*vecY_newt)) );
	endif
	if ( isempty(bTrgt) )
		vecY = vecY_newt;
		return;
	endif
	vecBeta_newt = matB*vecY_newt;
	b_newt = norm(vecBeta_newt);
	if ( b_newt <= bTrgt * (1.0+prm.candStepRelTol) )
		vecY = vecY_newt;
		return;
	endif
	bPrime_newt = (matB*vecY_newt)'*(matB*vecYPrime_newt)/b_newt;
	%
	pt0.mu = 0.0;
	pt0.vecY = vecY_newt;
	pt0.b = b_newt;
	pt0.bPrime = bPrime_newt;
	havePt1 = false;
	%
	levIterCount = 0;
	while (~havePt1)
		levIterCount++;
		assert( levIterCount <= prm.levIterMax );
		%
		mu0 = pt0.mu;
		b0 = pt0.b;
		bPrime0 = pt0.bPrime;
		% Intentionally overshoot in mu.
		mu = mu0 + 10.0 * ( (b0/bTrgt) - 1.0 ) * b0 / (-bPrime0);
		%
		assert( mu0 < mu );
		pt = __calcLevPtFromMu( matH, vecG, matC, matB, mu );
		assert( pt.b < pt0.b );
		%
		if ( abs(pt.b-bTrgt) < prm.candStepRelTol * bTrgt )
			vecY = pt.vecY;
			return;
		elseif ( pt.b > bTrgt )
			pt0 = pt;
		elseif ( pt.b < bTrgt )
			havePt1 = true;
			pt1 = pt;
		endif
	endwhile
	%
	applyConstraints = false;
	while (1)
		levIterCount++;
		assert( levIterCount <= prm.levIterMax );
		%
		mu0 = pt0.mu;
		b0 = pt0.b;
		bPrime0 = pt0.bPrime;
		mu1 = pt1.mu;
		b1 = pt1.b;
		bPrime1 = pt1.bPrime;
		mu_from0 = mu0 + ( (b0/bTrgt) - 1.0 ) * b0 / (-bPrime0);
		mu_from1 = mu1 + ( (b1/bTrgt) - 1.0 ) * b1 / (-bPrime1);
		mu_fromX = mu0 + ( (b0/bTrgt) - 1.0 ) * ( mu1 - mu0 ) / ( (b0/b1) - 1.0 );
		if ( abs(b0-bTrgt) < abs(b1-bTrgt) && mu_from0 < mu1 )
			mu = mu_from0;
		elseif ( abs(b0-bTrgt) > abs(b1-bTrgt) && mu_from1 > mu0 )
			mu = mu_from1;
		else
			mu = mu_fromX;
		endif
		if ( applyConstraints )
			mu = median([ mu0+0.1*(mu1-mu0), mu, mu1-0.1*(mu1-mu0) ]);
		endif
		%
		assert( mu0 < mu && mu < mu1 );
		pt = __calcLevPtFromMu( matH, vecG, matC, matB, mu );
		assert( pt1.b < pt.b && pt.b < pt0.b );
		%
		%
		if ( ~applyConstraints && abs(pt.b-bTrgt) > 0.5 * min([ abs(pt0.b-bTrgt), abs(pt1.b-bTrgt) ]) )
			applyConstraints = true;
		else
			applyConstraints = false;
		endif
		if ( abs(pt.b-bTrgt) < prm.candStepRelTol * bTrgt )
			vecY = pt.vecY;
			return;
		elseif ( pt.b > bTrgt )
			pt0 = pt;
		elseif ( pt.b < bTrgt )
			pt1 = pt;
		endif
	endwhile
	%
	return;
endfunction
function pt = __calcLevPtFromMu( matH, vecG, matC, matB, mu );
	pt.mu = mu;
	matR = chol( matH + mu*matC );
	pt.vecY = -( matR \ (matR'\vecG) );
	vecBeta = matB*pt.vecY;
	pt.vecYPrime = -( matR \ (matR'\(matC*pt.vecY)) );
	pt.b = norm(vecBeta);
	pt.bPrime = (matB*pt.vecY)'*(matB*pt.vecYPrime)/pt.b;
	return;
endfunction


function vecY = __findCandStep_pow( matH, vecG, matC, matB, bTrgt, prm )
	hScale = norm(diag(matH));
	cScale = norm(diag(matC));
	[ matR, cholFlag ] = chol( matH );
	newtNeedsExtrap = ( 0 ~= cholFlag || min(diag(matR)) < prm.cholRelTol * max(abs(diag(matR))) );
	if ( newtNeedsExtrap )
		matR1 = chol( matH + matC * (1.0*prm.epsRegu*hScale/cScale) );
		matR2 = chol( matH + matC * (2.0*prm.epsRegu*hScale/cScale) );
		vecY1 = matR1 \ (matR1'\(-vecG));
		vecY2 = matR2 \ (matR2'\(-vecG));
		vecY_newt = (2.0*vecY1) - vecY2;
		vecYPrime1 = -( matR1 \ (matR1'\(matC*vecY1)) );
		vecYPrime2 = -( matR2 \ (matR2'\(matC*vecY1)) );
		vecYPime_newt = (2.0*vecYPrime1) - vecYPrime2;
	else
		vecY_newt = matR \ (matR'\(-vecG));
		vecYPrime_newt = -( matR \ (matR'\(matC*vecY_newt)) );
	endif
	if ( isempty(bTrgt) )
		vecY = vecY_newt;
		return;
	endif
	vecBeta_newt = matB*vecY_newt;
	b_newt = norm(vecBeta_newt);
	if ( b_newt <= bTrgt * (1.0+prm.candStepRelTol) )
		vecY = vecY_newt;
		return;
	endif
	%
	matR = chol( matC );
	vecD = matR \ (matR'\(-vecG));
	gtd = vecG'*vecD;
	dthd = vecD'*matH*vecD;
	assert( dthd > 0.0 );
	assert( gtd <= 0.0 );
	s = -gtd/dthd;
	vecY_cauchy = s*vecD;
	vecBeta_cauchy = matB*vecY_cauchy;
	b_cauchy = norm(vecBeta_cauchy);
	%
	if ( b_cauchy >= bTrgt )
		vecY = vecY_cauchy * (bTrgt/b_cauchy);
		return;
	endif
	%
	% Find where the second leg intersects the boundary.
	vecY2 = vecY_newt - vecY_cauchy;
	vecBeta2 = matB*vecDelta2;
	% We can't use calcLinishRootOfQuad() here!
	% We want positive root of quad!
	%   Using y = yCauchy + t * y2,
	%   ||B*y||^2 = (t^2)*||beta2||^2 + t*(2.0*beta2'*betaCauchy) + bCauchy^2
	%   where beta2 = B*vecY2 and betaCauchy = B*vecYCauchy.
	a = sumsq(vecBeta2);
	b = 2.0*(vecBeta2'*vecBeta_cauchy);
	c = (b_cauchy^2) - (bTrgt^2);
	discrim = (b^2) - (4.0*a*c);
	assert( discrim >= 0.0 );
	assert( 0.0 < a );
	t = (-b+sqrt(discrim))/(2.0*a); % Because a must be positive.
	vecY = vecY_cauchy + (t*Y2);
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __moveTo( funchF, vecX_next, vecF_next, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __moveTo( vecX_next, vecF_next, fModelDat, prm )
	error( "To-do!" );
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __expandModel( funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	matV = fModelDat.matV; % Subspace basis matrix.
	matW_frozen = fModelDat.matW_frozen; % Only for dev/debug?
	matALo_frozen = fModelDat.matALo_frozen; % Only for dev/debug?
	matAHi_frozen = fModelDat.matAHi_frozen; % Only for dev/debug?
	matB_frozen = fModelDat.matB_frozen; % Only for dev/debug?
	%
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Only for dev/debug?
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	%
	vecRhoF = fModelDat.expand_vecU_suggested;
	vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, matV, prm );
	if ( norm(vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "__applyPrecon() failed to generate a linearly independent vector." );
		for n=1:sizeVLocal
			vecRhoF = fModelDat.matWLocal(:,1+sizeVLocal-n);
			vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
			vecV = __calcOrthonorm( vecU, matV, prm );
			if ( norm(vecV) > sqrt(eps) )
				break;
			endif
		endfor
	endif
	if ( norm(vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "__applyPrecon() still failed to generate a linearly independent vector." );
		for n=1:sizeV
			vecRhoF = fModelDat.matW(:,1+sizeV-n);
			vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
			vecV = __calcOrthonorm( vecU, matV, prm );
			if ( norm(vecV) > sqrt(eps) )
				break;
			endif
		endfor
	endif
	if ( norm(vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "__applyPrecon() yet again failed to generate a linearly independent vector." );
		for n=1:sizeX
			vecU = zeros(size(vecX));
			vecU(n) = 1.0;
			vecV = __calcOrthonorm( vecU, matV, prm );
			if ( norm(vecV) > sqrt(eps) )
				break;
			endif
		endfor
	endif
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalIncr++;
	%
	fModelDat.matV = [ matV, vecV ];
	fModelDat.matW_frozen = [ matW_frozen, vecW ];
	fModelDat.matALo_frozen = zeros(sizeV+1,sizeV+1);
	fModelDat.matALo_frozen(1:sizeV,1:sizeV) = matALo_frozen;
	fModelDat.matAHi_frozen = zeros(sizeV+1,sizeV+1);
	fModelDat.matAHi_frozen(1:sizeV,1:sizeV) = matAHi_frozen;
	fModelDat.matB_frozen = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = matB_frozen;
	fModelDat.matB_frozen(sizeV,sizeV) = norm(diag(matB))*prm.epsB; % Note: expand frozen per non-frozen.
	%
	fModelDat.matW = [ matW, vecW ];
	fModelDat.matALo = zeros(sizeV+1,sizeV+1);
	fModelDat.matALo(1:sizeV,1:sizeV) = matALo;
	fModelDat.matAHi = zeros(sizeV+1,sizeV+1);
	fModelDat.matAHi(1:sizeV,1:sizeV) = matAHi;
	fModelDat.matB = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = matB;
	fModelDat.matB(sizeV,sizeV) = norm(diag(matB))*prm.epsB;
	%
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	fModelDat.matWLocal = [ matWLocal, vecW ]; % Only for dev/debug?
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeB == sizeV ); % Not strictly required, but it's what we're doing.
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW_frozen,[sizeF,sizeV]) );
		assert( isrealarray(matALo_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matAHi_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matB_frozen,[sizeB,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( isrealarray(matWLocal,[sizeF,sizeVLocal]) );
		%
		assert( issymmetric(fModelDat.matALo) );
		assert( issymmetric(fModelDat.matAHi) );
		assert( issymmetric(fModelDat.matALo_frozen) );
		assert( issymmetric(fModelDat.matAHi_frozen) );
		% We're making matB sym too...
		assert( issymmetric(fModelDat.matB) );
		assert( issymmetric(fModelDat.matB_frozen) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(fModelDat.matV'*fModelDat.matV,eye(sizeV+1,sizeV+1)) < sqrt(eps) );
		assert( reldiff(fModelDat.matVLocal'*fModelDat.matVLocal,eye(sizeVLocal+1,sizeVLocal+1)) < sqrt(eps) );
	endif
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __refreshGivenDirection( vecY, funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	matV = fModelDat.matV; % Subspace basis matrix.
	matW_frozen = fModelDat.matW_frozen;
	matALo_frozen = fModelDat.matALo_frozen; % Only for dev/debug?
	matAHi_frozen = fModelDat.matAHi_frozen; % Only for dev/debug?
	matB_frozen = fModelDat.matB_frozen; % Only for dev/debug?
	%
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Only for dev/debug?
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( sizeVLocal < sizeV );
		assert( norm(vecY) > 0.0 );
		assert( norm(matV'*vecY) > (1.0-sqrt(eps))*norm(vecY) );
	endif
	vecY = matV*(matV'*vecY);
	%
	vecU = matV*vecY;
	vecV = __calcOrthonorm( vecU, matV_local, prm );
	if ( norm(vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__ERROR, __FILE__, __LINE__, "Given direction was clobbered by locally evaluated subspace." );
		retCode = RETCODE__INTERNAL_INCONSISTENCY;
		return;
	elseif ( norm(matVLocal'*vecV) > sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__ERROR, __FILE__, __LINE__, "Given direction is already in locally evaluated subspace." );
		retCode = RETCODE__INTERNAL_INCONSISTENCY;
		return;
	elseif ( norm(matV'*vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__ERROR, __FILE__, __LINE__, "Given is outside of allowed subspace." );
		retCode = RETCODE__INTERNAL_INCONSISTENCY;
		return;
	endif
	vecV = matV*(matV'*vecV);
	vecV /= norm(vecV);
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalIncr++;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	fModelDat.matWLocal = [ matWLocal, vecW ];
	%
	error( "To do!" );
	%
	% No update to matV nor anything _frozen.
	% But, update other quantites compared to frozen.
	%
	% Note that the vecY vectors are perhaps not perpendicular.
	% This won't work: fModelDat.matW = matW + (vecW - matW*vecY)*(vecY')/sumsq(vecY);
	fModelDat.matW = matW_frozen + ( matWLocal - matW_frozen*(matV')*matVLocal)*(matVLocal'*matV); % Yeah?
	matEVLocal = eye(sizeV,sizeV) - (matV'*matVLocal)*(matVLocal'*matV);
	matALo = matEVLocal'*matALo_frozen*matEVLocal;
	matAHi = matEVLocal'*matAHi_frozen*matEVLocal;
	fModelDat.matALo = (matALo'+matALo)/2.0;
	fModelDat.matAHi = (matAHi'+matAHi)/2.0;
	%fModelDat.matB = matB;
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeB == sizeV ); % Not strictly required, but it's what we're doing.
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( isrealarray(matWLocal,[sizeF,sizeVLocal]) );
		assert( isrealarray(matALo_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matAHi_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matB_frozen,[sizeB,sizeV]) );
		%
		assert( issymmetric(fModelDat.matALo) );
		assert( issymmetric(fModelDat.matAHi) );
		assert( issymmetric(fModelDat.matALo_frozen) );
		assert( issymmetric(fModelDat.matAHi_frozen) );
		% We're making matB sym too...
		assert( issymmetric(fModelDat.matB) );
		assert( issymmetric(fModelDat.matB_frozen) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(fModelDat.matV'*fModelDat.matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(fModelDat.matVLocal'*fModelDat.matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
	endif
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction
