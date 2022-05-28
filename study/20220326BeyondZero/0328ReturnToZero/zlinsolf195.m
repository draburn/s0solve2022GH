function [ vecX_best, vecF_best, retCode, fevalCount, stepsCount, datOut ] = zlinsolf195( funchF, vecX_initial, vecF_initial=[], prmIn=[] )
	startTime = time();
	if ( stopsignalpresent() )
		msg(__FILE__, __LINE__, "ERROR: Stop signal already present." );
		retCode = RETCODE__IMPOSED_STOP;
		break;
	endif
	%
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
			  "  time: %10.3e/%0.3e;  iter: %5d/%d;  feval: %5d/%d;  steps: %5d/%d;  size: %5d/%5d/%dx%d;  ||F||: %10.3e/%10.3e.", ...
			  time()-startTime, prm.timeMax, ...
			  iterCount, prm.iterMax, ...
			  fevalCount, prm.fevalMax, ...
			  stepsCount, prm.stepsMax, ...
			  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX,1), size(vecF,1), ...
			  norm(vecF_best), norm(vecF_initial) ) );
		endif
		[ retCode, fevalIncr, fModelDat ] = __analyzeModel( funchF, fModelDat, prm );
		fevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		
		error( "END OF VALID CODE." );
		
		msgif( prm.verbLev >= VERBLEV__ERROR, __FILE__, __LINE__, "ERROR: No acceptable action." );
		break;
	endwhile
	if ( prm.verbLev >= VERBLEV__PROGRESS )
		msg( __FILE__, __LINE__, "Final..." );
		msg( __FILE__, __LINE__, sprintf( ...
		  "  time: %10.3e/%0.3e;  iter: %5d/%d;  feval: %5d/%d;  steps: %5d/%d;  size: %5d/%5d/%dx%d;  ||F||: %10.3e/%10.3e.", ...
		  time()-startTime, prm.timeMax, ...
		  iterCount, prm.iterMax, ...
		  fevalCount, prm.fevalMax, ...
		  stepsCount, prm.stepsMax, ...
		  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX,1), size(vecF,1), ...
		  norm(vecF_best), norm(vecF_initial) ) );
	endif
return;
endfunction


function [ retCode, fevalIncr, vecF_initial, fModelDat, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn )
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
	prm.epsRelRegu = sqrt(eps);
	prm.candStepRelTol = 0.2;
	%%%prm.candStepRelTol = 1000.0*eps; %%%
	prm.findLevPrm = [];
	prm.findLevPrm.cholRelTol = prm.cholRelTol;
	prm.findLevPrm.epsRelRegu = prm.epsRelRegu;
	prm.findLevPrm.bRelRegu = prm.candStepRelTol;
	%
	prm.moveToELoCoeff = 0.9;
	prm.moveToEHiCoeff = 0.1;
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
	%
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matALo = [ 0.0 ];
	fModelDat.matAHi = [ 0.0 ];
	fModelDat.matB = [ 0.0 ];
	fModelDat.matVLocal = zeros(sizeX,0);
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
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeV <= sizeX );
		assert( sizeVLocal <= sizeV );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( issymmetric(matALo) );
		assert( issymmetric(matAHi) );
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
				matC += cScale * prm.epsRelRegu * matIV;
			endif
			cScale = norm(diag(matC));
			matRC = chol(matC);
		endif
	endif
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( isrealarray(matC,[sizeV,sizeV]) );
		assert( issymmetric(matC) );
		assert( min(diag(matC)) > 0.0 );
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
		if ( prm.candStepRelTol <= sqrt(eps) );
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
	switch ( prm.curveType )
	case { "l", "lev", "levenberg" }
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
	% Have to think about this...
	%fModelDat.expand_vecU_suggested = vecF - matW*vecY_ideal;
	%fModelDat.refresh_vecY_suggested = vecY_loVar;
	%fModelDat.tryStep_vecY_suggested = vecY_hiVar;
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
		vecY = findLevPt_0527( vecG, matH, matC, matB, bTrgt, prm.findLevPrm );
	case { "p", "powell", "dog", "dog leg", "dog-leg", "dl" }
		vecY = __findCandStep_pow( vecG, matH, matC, matB, bTrgt, prm );
	case { "g", "grad", "gradient", "gradient descent", "gradient-descent", "gradescent" }
		error( "Gradient curve is not (yet?) supported." );
	otherwise
		error( "Unsupported value of curveType." );
	endswitch
	%
	return;
endfunction


function vecY = __findCandStep_pow( vecG, matH, matC, matB, bTrgt, prm )
	hScale = norm(diag(matH));
	cScale = norm(diag(matC));
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag && min(diag(matR)) > prm.cholRelTol * max(abs(diag(matR))) )
		vecY_newt = matR \ (matR'\(-vecG));
	else
		matR1 = chol( matH + matC * (1.0*prm.epsRelRegu*hScale/cScale) );
		matR2 = chol( matH + matC * (2.0*prm.epsRelRegu*hScale/cScale) );
		vecY1 = matR1 \ (matR1'\(-vecG));
		vecY2 = matR2 \ (matR2'\(-vecG));
		vecY_newt = (2.0*vecY1) - vecY2;
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


function [ retCode, fevalIncr, fModelDat ] = __addDimensionToModel( funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	error( "Should vecRhoF be an argument? Should we pre-calc vecV in __analyzeModel()?" );
	%%%vecRhoF = fModelDat.expand_vecU_suggested;
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
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	fModelDat.matV = [ matV, vecV ];
	fModelDat.matW = [ matW, vecW ];
	fModelDat.matALo = zeros(sizeV+1,sizeV+1);
	fModelDat.matALo(1:sizeV,1:sizeV) = matALo;
	fModelDat.matAHi = zeros(sizeV+1,sizeV+1);
	fModelDat.matAHi(1:sizeV,1:sizeV) = matAHi;
	fModelDat.matB = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = matB;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeV <= sizeX );
		assert( sizeVLocal <= sizeV );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( issymmetric(matALo) );
		assert( issymmetric(matAHi) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(fModelDat.matV'*fModelDat.matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(fModelDat.matVLocal'*fModelDat.matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
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
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( isrealarray(vecY,[sizeV,1]) );
		assert( sizeVLocal < sizeV );
		assert( norm(vecY) > 0.0 );
	endif
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
		msgif( prm.verbLev >= VERBLEV__ERROR, __FILE__, __LINE__, "Given direction is somehow outside of allowed subspace." );
		retCode = RETCODE__INTERNAL_INCONSISTENCY;
		return;
	endif
	vecV = matV*(matV'*vecV); % Make sure it's fully in matV, for numerical reasons.
	vecV /= norm(vecV);
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalIncr++;
	%
	vecZ = matV'*vecV;
	matIV = eye(sizeV,sizeV);
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	%fModelDat.matV = matV;
	fModelDat.matW = matW + ( vecW - matW*vecZ )*(vecZ');
	fModelDat.matALo = ( matIV - vecZ*(vecZ') ) * matALo * ( matIV - vecZ*(vecZ') );
	fModelDat.matAHi = ( matIV - vecZ*(vecZ') ) * matAHi * ( matIV - vecZ*(vecZ') );
	%fModelDat.matB = matB;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	%
	% Explicitly symmetrize.
	fModelDat.matALo = (fModelDat.matALo'+fModelDat.matALo)/2.0;
	fModelDat.matAHi = (fModelDat.matAHi'+fModelDat.matAHi)/2.0;
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeV <= sizeX );
		assert( sizeVLocal <= sizeV );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( issymmetric(matALo) );
		assert( issymmetric(matAHi) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(fModelDat.matV'*fModelDat.matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(fModelDat.matVLocal'*fModelDat.matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
	endif
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __moveTo( vecY, vecF_next, funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	matIV = eye(sizeV,sizeV);
	%
	yNorm = norm(vecY);
	vecFModel_next = vecF + matW*vecY;
	vecRhoF = vecF_next - vecFModel_next;
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < yNorm )
		assert( norm(vecFModel_next) <= norm(vecF) );
		assert( norm(vecF_next) <= norm(vecF) );
	endif
	useQuadUpdate = true;
	if (~useQuadUpdate)
		matW_updated = matW + vecRhoF * (vecY')/(yNorm^2);
	else
		rhoSumsq = sumsq(vecRho);
		ytay = vecY'*matAHi*vecY;
		if ( rhoSumsq <= ytay )
			s = 1.0;
		else
			s = ytay/rhoSumsq;
		endif
		matW_updated = matW + (2.0-s) * vecRho * (vecY')/(yNorm^2);
	endif
	%
	vecDW = sum(matW.^2,1)';
	vecD1 = ones(sizeV,1)/sumsq(vecY);
	vecD2 = vecDW/sumsq(matW*vecY);
	vecD3  = vecDW/(veCY'*diag(vecDW)*vecY);
	vecDLo = min( [ vecD1, vecD2, vecD3 ], [], 2 );
	vecDHi = max( [ vecD1, vecD2, vecD3 ], [], 2 );
	vecYHat = vecY/norm(vecY);
	matELo = matIV - prm.moveToELoCoeff*vecYHat*(vecYHat');
	matEHi = matIV - prm.moveToEHiCoeff*vecYHat*(vecYHat');
	matALo_updated = matELo' * ( matALo + vecDLo ) * matELo;
	matAHi_updated = matEHi' * ( matAHi + vecDHi ) * matEHi;
	%
	fModelDat.vecX = vecX + matV*vecY;
	fModelDat.vecF = vecF_next;
	%fModelDat.matV = matV;
	fModelDat.matW = matW_updated;
	fModelDat.matALo = matALo_updated;
	fModelDat.matAHi = matAHi_updated;
	%fModelDat.matB = matB;
	fModelDat.matVLocal = zeros(sizeX,0);
	%
	% Explicitly symmetrize.
	fModelDat.matALo = (fModelDat.matALo'+fModelDat.matALo)/2.0;
	fModelDat.matAHi = (fModelDat.matAHi'+fModelDat.matAHi)/2.0;
	%
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeV <= sizeX );
		assert( sizeVLocal <= sizeV );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( issymmetric(matALo) );
		assert( issymmetric(matAHi) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(fModelDat.matV'*fModelDat.matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(fModelDat.matVLocal'*fModelDat.matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
	endif
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __shrinkTR( funchF, vecY, fModelDat, prm )
	[ retCode, fevalIncr, fModelDat ] = __modifyTR( funchF, vecY, fModelDat, prm );
	return;
endfunction
function [ retCode, fevalIncr, fModelDat ] = __expandTR( funchF, vecY, fModelDat, prm )
	[ retCode, fevalIncr, fModelDat ] = __modifyTR( funchF, vecY, fModelDat, prm );
	return;
endfunction
function [ retCode, fevalIncr, fModelDat ] = __modifyTR( funchF, vecY, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	matB = fModelDat.matB;
	sizeV = size(vecY,1);
	yNorm = norm(vecY);
	if ( valdLev >= VALDLEV__MEDIUM )
		sizeB = size(matB,1);
		assert( isrealarray(vecY,[sizeV,1]) );
		assert( 0.0 < yNorm );
		assert( isrealarray(matB,[sizeB,sizeV]) );
	endif
	vecYHat = vecY/yNorm;
	matEY = eye(sizeV,sizeV) - vecYHat*(vecYHat');
	matB = matEY'*matB*matEY + vecYHat*(vecYHat');
	fModelDat.matB = (matB'+matB)/2.0;
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction
