function [ vecX, vecF, retCode, fevalCount, stepsCount, datOut ] = sxsolf100( funchF, vecX_initial, vecF_initial=[], prmIn=[] )
	mydefs;
	startTime = time();
	vecX = [];
	vecF = [];
	retCode = RETCODE__NOT_SET;
	fevalCount = 0;
	stepsCount = 0;
	datOut = [];
	if ( stopsignalpresent() )
		msg(__FILE__, __LINE__, "ERROR: Stop signal already present." );
		retCode = RETCODE__IMPOSED_STOP;
		return;
	endif
	%
	[ retCode, fevalIncr, vecF_initial, fModelDat, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn );
	fevalCount += fevalIncr; clear fevalIncr;
	if ( 0 ~= retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	%
	iterCount = 0;
	vecX = vecX_initial;
	vecF = vecF_initial;
	omega = sumsq(vecF)/2.0;
	while (1)
		if ( norm(vecF) <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: norm(vecF) <= prm.fTol." );
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
		elseif ( isempty(fModelDat) )
			[ retCode, fevalIncr, fModelDat ] = __initFModel( funchF, vecX, vecF, prm );
			fevalCount += fevalIncr; clear fevalIncr;
			if ( 0~=retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
		endif
		%
		%
		iterCount++;
		if ( prm.verbLev >= VERBLEV__DETAILS )
			msg( __FILE__, __LINE__, sprintf( ...
			  "   time: %9.2e;  iter: %3d;  feval: %3d;  steps: %3d;  size: %3d / %3d ( / %d x %d );  omega: %8.2e.", ...
			  time()-startTime, ...
			  iterCount, ...
			  fevalCount, ...
			  stepsCount, ...
			  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX_initial,1), size(vecF_initial,1), ...
			  sumsq(vecF)/2.0 ) );
		endif
		[ retCode, fevalIncr, studyDat ] = __studyFModel( funchF, fModelDat, prm );
		fevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		%
		%
	endwhile
	%
	error( "END OF VALD CODE." );
return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAJOR SEQUENTIAL FUNCTIONS
%

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
	assert( isrealarray(vecX_initial,[sizeX,1]) );
	assert( isrealarray(vecF_initial,[sizeF,1]) );
	%
	%prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__LOW; % "Production / optimization".
	%prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration testing.
	%prm.verbLev = VERBLEV__PROGRESS; prm.valdLev = VALDLEV__HIGH; % Integration dev.
	%prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Feature refinement dev.
	%prm.verbLev = VERBLEV__COPIOUS; prm.valdLev = VALDLEV__VERY_HIGH; % New feature dev.
	prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Refactor / debug.
	%
	prm.timeMax = -1.0;
	prm.iterMax = ceil( 100 + 10*sqrt(sizeX+sizeF) + sizeX );
	prm.fevalMax = prm.iterMax;
	prm.stepsMax = 100;
	prm.fTol = sizeF*100.0*eps;
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
	%%%prm.candStepRelTol = 0.2;
	%%%prm.candStepRelTol = sqrt(eps);
	prm.candStepRelTol = 1.0e-4;
	prm.findLevPrm = [];
	prm.findLevPrm.cholRelTol = prm.cholRelTol;
	prm.findLevPrm.epsRelRegu = prm.epsRelRegu;
	prm.findLevPrm.bRelTol = prm.candStepRelTol;
	prm.findLevPrm.verbLev = VERBLEV__WARNING;
	prm.findLevPrm.valdLev = prm.valdLev;
	%
	prm.fModelDat_initial = [];
	%
	prm.precon_funchPrecon = [];
	prm.precon_matL = [];
	prm.precon_matU = [];
	prm.precon_matJA = [];
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
	fModelDat = prm.fModelDat_initial;
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validatePrm( prm );
	endif
	if ( prm.verbLev >= VERBLEV__INFO )
		msg( __FILE__, __LINE__, sprintf( ...
		  " Limits = { time: %0.2e;  iter: %d;  feval: %d;  steps: %d;  size: %d x %d;  omega: %0.2e }.", ...
		  prm.timeMax, ...
		  prm.iterMax, ...
		  prm.fevalMax, ...
		  prm.stepsMax, ...
		  size(vecX_initial,1), size(vecF_initial,1), ...
		  prm.omegaTol ) );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __initFModel( funchF, vecX, vecF, prm )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		sizeX = size(vecX,1);
		sizeF = size(vecF,1);
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
	endif
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	fModelDat = [];
	%
	vecU = vecF;
	vecV = __calcOrthonorm( vecU, [], prm );
	if ( norm(vecV) < 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "__calcOrthonorm() ." );
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
	fModelDat.vecX_initial = vecX;
	fModelDat.vecF_initial = vecF;
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matB = [ 0.0 ];
	fModelDat.matVLocal = [ vecV ];
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, studyDat ] = __studyFModel( funchF, fModelDat, prm )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	studyDat = [];
	%
	%
	%vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	omega = sumsq(vecF)/2.0;
	matIV = eye(sizeV,sizeV);
	matWTW = matW'*matW;
	matD = diag(diag(matWTW));
	vecWTF = matW'*vecF;
	matBTB = matB'*matB;
	%
	if ( ~isempty(prm.matC) )
		matC = matV'*prm.matC*matV;
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
			matC = matWTW;
		case { "m", "marq", "marquardt", "ddwtw" }
			matC = matD;
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
	%
	vecY_unb = findLevPt_0527( vecWTF, matWTW, matC, [], [], prm.findLevPrm );
	vecY_bnd = findLevPt_0527( vecWTF, matWTW, matC, matB, 1.0, prm.findLevPrm );
	eta_unb = max([ 0.0, omega + vecWTF'*vecY_unb + abs((vecY_unb'*matWTW*vecY_unb)/2.0) ]);
	eta_bnd = max([ 0.0, omega + vecWTF'*vecY_bnd + abs((vecY_bnd'*matWTW*vecY_bnd)/2.0) ]);
	b_unb = norm(matB*vecY_unb);
	b_bnd = norm(matB*vecY_bnd);
	if ( isempty(matVLocal) )
		vecY_loc = [];
		eta_loc = omega;
		b_loc = 0.0;
	else
		matWLocal = matW*(matV'*matVLocal);
		matWTWLocal = matWLocal'*matWLocal;
		vecWTFLocal = matWLocal'*vecF;
		matCLocal = matVLocal'*matV*matC*matV'*matVLocal;
		matBLocal = matB*(matV'*matVLocal);
		vecY_loc = findLevPt_0527( vecWTFLocal, matWTWLocal, matCLocal, matBLocal, 1.0, prm.findLevPrm );
		eta_loc = max([ 0.0, omega + vecWTF'*vecY_loc + abs((vecY_loc'*matWTW*vecY_loc)/2.0) ]);
		b_loc = norm(matBLocal*vecY_loc);
	endif
	%
	%
	studyDat.vecY_unb = vecY_unb;
	studyDat.vecY_bnd = vecY_bnd;
	studyDat.vecY_loc = vecY_loc;
	studyDat.eta_unb = eta_unb;
	studyDat.eta_bnd = eta_bnd;
	studyDat.eta_loc = eta_loc;
	%
	if ( prm.verbLev >= VERBLEV__DETAILS )
		msg( __FILE__, __LINE__, sprintf( ...
		  "   size: %3d / %3d ( / %d x %d );  omega: %8.2e / %8.2e / %8.2e / %8.2e (/ %8.2e);  b: %8.2e / %8.2e / %8.2e.", ...
		  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), ...
		  size(fModelDat.vecX_initial,1), size(fModelDat.vecF_initial,1), ...
		  omega, eta_loc, eta_bnd, eta_unb, prm.omegaTol, ...
		  b_loc, b_bnd, b_unb ) );
	endif
	%
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RANDOM ACCESS FUNCTIONS
%

function vecU = __applyPrecon( vecRhoF, prm, vecX, vecF )
	% TODO: Oh boy...
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
	mydefs;
	u0 = norm(vecU);
	if (0.0==u0)
		vecV = zeros(size(vecU));
		return;
	elseif (isempty(matV))
		vecV = vecU/u0;
		return;
	elseif ( size(matV,2) >= size(matV,1) )
		vecV = zeros(size(vecU));
		return;
	endif
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( reldiff(matV'*matV,eye(size(matV,2),size(matV,2))) < sqrt(eps) );
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
	% TODO: Update epsFD based on steps taken.... yeah, right.
	mydefs;
	v = norm(vecV);
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < v );
	endif
	vecFP = funchF( vecX + prm.epsFD*vecV );
	vecW = ( vecFP - vecF ) / (prm.epsFD * v);
	return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEV FUNCTIONS
%

function __validatePrm( prm )
	mydefs;
	if ( prm.valdLev < VALDLEV__LOW )
		return;
	endif
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( isrealscalar(prm.verbLev) );
		assert( isrealscalar(prm.valdLev) );
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
	return;
endfunction


function __validateFModelDat( fModelDat, prm )
	mydefs;
	if ( prm.valdLev < VALDLEV__LOW )
		return;
	endif
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	sizeX = size(fModelDat.vecX,1);
	sizeF = size(fModelDat.vecF,1);
	sizeV = size(fModelDat.matV,2);
	sizeB = size(fModelDat.matB,1);
	sizeVLocal = size(fModelDat.matVLocal,2);
	%
	%
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		assert( sizeV <= sizeX );
		assert( sizeVLocal <= sizeV );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
	if ( ~isempty(matVLocal) )
		assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(matVLocal'*matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
		assert( reldiff(matV*(matV'*matVLocal),matVLocal) < sqrt(eps) );
	endif
	endif
	%
	%
	return;
endfunction
