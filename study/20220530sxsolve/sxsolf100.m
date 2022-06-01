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
	datOut.fevalCountOfSteps(stepsCount+1) = fevalCount;
	datOut.fNormOfSteps(stepsCount+1) = norm(vecF);
	datOut.vecXOfSteps(:,stepsCount+1) = vecX;
	datOut.vecFOfSteps(:,stepsCount+1) = vecF;
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
			  "   time: %8.2e;  iter: %3d;  feval: %3d;  steps: %3d;  size: %3d / %3d ( / %d x %d );  omega: %8.2e.", ...
			  time()-startTime, ...
			  iterCount, ...
			  fevalCount, ...
			  stepsCount, ...
			  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX_initial,1), size(vecF_initial,1), ...
			  sumsq(vecF)/2.0 ) );
		endif
		%
		[ retCode, studyDat ] = __studyFModel( fModelDat, prm );
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		%
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __takeAction( funchF, fModelDat, studyDat, prm );
		fevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		if ( ~isempty(vecX_next) )
			if ( prm.valdLev >= VALDLEV__LOW )
				assert( ~isempty(vecF_next) );
				assert( isrealarray(vecX_next,[size(vecX,1),1]) );
				assert( isrealarray(vecF_next,[size(vecF,1),1]) );
				assert( norm(vecF_next) < norm(vecF) );
			endif
			omega_next = sumsq(vecF_next)/2.0;
			stepsCount++;
			if ( prm.verbLev >= VERBLEV__PROGRESS )
				msg( __FILE__, __LINE__, sprintf( ...
				  " Step %3d ( at t %8.2e, i %3d, f %3d with %8.2e x %8.2e ):  %8.2e -> %8.2e ( down %8.2e; tol %8.2e ).", ...
				  stepsCount, time()-startTime, iterCount, fevalCount, ...
				  norm(vecX_next-vecX) , norm(vecF_next-vecF), ...
				  omega, omega_next, omega - omega_next, prm.omegaTol ) );
			endif
			vecX = vecX_next;
			vecF = vecF_next;
			omega = omega_next;
			datOut.fevalCountOfSteps(stepsCount+1) = fevalCount;
			datOut.fNormOfSteps(stepsCount+1) = norm(vecF);
			datOut.vecXOfSteps(:,stepsCount+1) = vecX;
			datOut.vecFOfSteps(:,stepsCount+1) = vecF;
			clear vecX_next;
			clear vecF_next;
			clear omega_next;
		endif
		%
		continue;
	endwhile
	%
	if ( fevalCount > datOut.fevalCountOfSteps(end) )
		datOut.fevalCountOfSteps(stepsCount+1) = fevalCount;
		datOut.fNormOfSteps(stepsCount+1) = norm(vecF);
		datOut.vecXOfSteps(:,stepsCount+1) = vecX;
		datOut.vecFOfSteps(:,stepsCount+1) = vecF;
	endif
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
	prm.verbLev = VERBLEV__PROGRESS; prm.valdLev = VALDLEV__HIGH; % Integration dev.
	%prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Feature refinement dev.
	%prm.verbLev = VERBLEV__COPIOUS; prm.valdLev = VALDLEV__VERY_HIGH; % New feature dev.
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Refactor / debug.
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
	prm.curveType = "lev";
	prm.curveScaling = "b";
	prm.matC = [];
	prm.cholRelTol = sqrt(eps);
	prm.epsRelRegu = sqrt(eps);
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
	% Add other parameters once code more settled.
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
	fModelDat.matWLocal = [ vecW ];
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, taFevalCount, fModelDat, vecX_next, vecF_next ] = __takeAction( funchF, fModelDat, studyDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	taFevalCount = 0;
	vecX_next = [];
	vecF_next = [];
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Projection of locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	vecY_unb = studyDat.vecY_unb;
	vecY_bnd = studyDat.vecY_bnd;
	vecY_loc = studyDat.vecY_loc;
	eta_unb = studyDat.eta_unb;
	eta_bnd = studyDat.eta_bnd;
	eta_loc = studyDat.eta_loc;
	%
	%
	%
	omega = sumsq(vecF)/2.0;
	omega_initial = sumsq(vecF)/2.0;
	%
	omegaThresh = 0.1 * omega * ( omega / omega_initial )^0.5;
	if ( eta_unb > max([ 0.1*omegaThresh, prm.omegaTol ]) )
		vecR = matW(:,end);
		vecS = vecR;
		vecU = __applyPrecon( vecS, prm, vecX, vecF );
		vecV = __calcOrthonorm( vecU, matV, prm );
		if ( norm(vecV) > sqrt(eps) )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Action: Expand subspace (per actual residual)." );
			[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif
		%
		% Although f - W*y is the residual, if the latest w is perp to f, this quantity will be repeated
		%  between iterations. So, we should also properly consider the Kryov subspace...
		vecR = matW(:,end);
		vecS = vecR;
		vecU = __applyPrecon( vecS, prm, vecX, vecF );
		vecV = __calcOrthonorm( vecU, matV, prm );
		if ( norm(vecV) > sqrt(eps) )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Action: Expand subspace (per Krylov)." );
			[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif
	endif
	%
	%
	if ( eta_loc <= max([ omegaThresh, prm.omegaTol ]) )
		assert( ~isempty(vecY_loc) );
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( "  Action: Try step ( %8.2e -> %8.2e ).", omega, eta_loc ) );
		vecY = matV'*(matVLocal*vecY_loc);
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	vecV = __calcOrthonorm( matV*vecY_bnd, matVLocal, prm );
	if ( norm(vecV) > sqrt(eps) )
		vecV = matV*(matV'*vecV); % Force in subspace for numerical stability..
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Action: Re-evaluate direction." );
		[ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	% Desperation.
	% See if we can find any w outside of matV.
	for n=1:sizeV
		vecR = matW(:,n);
		vecS = vecR;
		vecU = __applyPrecon( vecS, prm, vecX, vecF );
		vecV = __calcOrthonorm( vecU, matV, prm );
		if ( norm(vecV) > sqrt(eps) )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Action: Expand subspace (desperation)." );
			[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif
	endfor
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "ACTION" FUNCTIONS
%

function [ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Projection of locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( sizeV < sizeX );
		assert( isrealarray(vecV,[sizeX,1]) );
		assert( abs(norm(vecV)-1.0) < sqrt(eps) );
		assert( norm(matV'*vecV) < sqrt(eps) );
		assert( norm(matVLocal'*vecV) < sqrt(eps) );
	endif
	%
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalIncr++;
	%
	%
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	fModelDat.matV = [ matV, vecV ];
	fModelDat.matW = [ matW, vecW ];
	fModelDat.matB = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = matB;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	fModelDat.matWLocal = [ matWLocal, vecW ];
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, tsFevalCount, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	tsFevalCount = 0;
	vecX_next = [];
	vecF_next = [];
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Projection of locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	vecY_unb = studyDat.vecY_unb;
	vecY_bnd = studyDat.vecY_bnd;
	vecY_loc = studyDat.vecY_loc;
	eta_unb = studyDat.eta_unb;
	eta_bnd = studyDat.eta_bnd;
	eta_loc = studyDat.eta_loc;
	%
	assert( norm(vecY) > 0.0 );
	%
	vecX_trial = vecX + matV * vecY;
	vecF_trial = funchF( vecX_trial );
	tsFevalCount++;
	%
	vecFModel = vecF + matW * vecY;
	vecRho = vecF_trial - vecFModel;
	omega = sumsq(vecF)/2.0;
	omega_trial = sumsq(vecF_trial)/2.0;
	%
	relFallThresh = 0.5;  %%% Make param.
	omegaThresh = omega - relFallThresh * ( omega - eta_bnd ); % Require fall to be at least half of bound ideal.
	if ( omega_trial <= omegaThresh )
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		  "  Accepting step: %8.2e -> %8.2e, leq %8.2e ( down frac %8.2e, remain frac %8.2e ).", ...
		  omega, omega_trial, omegaThresh, 1.0 - omega_trial/omega, omega_trial/omega ) );
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		%
		trExpandCoeff = 1.5; %%% Make param.
		if ( trExpandCoeff*norm(matB*vecY) >= 1.0 )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Expanding trust region." );
			[ retCode, fModelDat ] = __expandTR( trExpandCoeff*vecY, fModelDat, prm );
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
		endif
		%
		[ retCode,  fModelDat ] = __moveTo( vecY, vecF_next, fModelDat, prm );
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		%
		%
		%
		useCoasting = true; %%% Make param.
		if (useCoasting)
			p = 0.0;
			%
			coastMax = 5; %%% Make param.
			for coastCount=1:coastMax
				%
				% Keep vecY for now, but overwrite.
				% Keep fModelDat, but overwrite.
				clear vecX;
				clear vecF;
				clear matV;
				clear matW;
				clear matB;
				clear matVLocal;
				clear matWLocal;
				clear sizeX;
				clear sizeF;
				clear sizeV;
				clear sizeB;
				clear sizeVLocal;
				%
				clear studyDat;
				clear vecY_unb;
				clear vecY_bnd;
				clear vecY_loc;
				clear eta_unb;
				clear eta_bnd;
				clear eta_loc;
				%
				clear vecX_trial;
				clear vecF_trial;
				clear vecFModel;
				% Keep vecRho, but overwrite.
				clear omega;
				clear omega_trial;
				%
				% Keep vecX_next, but overwrite.
				% Keep vecF_next, but overwrite.
				%
				p += sumsq( vecRho ) / sumsq( vecY );
				vecX = vecX_next;
				vecF = vecF_next;
				omega = sumsq(vecF)/2.0;
				[ retCode, studyDat ] = __studyFModel( fModelDat, prm );
				if ( 0~= retCode )
					msgretcodeif( true, __FILE__, __LINE__, retCode );
					return;
				endif
				matV = fModelDat.matV;
				matW = fModelDat.matW;
				vecY = studyDat.vecY_bnd;
				eta = studyDat.eta_bnd;
				eta_coast = eta + p*sumsq( vecY )/2.0;
				coastingFallRelThresh = 0.1; %%% Make param.
				coastingC1 = 0.5;
				if ( eta_coast < (1.0-coastingFallRelThresh)*omega && eta < (1.0-coastingC1)*omega )
					vecX_trial = vecX + matV*vecY;
					vecF_trial = funchF( vecX_trial );
					tsFevalCount++;
					%
					vecFModel = vecF + matW * vecY;
					vecRho = vecF_trial - vecFModel;
					omega_trial = sumsq(vecF_trial)/2.0;
					%
					if ( omega_trial < omega )
						msgif( prm.verbLev >= VERBLEV__DETAILED, __FILE__, __LINE__, sprintf( ...
						  "   Coasting was good ( %0.3e -> %0.3e ~ %0.3e? -> %0.3e ).", omega, eta, eta_coast, omega_trial ) );
						vecX_next = vecX_trial;
						vecF_next = vecF_trial;
						[ retCode, fModelDat ] = __moveTo( vecY, vecF_next, fModelDat, prm );
						if ( 0~= retCode )
							msgretcodeif( true, __FILE__, __LINE__, retCode );
							return;
						endif
						continue;
					else
						msgif( prm.verbLev >= VERBLEV__DETAILED, __FILE__, __LINE__, sprintf( ...
						  "   Coasting was bad ( %0.3e -> %0.3e ~ %0.3e? -> %0.3e ).", omega, eta, eta_coast, omega_trial ) );
						break;
					endif
				else
					msgif( prm.verbLev >= VERBLEV__DETAILED, __FILE__, __LINE__, sprintf( ...
					  "   Not attempting to coast ( %0.3e -> %0.3e ~ %0.3e? ).", omega, eta, eta_coast ) );
					break;
				endif
			endfor
		endif
		%
		%
		retCode = RETCODE__SUCCESS;
		return;
	endif
	%
	%
	if ( omega_trial < omega )
		error( "Omega decreased, but not as much as desired; we need to keep track of this (x,F), in case we can't do better!" );
	endif
	%
	msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
	  "  Rejecting step: %8.2e -> %8.2e gt %8.2e ( up frac %8.2e, scale frac %8.2e ).", ...
	  omega, omega_trial, omegaThresh, omega_trial/omega - 1.0, omega_trial/omega ) );
	%
	msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Shrinking trust region." );
	trShrinkCoeff = 0.5; %%% Make param.
	[ retCode, fModelDat ] = __shrinkTR( trShrinkCoeff*vecY, fModelDat, prm );
	if ( 0~= retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Projection of locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( sizeVLocal < sizeV );
		assert( isrealarray(vecV,[sizeX,1]) );
		assert( abs(norm(vecV)-1.0) < sqrt(eps) );
		assert( abs(norm(matV'*vecV)-1.0) < sqrt(eps) );
		assert( norm(matVLocal'*vecV) < sqrt(eps) );
	endif
	% Do an additional normalization to be safe...
	vecV = matV*(matV'*vecV);
	assert( norm(vecV)>0.0 );
	vecV /= norm(vecV);
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalIncr++;
	%
	vecY = matV'*vecV;
	assert( norm(vecY)>0.0 );
	vecYHat = vecY/norm(vecY);
	matW_updated = matW + ( vecW - matW*vecYHat )*(vecYHat');
	%
	%
	%
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	%fModelDat.matV = matV;
	fModelDat.matW = matW_updated;
	%fModelDat.matB = matB;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	fModelDat.matWLocal = [ matWLocal, vecW ];
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RANDOM ACCESS FUNCTIONS
%

function [ retCode, studyDat ] = __studyFModel( fModelDat, prm )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__NOT_SET;
	studyDat = [];
	%
	%
	%vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal;  % Projection of locally evaluated subspace basis matrix.
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
			msgif( prm.verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, "  Curve scaling matrix was zero; setting to I." );
			matC = matIV;
		else
			[ matRC, cholFlag ] = chol(matC);
			if ( 0 ~= cholFlag || min(diag(matRC)) < prm.cholRelTol * max(abs(diag(matRC))) )
				msgif( prm.verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "  Curve scaling matrix was non positive-definite; applying regularization." );
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
	eta_unb = sumsq( vecF + matW*vecY_unb )/2.0;
	eta_bnd = sumsq( vecF + matW*vecY_bnd )/2.0;
	b_unb = norm(matB*vecY_unb);
	b_bnd = norm(matB*vecY_bnd);
	if ( isempty(matVLocal) )
		vecY_loc = [];
		eta_loc = omega;
		b_loc = 0.0;
	else
		matWTWLocal = matWLocal'*matWLocal;
		vecWTFLocal = matWLocal'*vecF;
		matBLocal = matB*(matV'*matVLocal);
		matCLocal = (matVLocal'*matV)*matC*(matV'*matVLocal);
		matCLocal = (matCLocal'+matCLocal);
		vecY_loc = findLevPt_0527( vecWTFLocal, matWTWLocal, matCLocal, matBLocal, 1.0, prm.findLevPrm );
		eta_loc = max([ 0.0, omega + vecWTFLocal'*vecY_loc + abs((vecY_loc'*matWTWLocal*vecY_loc)/2.0) ]);
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
		%msg( __FILE__, __LINE__, sprintf( ...
		%  "   omega: %8.2e / %8.2e / %8.2e / %8.2e (/ %8.2e);  y: %8.2e / %8.2e / %8.2e;  b: %8.2e / %8.2e / %8.2e.", ...
		%  omega, eta_loc, eta_bnd, eta_unb, prm.omegaTol, ...
		%  norm(vecY_loc), norm(vecY_bnd), norm(vecY_unb), ...
		%  b_loc, b_bnd, b_unb ) );
		msg( __FILE__, __LINE__, sprintf( ...
		  "   eta: %8.2e / %8.2e / %8.2e / %8.2e (/ %8.2e);  y: %8.2e / %8.2e / %8.2e;  b: %8.2e / %8.2e / %8.2e.", ...
		  eta_unb, eta_bnd, eta_loc, omega, prm.omegaTol, ...
		  norm(vecY_unb), norm(vecY_bnd), norm(vecY_loc), ...
		  b_bnd, b_unb, b_loc ) );
	endif
	%
	%if ( prm.valdLev >= VALDLEV__LOW )
	%	__validateStudyDat( fModelDat, studyDat, prm );
	%endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


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


function [ retCode, fModelDat ] = __moveTo( vecY, vecF_next, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	%
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Projection of locally evaluated subspace basis matrix.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	%
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
	if (useQuadUpdate)
		matW_updated = matW + 2.0 * vecRhoF * (vecY')/(yNorm^2);
	else
		matW_updated = matW + vecRhoF * (vecY')/(yNorm^2);
	endif
	%
	%
	%
	fModelDat.vecX = vecX + matV*vecY;
	fModelDat.vecF = vecF_next;
	%fModelDat.matV = matV;
	fModelDat.matW = matW_updated;
	%fModelDat.matB = matB;
	fModelDat.matVLocal = zeros(sizeX,0);
	fModelDat.matWLocal = zeros(sizeF,0);
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
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


function [ retCode, fModelDat ] = __shrinkTR( vecY, fModelDat, prm )
	[ retCode, fModelDat ] = __modifyB( vecY, fModelDat, prm );
	return;
endfunction
function [ retCode, fModelDat ] = __expandTR( vecY, fModelDat, prm )
	[ retCode, fModelDat ] = __modifyB( vecY, fModelDat, prm );
	return;
endfunction
function [ retCode, fModelDat ] = __modifyB( vecY, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	%
	matB = fModelDat.matB;
	sizeV = size(vecY,1);
	yNorm = norm(vecY);
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		sizeB = size(matB,1);
		assert( isrealarray(vecY,[sizeV,1]) );
		assert( 0.0 < yNorm );
		assert( isrealarray(matB,[sizeB,sizeV]) );
	endif
	vecYHat = vecY/yNorm;
	matEY = eye(sizeV,sizeV) - vecYHat*(vecYHat')/yNorm;
	matB = matEY'*matB*matEY + vecYHat*(vecYHat')/yNorm;
	fModelDat.matB = (matB'+matB)/2.0;
	%
	retCode = RETCODE__SUCCESS;
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
	matWLocal = fModelDat.matWLocal;  % Projection of locally evaluated subspace basis matrix.
	%
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
		assert( isrealarray(matWLocal,[sizeF,sizeVLocal]) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
	if ( ~isempty(matVLocal) )
		assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(matVLocal'*matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
		assert( reldiff(matV*(matV'*matVLocal),matVLocal) < sqrt(eps) );
		assert( reldiff(matWLocal,matW*matV'*matVLocal) < sqrt(eps) );
	endif
	endif
	%
	%
	return;
endfunction
