% sxsolf100, but looking at revised actions, per ideas durign precon integration.
% sxsolf180, but with hacks removed(?)
% sxsolf181, but more hacks to improve actions (ish).
% sxsolf182 revised -- not so hacky, but keeping history.

function [ vecX, vecF, retCode, fevalCount, stepsCount, datOut ] = sxsolf183( funchF, vecX_initial, vecF_initial=[], prmIn=[] )
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
	vecX = vecX_initial; % These represent the best seen so far;
	vecF = vecF_initial; %  our current "location" is per the values in fModelDat.
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
			[ retCode, fevalIncr, fModelDat ] = __initFModel( funchF, vecX_initial, vecF_initial, prm );
			fevalCount += fevalIncr; clear fevalIncr;
			if ( 0~=retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
		endif
		%
		%
		iterCount++;
		if ( prm.verbLev >= VERBLEV__COPIOUS )
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
		
		fModelDat.iterCount = iterCount;
		fModelDat.stepsCount = stepsCount; % Not the same as fModelDat's localeIndex; "stepsCount" measures changes to "best".
		fModelDat.fevalCount = fevalCount;
		
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __takeAction_hist( funchF, fModelDat, studyDat, prm );
		fevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		if ( ~isempty(vecX_next) )
			if ( prm.valdLev >= VALDLEV__LOW )
				assert( ~isempty(vecF_next) );
				assert( isrealarray(vecX_next,[size(vecX_initial,1),1]) );
				assert( isrealarray(vecF_next,[size(vecF_initial,1),1]) );
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
	datOut.fModelDat = fModelDat;
	if ( fevalCount > datOut.fevalCountOfSteps(end) )
		stepsCount++;
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
	%prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__LOW; % Post-establishment.
	%prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__ZERO; % Performance testing.
	%prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__LOW; % Routine use.
	%prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration testing.
	%prm.verbLev = VERBLEV__PROGRESS+10; prm.valdLev = VALDLEV__VERY_HIGH; % Integration dev.
	prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__VERY_HIGH; % Feature refinement dev.
	%prm.verbLev = VERBLEV__COPIOUS; prm.valdLev = VALDLEV__VERY_HIGH; % New feature dev.
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Refactor / debug.
	%
	prm.timeMax = -1.0;
	prm.iterMax = ceil( 100 + 10*sqrt(sizeX+sizeF) + sizeX );
	prm.fevalMax = prm.iterMax;
	prm.stepsMax = 100;
	prm.fTol = sqrt(sizeF)*100.0*eps;
	%
	prm.epsFD = 1.0e-4;
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
	prm.findLevPrm.valdLev = VALDLEV__ZERO;
	%
	prm.fModelDat_initial = [];
	%
	prm.precon_funchPrecon = [];
	prm.precon_matL = [];
	prm.precon_matU = [];
	prm.precon_matJA = [];
	%
	% Add other parameters once code more settled.
	prm.useDogLeg = false;
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
	fModelDat.vecX = vecX; % Current "local" point.
	fModelDat.vecF = vecF;
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matB = [ 0.0 ];
	fModelDat.matVLocal = [ vecV ];
	fModelDat.matWLocal = [ vecW ];
	fModelDat.vecX_cand = []; % Candidate for next "local" point.
	fModelDat.vecF_cand = [];
	fModelDat.vecX_prev = [];
	fModelDat.vecF_prev = [];
	%
	
	% HISTDAT
	fModelDat.localeCount = 1;
	
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, taFevalCount, fModelDat, vecX_next, vecF_next ] = __takeAction_hist( funchF, fModelDat, studyDat, prm )
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
	omega = sumsq(vecF)/2.0;
	omega_initial = sumsq(fModelDat.vecF_initial)/2.0;
	if (isempty(fModelDat.vecF_prev))
		omega_prev = 0.0;
		omegaTolTemp = max([ 0.7*prm.omegaTol, 0.01*omega ]);
	else
		omega_prev = sumsq(fModelDat.vecF_prev)/2.0;
		omegaTolTemp = max([ 0.7*prm.omegaTol, 0.01*omega*min([ 1.0, omega/omega_initial ]) ]);
	endif
	%
	vecY_unb = studyDat.vecY_unb;
	vecY_bnd = studyDat.vecY_bnd;
	vecY_loc = studyDat.vecY_loc;
	vecY_locunb = studyDat.vecY_locunb;
	eta_unb = studyDat.eta_unb;
	eta_bnd = studyDat.eta_bnd;
	eta_loc = studyDat.eta_loc;
	eta_locunb = studyDat.eta_locunb;
	vecRho_locunb = vecF + matWLocal*vecY_locunb;
	vecRho_loc = vecF + matWLocal*vecY_loc;
	vecRho_bnd = vecF + matW*vecY_bnd;
	vecRho_unb = vecF + matW*vecY_unb;
	
	% HISTDAT...
	ACTION_TYPE__TRY_STEP = 100;
	ACTION_TYPE__PULL_RECORD = 200;
	ACTION_TYPE__EXPAND_SUBSPACE = 300;
	if ( 0 == sizeVLocal )
		fModelDat.localeCount++;
	endif
	iterCount = fModelDat.iterCount;
	localeCount = fModelDat.localeCount;
	
	
	% PROBE HISTDAT.
	histDatSays_weShouldCutAndRun = false;
	histDatSays_weShouldNotTrustRecord = false;
	lookBack_iterDistance = 0;
	if ( 1 <= sizeVLocal && 2 <= iterCount )
		% This algorithm is probably not quite sensible but it's a start.
		badPullCount = 0;
		for n = iterCount - 1 : -1 : 1
			if ( 0 == fModelDat.histDat(n).iterDat.sizeVLocal )
				break;
			elseif ( abs( 1.0 - fModelDat.histDat(n).iterDat.eta_loc / eta_loc ) > 0.01 )
				n++;
				break;
			endif
			if ( fModelDat.histDat(n).iterDat.actionType == ACTION_TYPE__PULL_RECORD )
				badPullCount++;
			endif
		endfor
		lookBack_iterDistance = iterCount - n;
		%
		if ( badPullCount >= 2 )
			histDatSays_weShouldNotTrustRecord = true;
			msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
			  "      History says we should not trust records (%d).", badPullCount ) );
		endif
	endif
	if ( 3 <= sizeVLocal )
		for n = iterCount-1 : -1 : 1
		if ( fModelDat.histDat(n).iterDat.sizeVLocal <= sizeVLocal - 2 )
			break;
		endif
		endfor
		slopeFromStart = ( log(eta_loc) - log(omega) ) / ( sizeVLocal - 0.0 );
		lookBack_eta_loc = fModelDat.histDat(n).iterDat.eta_loc;
		deltaK = sizeVLocal - fModelDat.histDat(n).iterDat.sizeVLocal;
		assert( deltaK > 0 );
		slopeRecent = ( log(eta_loc) - log(lookBack_eta_loc) ) / deltaK;
		if ( abs(slopeRecent) < 0.01 * abs(slopeFromStart) )
			histDatSays_weShouldCutAndRun = true;
			msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
			  "      History says we should cut and run: %0.3e = ( l(%0.3e) - l(%0.3e) ) / %d  vs   %0.3e = ( l(%0.3e) - l(%0.3e) ) / %d.", ...
			  slopeRecent, eta_loc, lookBack_eta_loc, deltaK, ...
			  slopeFromStart, eta_loc, omega, sizeVLocal ) );
		endif
	endif
	
	
	% UPDATE HISTDAT.
	fModelDat.histDat(iterCount).iterDat.localeCount = localeCount;
	fModelDat.histDat(iterCount).iterDat.sizeV = sizeV;
	fModelDat.histDat(iterCount).iterDat.sizeVLocal = sizeVLocal;
	fModelDat.histDat(iterCount).iterDat.eta_unb = eta_unb; % Record, unbounded.
	fModelDat.histDat(iterCount).iterDat.eta_bnd = eta_bnd; % Record, bounded.
	fModelDat.histDat(iterCount).iterDat.eta_loc = eta_loc; % Local, bounded.
	fModelDat.histDat(iterCount).iterDat.eta_locunb = eta_locunb; % Local, unbounded.
	fModelDat.histDat(iterCount).iterDat.omega = omega; % Local, unbounded.
	fModelDat.histDat(iterCount).iterDat.stepsCount = fModelDat.stepsCount;
	fModelDat.histDat(iterCount).iterDat.fevalCount = fModelDat.fevalCount;
	fModelDat.histDat(iterCount).iterDat.histDatSays_weShouldCutAndRun = histDatSays_weShouldCutAndRun;
	fModelDat.histDat(iterCount).iterDat.histDatSays_weShouldNotTrustRecord = histDatSays_weShouldNotTrustRecord;
	fModelDat.histDat(iterCount).iterDat.lookBack_iterDistance = lookBack_iterDistance;
	if ( 0 == sizeVLocal )
		fModelDat.localeDat(localeCount).iterCount = iterCount;
		fModelDat.localeDat(localeCount).vecX = vecX;
		fModelDat.localeDat(localeCount).vecF = vecF;
	endif
	
	
	if (  1 <= sizeVLocal  &&  ( eta_loc < omegaTolTemp || histDatSays_weShouldCutAndRun || sizeVLocal == sizeX )  )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Try step (strong cnvg)" ) );
		vecY = matV'*(matVLocal*vecY_loc);
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__TRY_STEP;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 0;
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	if (~histDatSays_weShouldNotTrustRecord)
	if (  0 == sizeVLocal || eta_bnd < 0.1 * eta_loc  ||  eta_bnd < 0.5 * omegaTolTemp  )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Pull VR*yRB" ) );
		vecV = __calcOrthonorm( matV*vecY_bnd, matVLocal, prm );
		if ( norm(vecV) > sqrt(eps) )
			vecV = matV*(matV'*vecV); % Force in subspace for numerical stability.
			fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__PULL_RECORD;
			fModelDat.histDat(iterCount).iterDat.actionSubType = 0;
			[ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "  WHOOPS, WE'VE ALREADY TRIED THIS!" );
	endif
	endif
	%
	vecU = __applyPrecon( vecRho_loc, prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, matV, prm );
	if ( norm(vecV) > sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Expand (rhoLB vs VR)" ) );
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__EXPAND_SUBSPACE;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 0;
		[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	vecU = __applyPrecon( vecRho_locunb, prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, matV, prm );
	if ( norm(vecV) > sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Expand (rhoLU vs VR)" ) );
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__EXPAND_SUBSPACE;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 10;
		[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	if (  1 <= sizeVLocal  &&  ( eta_loc < 0.5*omega || histDatSays_weShouldCutAndRun || sizeVLocal == sizeX )  )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Try step (weak cnvg)" ) );
		vecY = matV'*(matVLocal*vecY_loc);
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__TRY_STEP;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 10;
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	for n=1:sizeV
	vecU = __applyPrecon( matW(:,n), prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, matV, prm );
	if ( norm(vecV) > sqrt(eps) )
		strAction = sprintf( "Expand (wR(%d) vs VR)", n );
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), strAction ) );
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__EXPAND_SUBSPACE;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 20;
		[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	endfor
	%
	for n=1:sizeVLocal
	vecU = __applyPrecon( matWLocal(:,n), prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, matVLocal, prm );
	if ( norm(vecV) > sqrt(eps) )
		strAction = sprintf( "Expand (wL(%d) vs VL)", n );
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), strAction ) );
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__EXPAND_SUBSPACE;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 30;
		[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	endfor
	%
	for n=1:sizeX
	vecU = zeros(sizeX,1);
	vecU(n) = 1.0;
	vecV = __calcOrthonorm( vecU, matVLocal, prm );
	if ( norm(vecV) > sqrt(eps) )
		strAction = sprintf( "Expand (e(%d)) vs VL", n );
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), strAction ) );
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__EXPAND_SUBSPACE;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 40;
		[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	endfor
	%
	%
	if (  1 <= sizeVLocal  &&  ( eta_loc < 0.99*omega || histDatSays_weShouldCutAndRun || sizeVLocal == sizeX )  )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Try VL*yLB" ) );
		vecY = matV'*(matVLocal*vecY_loc);
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__TRY_STEP;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 20;
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	if (  1 <= sizeVLocal  &&  ( eta_loc < 0.99*omega || histDatSays_weShouldCutAndRun || sizeVLocal == sizeX )  )
		msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( ...
		  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
		  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
		  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
		  sum(sumsq(matB)), "Try step (marginal cnvg)" ) );
		vecY = matV'*(matVLocal*vecY_loc);
		fModelDat.histDat(iterCount).iterDat.actionType = ACTION_TYPE__TRY_STEP;
		fModelDat.histDat(iterCount).iterDat.actionSubType = 20;
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
	  "  %4d: %4d / %4d / %dx%d;  %8.2e // (%8.2e) %8.2e / (%8.2e) %8.2e // %8.2e / %8.2e / %8.2e;  |%8.2e|:  %s.", ...
	  fModelDat.fevalCount, sizeVLocal, sizeV, sizeF, sizeX, ...
	  prm.omegaTol, eta_unb, eta_bnd, eta_locunb, eta_loc, omega, omega_prev, omega_initial, ...
	  sum(sumsq(matB)), "Give up (exhausted)" ) );
	retCode = RETCODE__ALGORITHM_BREAKDOWN;
	%
	return;
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
	if ( sizeV == sizeX )
	
		%vecY = -(matW'*matW)\(matW'*vecF);
		%vecRho = vecF + matW*vecY;
		%eta = sumsq(vecRho)/2.0
		%error( "HALT");
	
	
		% 2022-06-07: Another case to handle in light of preconditioner integration.
		if ( prm.valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecV,[sizeX,1]) );
			assert( abs(norm(vecV)-1.0) < sqrt(eps) );
			assert( norm(matVLocal'*vecV) < sqrt(eps) ); % Any reason to ever defy this?
			assert( sizeVLocal < sizeX );
		endif
		vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
		fevalIncr++;
		matW_updated = matW + ( vecW - matW*(matV'*vecV) ) * ( vecV'*matV );
		matVLocal_updated = [ matVLocal, vecV ];
		matWLocal_updated = [ matWLocal, vecW ];
		fModelDat.matW = matW_updated;
		fModelDat.matVLocal = matVLocal_updated;
		fModelDat.matWLocal = matWLocal_updated;
		%
		if ( prm.valdLev >= VALDLEV__MEDIUM )
			assert( reldiff(fModelDat.matW*fModelDat.matV'*vecV,vecW) < sqrt(eps) );
			assert( reldiff(fModelDat.matWLocal*fModelDat.matVLocal'*vecV,vecW) < sqrt(eps) );
			__validateFModelDat( fModelDat, prm );
		endif
		%
		retCode = RETCODE__SUCCESS;
		return;
		
	endif
	%
	%
	%
	if ( norm(matV'*vecV) > sqrt(eps) )
		% 2022-06-07: Lacking this is why 180 was doing funny things.
		% Doing something new.
		% Later, we can merge this with the two cases.
		if ( prm.valdLev >= VALDLEV__MEDIUM )
			assert( sizeV < sizeX );
			assert( isrealarray(vecV,[sizeX,1]) );
			assert( abs(norm(vecV)-1.0) < sqrt(eps) );
			assert( norm(matVLocal'*vecV) < sqrt(eps) ); %Don't require this either?
			assert( sizeVLocal <= sizeV );
		endif
		vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
		fevalIncr++;
		%
		vecVPerp = __calcOrthonorm( vecV, matV, prm );
		assert( abs(norm(vecVPerp)-1.0) < sqrt(eps) );
		matV_updated = [ matV, vecVPerp ];
		matW_temp = [ matW, zeros(sizeF,1) ];
		matW_updated = matW_temp + ( vecW - matW_temp*(matV_updated')*vecV ) * ( vecV'*matV_updated );
		fModelDat.matW = matW_updated;
		%
		%vecVPerpLocal = __calcOrthonorm( vecV, matVLocal, prm );
		%assert( abs(norm(vecVPerpLocal)-1.0) < sqrt(eps) );
		%matVLocal_updated = [ matVLocal, vecVPerpLocal ];
		%matWLocal_temp = [ matWLocal, zeros(sizeF,1) ];
		%matWLocal_updated = matWLocal_temp + ( vecW - matWLocal_temp*(matVLocal_updated')*vecV ) * ( vecV'*matVLocal_updated );
		matVLocal_updated = [ matVLocal, vecV ];
		matWLocal_updated = [ matWLocal, vecW ];
		%
		fModelDat.matV = matV_updated;
		fModelDat.matW = matW_updated;
		fModelDat.matB = zeros(sizeV+1,sizeV+1);
		fModelDat.matB(1:sizeV,1:sizeV) = matB;
		fModelDat.matVLocal = matVLocal_updated;
		fModelDat.matWLocal = matWLocal_updated;
		%
		if ( prm.valdLev >= VALDLEV__MEDIUM )
			assert( reldiff(fModelDat.matW*fModelDat.matV'*vecV,vecW) < sqrt(eps) );
			assert( reldiff(fModelDat.matWLocal*fModelDat.matVLocal'*vecV,vecW) < sqrt(eps) );
			__validateFModelDat( fModelDat, prm );
		endif
		%
		retCode = RETCODE__SUCCESS;
		return;
	endif
	%
	%
	if ( prm.valdLev >= VALDLEV__MEDIUM )
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
	vecX_cand = fModelDat.vecX_cand; % Earlier candidate for next point.
	vecF_cand = fModelDat.vecF_cand;
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
	eta = sumsq(vecFModel)/2.0;
	vecRho = vecF_trial - vecFModel;
	omega = sumsq(vecF)/2.0;
	omega_trial = sumsq(vecF_trial)/2.0;
	%
	if ( omega_trial < 0.1*eta )
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "       Trial was SIGNIFICANTLY better than model!" );
	elseif ( omega_trial < eta )
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "       Trial was better than model!" );
	endif
	%
	%
	if ( omega_trial > omega )
		sgnChar = "+";
	else
		sgnChar = "";
	endif
	if ( isempty(vecF_cand) )
		omega_cand_msg = 0.0;
	else
		omega_cand_msg = sumsq(vecF_cand)/2.0;
	endif
	%
	%
	if ( omega_trial >= omega )
		msgif ( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		  "       Try step: %8.2e -> %8.2e ( delta %c%8.2e, %c%8.2e rel ) vs model %8.2e, cand %8.2e: %s.", ...
		  omega, omega_trial, sgnChar, omega_trial-omega, sgnChar, omega_trial/omega - 1.0, eta, omega_cand_msg, ...
		  "Wholly rejecting trial; increasing B (reducing TR)" ) );
		%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		%  "  Wholly rejecting trial: %8.2e -> %8.2e ( up frac %8.2e, scale frac %8.2e ).", ...
		%  omega, omega_trial, omega_trial/omega - 1.0, omega_trial/omega ) );
		%
		%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Shrinking trust region." );
		trShrinkCoeff = 0.5; %%% Make param.
		[ retCode, fModelDat ] = __shrinkTR( trShrinkCoeff*vecY, fModelDat, prm );
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		return;
	endif
	%
	%
	if ( ~isempty(vecF_cand) )
		omega_cand = sumsq(vecF_cand)/2.0;
		assert( omega_cand < omega );
		if ( omega_trial > omega_cand )
			msgif ( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			  "       Try step: %8.2e -> %8.2e ( delta %c%8.2e, %c%8.2e rel ) vs model %8.2e, cand %8.2e: %s.", ...
			  omega, omega_trial, sgnChar, omega_trial-omega, sgnChar, omega_trial/omega - 1.0, eta, omega_cand_msg, ...
			  "Retroactively accepting earlier candidate" ) );
			%msgif( prm.verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, "  Trial is worse than earlier candidate." );
			%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			%  "  Accepting earlier candidate: %8.2e -> %8.2e vs %8.2e ( down frac %8.2e, remain frac %8.2e ).", ...
			%  omega, omega_cand, omega_trial, 1.0 - omega_trial/omega, omega_trial/omega ) );
			vecY_cand = matV'*(vecX_cand-vecX); % This is hack-ish, to work with pre-existing code.
			assert( reldiff(matV*vecY_cand,vecX_cand-vecX) < sqrt(eps) );
			[ retCode,  fModelDat ] = __moveTo( vecY_cand, vecF_cand, fModelDat, prm );
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
			% These values should have been already returned as best.
			return;
		endif
	endif
	%
	%
	relFallThresh = 0.5;  %%% Make param.
	omegaThresh = omega - relFallThresh * ( omega - eta_bnd ); % Require fall to be at least half of bound ideal.
	if ( omega_trial <= omegaThresh )
		msgif ( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		  "       Try step: %8.2e -> %8.2e ( delta %c%8.2e, %c%8.2e rel ) vs model %8.2e, cand %8.2e: %s.", ...
		  omega, omega_trial, sgnChar, omega_trial-omega, sgnChar, omega_trial/omega - 1.0, eta, omega_cand_msg, ...
		  "Accepting trial" ) );
		%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		%  "  Accepting trial: %8.2e -> %8.2e, leq %8.2e ( down frac %8.2e, remain frac %8.2e ).", ...
		%  omega, omega_trial, omegaThresh, 1.0 - omega_trial/omega, omega_trial/omega ) );
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		%
		trExpandCoeff = 1.5; %%% Make param.
		if ( trExpandCoeff*norm(matB*vecY) >= 1.0 )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "       Reducing B (expanding TR)." );
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
				
				%coastingFallRelThresh = 0.1; %%% Make param.
				%coastingC1 = 0.5;
				%if ( eta_coast < (1.0-coastingFallRelThresh)*omega && eta < (1.0-coastingC1)*omega )
				omega_initial = sumsq(fModelDat.vecF_initial)/2.0;
				omegaTolTemp = max([ 0.9*prm.omegaTol, 0.05*omega*min([ 1.0, omega/omega_initial ]) ]);
				if (  eta < omegaTolTemp  &&  eta_coast < 0.1 * omega  )
				
					vecX_trial = vecX + matV*vecY;
					vecF_trial = funchF( vecX_trial );
					tsFevalCount++;
					%
					vecFModel = vecF + matW * vecY;
					vecRho = vecF_trial - vecFModel;
					omega_trial = sumsq(vecF_trial)/2.0;
					%
					if ( omega_trial < omega )
						msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
						  "       Coasting was good ( %0.3e -> %0.3e ~ %0.3e? -> %0.3e ).", omega, eta, eta_coast, omega_trial ) );
						vecX_next = vecX_trial;
						vecF_next = vecF_trial;
						[ retCode, fModelDat ] = __moveTo( vecY, vecF_next, fModelDat, prm );
						if ( 0~= retCode )
							msgretcodeif( true, __FILE__, __LINE__, retCode );
							return;
						endif
						continue;
					else
						msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
						  "       Coasting was bad ( %0.3e -> %0.3e ~ %0.3e? -> %0.3e ).", omega, eta, eta_coast, omega_trial ) );
						break;
					endif
				else
					msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
					  "       Not attempting to coast ( %0.3e -> %0.3e ~ %0.3e? ).", omega, eta, eta_coast ) );
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
		%msg( __FILE__, __LINE__, "*** DOIN' DA NU THANG! ***" );
		assert( omega_trial < omega );
		assert( isempty(vecF_cand) || omega_trial <= sumsq(fModelDat.vecF_cand)/2.0 );
		assert( omega_trial > omegaThresh );
		msgif ( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		  "       Try step: %8.2e -> %8.2e ( delta %c%8.2e, %c%8.2e rel ) vs model %8.2e, cand %8.2e: %s.", ...
		  omega, omega_trial, sgnChar, omega_trial-omega, sgnChar, omega_trial/omega - 1.0, eta, omega_cand_msg, ...
		  "Temporarily rejecting trial; increasing B (reducing TR)" ) );
		%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		%  "  Temporarily rejecting trial: %8.2e -> %8.2e, gt %8.2e ( down frac %8.2e, remain frac %8.2e ).", ...
		%  omega, omega_trial, omegaThresh, 1.0 - omega_trial/omega, omega_trial/omega ) );
		fModelDat.vecF_cand = vecF_trial;
		fModelDat.vecX_cand = vecX_trial;
		% DRaburn 2022-05-31-1910:
		%  I've redefined "vecX" and "vecF" in main loop to be the "best", not neccesarily the current point.
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		%
		%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Shrinking trust region." );
		trShrinkCoeff = 0.5; %%% Make param.
		[ retCode, fModelDat ] = __shrinkTR( trShrinkCoeff*vecY, fModelDat, prm );
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		%
		[ retCode, studyDat ] = __studyFModel( fModelDat, prm );
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		%
		if ( studyDat.eta_loc >= omega_trial )
			%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			%  "  On second thought, accepting trial: %8.2e -> %8.2e (vs updated loc %8.2e) ( down frac %8.2e, remain frac %8.2e ).", ...
			%  omega, omega_trial, studyDat.eta_loc, 1.0 - omega_trial/omega, omega_trial/omega ) );
			msgif ( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			  "       On second thought, accepting trial (updated eta_loc = %8.2e vs omega_trial = %8.2e, omega = %8.2e ).", ...
			  studyDat.eta_loc, omega_trial, omega ) );
			[ retCode, fModelDat ] = __moveTo( vecY, vecF_trial, fModelDat, prm );
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
		else
			msgif ( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			  "       On second thought... still temporarily rejecting trial (updated eta_loc = %8.2e (vs omega_trial = %8.2e, omega = %8.2e ).", ...
			  studyDat.eta_loc, omega_trial, omega ) );
			%msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			%  "  On second thought, still rejecting trial: %8.2e -> %8.2e (vs updated loc %8.2e) ( down frac %8.2e, remain frac %8.2e ).", ...
			%  omega, omega_trial, studyDat.eta_loc, 1.0 - omega_trial/omega, omega_trial/omega ) );
		endif
		%
		return;
	endif
	error( "This line should be impossible to reach." );
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
	if ( prm.valdLev >= VALDLEV__MEDIUM )
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
	
	if (0)
		msg( __FILE__, __LINE__, "Infodump..." );
	endif
	
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
	vecY_unb = __findLevPt( vecWTF, matWTW, matC, [], [], prm );
	vecY_bnd = __findLevPt( vecWTF, matWTW, matC, matB, 1.0, prm );
	eta_unb = sumsq( vecF + matW*vecY_unb )/2.0;
	eta_bnd = sumsq( vecF + matW*vecY_bnd )/2.0;
	b_unb = norm(matB*vecY_unb);
	b_bnd = norm(matB*vecY_bnd);
	if ( isempty(matVLocal) )
		vecY_locunb = [];
		eta_locunb = omega;
		b_locunb = 0.0;
		vecY_loc = [];
		eta_loc = omega;
		b_loc = 0.0;
	else
		matWTWLocal = matWLocal'*matWLocal;
		vecWTFLocal = matWLocal'*vecF;
		matBLocal = matB*(matV'*matVLocal);
		matCLocal = (matVLocal'*matV)*matC*(matV'*matVLocal);
		matCLocal = (matCLocal'+matCLocal);
		vecY_locunb = __findLevPt( vecWTFLocal, matWTWLocal, matCLocal, [], [], prm );
		eta_locunb = sumsq( vecF + matWLocal*vecY_locunb )/2.0;
		b_locunb = norm(matBLocal*vecY_locunb);
		vecY_loc = __findLevPt( vecWTFLocal, matWTWLocal, matCLocal, matBLocal, 1.0, prm );
		eta_loc = max([ 0.0, omega + vecWTFLocal'*vecY_loc + abs((vecY_loc'*matWTWLocal*vecY_loc)/2.0) ]);
		%eta_loc = sumsq( vecF + matWLocal*vecY_loc )/2.0;
		b_loc = norm(matBLocal*vecY_loc);
	endif
	%
	%
	studyDat.vecY_unb = vecY_unb; % Record unbound.
	studyDat.vecY_bnd = vecY_bnd; % Record bound.
	studyDat.vecY_locunb = vecY_locunb; % Local unbound.
	studyDat.vecY_loc = vecY_loc; % Local bound.
	studyDat.eta_unb = eta_unb;
	studyDat.eta_bnd = eta_bnd;
	studyDat.eta_locunb = eta_locunb;
	studyDat.eta_loc = eta_loc;
	%
	if ( prm.verbLev >= VERBLEV__COPIOUS )
		msg( __FILE__, __LINE__, sprintf( ...
		  "   eta: %8.2e / %8.2e / %8.2e / %8.2e (/ %8.2e);  y: %8.2e / %8.2e / %8.2e;  b: %8.2e / %8.2e / %8.2e.", ...
		  eta_unb, eta_bnd, eta_loc, omega, prm.omegaTol, ...
		  norm(vecY_unb), norm(vecY_bnd), norm(vecY_loc), ...
		  b_unb, b_bnd, b_loc ) );
	endif
	%
	%if ( prm.valdLev >= VALDLEV__LOW )
	%	__validateStudyDat( fModelDat, studyDat, prm );
	%endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function vecY = __findLevPt( vecG, matH, matC, matB, bTrgt, prm )
	if ( prm.useDogLeg )
		vecY = __dogLeg( vecG, matH, matC, matB, bTrgt, prm );
	else
		vecY = findLevPt_0527( vecG, matH, matC, matB, bTrgt, prm.findLevPrm );
	endif
	return;
endfunction
function vecY = __dogLeg( vecG, matH, matC, matB, bTrgt, prm )
	% I'm not 100% sure this Powell's dog leg is correct.
	% Regardless, this code does not seem to help.
	mydefs;
	DO_HACKS = false;
	if (DO_HACKS)
		vecY_wouldHaveBeen = findLevPt_0527( vecG, matH, matC, matB, bTrgt, prm.findLevPrm );
	endif
	vecYNewton = findLevPt_0527( vecG, matH, matC, [], [], prm.findLevPrm );
	%
	% If Newton step is in bounds, take it.
	if ( isempty(matB) || sumsq(matB*vecYNewton) <= 1.0 )
		vecY = vecYNewton;
		return;
	endif
	msgif( prm.verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Restricting step to bound (%g).",  sumsq(matB*vecYNewton) ) );
	%
	matRC = chol(matC);
	vecYCauchyDir = matRC \ (matRC'\(-vecG));
	%
	ythy = vecYCauchyDir'*matH*vecYCauchyDir;
	assert( ythy >= 0.0 );
	ythmg = -vecYCauchyDir'*vecG;
	assert( ythmg > 0.0 );
	pCauchy = ythmg / ythy;
	vecYCauchy = pCauchy*vecYCauchyDir;
	%
	% If Cauchy step goes OOB, take its intersection with boundary.
	vecBC = matB*vecYCauchy;
	bcsq = sumsq( vecBC );
	if ( bcsq >= 1.0 );
		vecY = vecYCauchy / sqrt(bcsq);
		assert( reldiff( sumsq(matB*vecY), 1.0 ) < sqrt(eps) );
		if (DO_HACKS)
			eta = vecG'*vecY + (vecY'*matH*vecY)/2.0;
			eta_wouldHaveBeen = vecG'*vecY_wouldHaveBeen + (vecY_wouldHaveBeen'*matH*vecY_wouldHaveBeen)/2.0;
			[ norm(vecY), norm(vecY_wouldHaveBeen), norm(vecY-vecY_wouldHaveBeen) ]
			[ eta, eta_wouldHaveBeen, eta-eta_wouldHaveBeen ]
			%error( "BEHOLD!" );
		endif
		return;
	endif
	%
	% Find where the second leg intersects the boundary.
	vecY2 = vecYNewton - vecYCauchy;
	vecB2 = matB*vecY2;
	% We can't use calcLinishRootOfQuad() here!
	% We want positive root of quad!
	%   Using y = yCauchy + t * y2,
	%   ||B*y|| = (t^2)*(b2'*b2) + t*(2.0*b2'*bc) + (bc'*bc) - 1.0,
	%   where b2 = B*y2 and bc = B*yCauchy.
	a = sumsq(vecB2);
	b = 2.0*(vecBC'*vecB2);
	c = bcsq-1.0;
	discrim = (b^2) - (4.0*a*c);
	assert( discrim >= 0.0 );
	assert( 0.0 < a );
	t = (-b+sqrt(discrim))/(2.0*a); % Because a must be positive.
	vecY = vecYCauchy + (t*vecY2);
	assert( reldiff( sumsq(matB*vecY), 1.0 ) < sqrt(eps) );
	assert( t >= 0.0 );
	assert( t <= 1.0 );
	if (DO_HACKS)
		eta = vecG'*vecY + (vecY'*matH*vecY)/2.0;
		eta_wouldHaveBeen = vecG'*vecY_wouldHaveBeen + (vecY_wouldHaveBeen'*matH*vecY_wouldHaveBeen)/2.0;
		[ norm(vecY), norm(vecY_wouldHaveBeen), norm(vecY-vecY_wouldHaveBeen) ]
		[ eta, eta_wouldHaveBeen, eta-eta_wouldHaveBeen ]
		%error( "BEHOLD!" );
	endif
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
	%vecX_cand = fModelDat.vecX_cand;
	vecF_cand = fModelDat.vecF_cand;
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
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		assert( 0.0 < yNorm )
		assert( norm(vecFModel_next) <= norm(vecF) );
		assert( norm(vecF_next) <= norm(vecF) );
		if ( ~isempty(vecF_cand) )
			assert( norm(vecF_next) <= norm(vecF_cand) );
		endif
	endif
	useQuadUpdate = true;
	if (useQuadUpdate)
		matW_updated = matW + 2.0 * vecRhoF * (vecY')/(yNorm^2);
	else
		matW_updated = matW + vecRhoF * (vecY')/(yNorm^2);
	endif
	%
	%
	fModelDat.vecX_prev = vecX;
	fModelDat.vecF_prev = vecF;
	%
	fModelDat.vecX = vecX + matV*vecY;
	fModelDat.vecF = vecF_next;
	%fModelDat.matV = matV;
	fModelDat.matW = matW_updated;
	%fModelDat.matB = matB;
	fModelDat.matVLocal = zeros(sizeX,0);
	fModelDat.matWLocal = zeros(sizeF,0);
	fModelDat.vecX_cand = [];
	fModelDat.vecF_cand = [];
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
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		assert( 0.0 < v );
	endif
	vecW = ( funchF(vecX+prm.epsFD*vecV) - vecF ) / prm.epsFD;
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
	endif
	%
	if ( prm.valdLev >= VALDLEV__MEDIUM )
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
	return;
endfunction


function __validateFModelDat( fModelDat, prm )
	mydefs;
	if ( prm.valdLev < VALDLEV__MEDIUM )
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
