% High Priority:
% TODO: Compare to zlinsolf100 and get working comprably well.
% TODO: Compare to findZero800 and get working comprably well.
% TODO:  Conventional AP benefits from J ~ I. We best sense and exploit J ~ D.
% TODO: Consider more desperation moves, such as expanding per any matW column or normal basis vector.

% Lower Priority:
% TODO: Write __tryStep() with re-try and "latest trial is worse than previous trial!" handling.
% TODO: Update epsFD based on steps taken.
% TODO: Add scaling-matrix preconditioning.
% TODO: Allow for evaluating a subspace vector that is not entirely outside matV if the preconditioner seems highly reliable.
% TODO: Bad-local min handling + jump.
% NEXT VER: Move "Expand subspace (init)" to __initFModel()?

% Consider:
% TODO: Avoid __studyFModel() redundancy in __tryStep(),
%       perhaps by making __moveTo() considered as a normal action
%       or passing the calculated studyDat out for the next iteration;
%       Or, update studyDat at END of loop?
% TODO: Sepratrix-domain-enumeration for basis bad local min.

function [ vecX, vecF, retCode, fevalCount, stepsCount, datOut ] = zlinsolf195( funchF, vecX_initial, vecF_initial=[], prmIn=[] )
	% INIT
	mydefs;
	startTime = time();
	if ( stopsignalpresent() )
		msg(__FILE__, __LINE__, "ERROR: Stop signal already present." );
		retCode = RETCODE__IMPOSED_STOP;
		return;
	endif
	vecX = [];
	vecF = [];
	retCode = RETCODE__NOT_SET;
	fevalCount = 0;
	stepsCount = 0;
	datOut = [];
	[ retCode, fevalIncr, vecF_initial, fModelDat, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn );
	fevalCount += fevalIncr; clear fevalIncr;
	if ( 0~= retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	if ( prm.verbLev >= VERBLEV__INFO )
		msg( __FILE__, __LINE__, "Limits..." );
		msg( __FILE__, __LINE__, sprintf( ...
		  "   time: %8.2e;  iter: %3d;  feval: %3d;  steps: %3d;  size: %3d x %3d;  omega: %8.2e.", ...
		  prm.timeMax, ...
		  prm.iterMax, ...
		  prm.fevalMax, ...
		  prm.stepsMax, ...
		  size(vecX_initial,1), size(vecF_initial,1), ...
		  prm.omegaTol ) );
	endif
	%
	%
	% PREP FOR MAIN LOOP
	iterCount = 0;
	vecX = vecX_initial;
	vecF = vecF_initial;
	omega = sumsq(vecF)/2.0;
	datOut.fevalCountOfSteps(stepsCount+1) = fevalCount;
	datOut.fNormOfSteps(stepsCount+1) = norm(vecF);
	datOut.vecXOfSteps(:,stepsCount+1) = vecX;
	datOut.vecFOfSteps(:,stepsCount+1) = vecF;
	%
	[ retCode, fevalIncr, fModelDat ] = __initFModel( funchF, vecX, vecF, prm );
	fevalCount += fevalIncr; clear fevalIncr;
	if ( 0~=retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	%
	%
	% MAIN LOOP
	while (1)
		% Tier 0a actions: simple stoping criteria.
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
		endif
		%
		%
		iterCount++;
		%
		%
		if ( prm.verbLev >= VERBLEV__DETAILS )
			%msg( __FILE__, __LINE__, "Progress..." );
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
				  " Step %3d ( at %8.2e, %3d, %3d with %8.2e x %8.2e ):  %8.2e -> %8.2e ( down %8.2e; %8.2e ).", ...
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
	if ( prm.verbLev >= VERBLEV__INFO )
		msg( __FILE__, __LINE__, "Final..." );
		msg( __FILE__, __LINE__, sprintf( ...
		  "   time: %9.2e;  iter: %3d;  feval: %3d;  steps: %3d;  size: %3d / %3d ( / %d x %d );  omega: %8.2e.", ...
		  time()-startTime, ...
		  iterCount, ...
		  fevalCount, ...
		  stepsCount, ...
		  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX_initial,1), size(vecF_initial,1), ...
		  sumsq(vecF)/2.0 ) );
	endif
	if (1)
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
	%prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__LOW; % "Production / optimization".
	%prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Late integration.
	prm.verbLev = VERBLEV__PROGRESS; prm.valdLev = VALDLEV__HIGH; % Early integration.
	%prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Dev / refine.
	%prm.verbLev = VERBLEV__COPIOUS; prm.valdLev = VALDLEV__VERY_HIGH; % Deep dev.
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Debug.
	%
	prm.timeMax = -1.0; %%%
	%%%prm.timeMax = 10.0; %%%
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
	prm.moveToELoCoeff = 0.9;
	prm.moveToEHiCoeff = 0.1;
	%
	prm.initExpandExp = 0.5;
	prm.expandRelThresh = 0.5;
	prm.stallRelThresh = 1.0E-4;
	prm.tryFallRelThresh = 0.9;
	prm.acceptFallRelThresh = 0.1;
	prm.reevalFallRelThresh = 0.1;
	%
	prm.trExpandCoeff = 1.5;
	prm.trShrinkCoeff = 0.5;
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
	fModelDat.vecX_initial = vecX;
	fModelDat.vecF_initial = vecF;
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matALo = [ 0.0 ];
	fModelDat.matAHi = [ 0.0 ];
	fModelDat.matB = [ 0.0 ];
	%%%fModelDat.matB = [ 10.0 ]; %%%
	fModelDat.matVLocal = [ vecV ];
	fModelDat.strState = "init";
	%
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
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%strState = fModelDat.strState; % String indicating "state" of fModelDat.
	%
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	matWTW = matW'*matW;
	matD = diag(diag(matWTW));
	vecWTF = matW'*vecF;
	matBTB = matB'*matB;
	matIV = eye(sizeV,sizeV);
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
	vecY_ideal = __findCandStep( matWTW, vecWTF, matC, [], [], prm );
	vecY_zeroV = __findCandStep( matWTW, vecWTF, matC, matB, 1.0, prm );
	
	
	hack0529 = true;
	if (hack0529)
	vecY_loVar = vecY_zeroV;
	vecY_hiVar = vecY_zeroV;
	funchEta_zeroV = @(y)( max([ 0.0, sumsq(vecF)/2.0 + vecWTF'*y + abs((y'*matWTW*y)/2.0) ]) );
	funchEta_loVar = @(y)( sumsq(vecF)*10.0 + norm(vecWTF)*norm(y)*10.0 + norm(y)*norm(matWTW*y)*10.0 + norm(y)*norm(matALo*y)*10.0 );
	funchEta_hiVar = @(y)( sumsq(vecF)*10.0 + norm(vecWTF)*norm(y)*10.0 + norm(y)*norm(matWTW*y)*10.0 + norm(y)*norm(matAHi*y)*10.0 );
	else
	vecY_loVar = __findCandStep( matWTW + matALo, vecWTF, matC, matB, 1.0, prm );
	vecY_hiVar = __findCandStep( matWTW + matAHi, vecWTF, matC, matB, 1.0, prm );
	%
	% "eta" is estimate for cost function;
	% "omega" is observed values.
	funchEta_zeroV = @(y)( max([ 0.0, sumsq(vecF)/2.0 + vecWTF'*y + abs((y'*matWTW*y)/2.0) ]) );
	funchEta_loVar = @(y)( max([ 0.0, sumsq(vecF)/2.0 + vecWTF'*y + abs((y'*matWTW*y)/2.0) + abs((y'*matALo*y)/2.0) ]) );
	funchEta_hiVar = @(y)( max([ 0.0, sumsq(vecF)/2.0 + vecWTF'*y + abs((y'*matWTW*y)/2.0) + abs((y'*matAHi*y)/2.0) ]) );
	endif
	
	
	%
	if ( prm.verbLev >= VERBLEV__DETAILS+10 )
		if ( 0 == sizeVLocal )
			funch_yLocNorm = @(y)( 0.0 );
			funch_bLocNorm = @(y)( 0.0 );
		else
			funch_yLocNorm = @(y)( norm(matVLocal'*(matV*y)) );
			funch_bLocNorm = @(y)( norm(matB*(matV'*(matVLocal*(matVLocal'*(matV*y))))) );
		endif
		msg( __FILE__, __LINE__, sprintf( "   vecY_ideal: %8.2e / %8.2e / %8.2e ( %8.2e, %8.2e;  %8.2e, %8.2e ).", ...
		  funchEta_zeroV(vecY_ideal), funchEta_loVar(vecY_ideal), funchEta_hiVar(vecY_ideal), ...
		  norm(vecY_ideal), norm(matB*vecY_ideal), funch_yLocNorm(vecY_ideal), funch_bLocNorm(vecY_ideal) ) );
		msg( __FILE__, __LINE__, sprintf( "   vecY_zeroV: %8.2e / %8.2e / %8.2e ( %8.2e, %8.2e;  %8.2e, %8.2e ).", ...
		  funchEta_zeroV(vecY_zeroV), funchEta_loVar(vecY_zeroV), funchEta_hiVar(vecY_zeroV), ...
		  norm(vecY_zeroV), norm(matB*vecY_zeroV), funch_yLocNorm(vecY_zeroV), funch_bLocNorm(vecY_zeroV) ) );
		msg( __FILE__, __LINE__, sprintf( "   vecY_loVar: %8.2e / %8.2e / %8.2e ( %8.2e, %8.2e;  %8.2e, %8.2e ).", ...
		  funchEta_zeroV(vecY_loVar), funchEta_loVar(vecY_loVar), funchEta_hiVar(vecY_loVar), ...
		  norm(vecY_loVar), norm(matB*vecY_loVar), funch_yLocNorm(vecY_loVar), funch_bLocNorm(vecY_loVar) ) );
		msg( __FILE__, __LINE__, sprintf( "   vecY_hiVar: %8.2e / %8.2e / %8.2e ( %8.2e, %8.2e;  %8.2e, %8.2e ).", ...
		  funchEta_zeroV(vecY_hiVar), funchEta_loVar(vecY_hiVar), funchEta_hiVar(vecY_hiVar), ...
		  norm(vecY_hiVar), norm(matB*vecY_hiVar), funch_yLocNorm(vecY_hiVar), funch_bLocNorm(vecY_hiVar) ) );
	endif
	%
	%
	%
	studyDat.matWTW = matWTW;
	studyDat.matD = matD;
	studyDat.vecWTF = vecWTF;
	studyDat.matBTB = matBTB;
	studyDat.matC = matC;
	%
	studyDat.vecY_ideal = vecY_ideal;
	studyDat.vecY_zeroV = vecY_zeroV;
	studyDat.vecY_loVar = vecY_loVar;
	studyDat.vecY_hiVar = vecY_hiVar;
	%
	% "eta" is estimate for cost function;
	% "omega" is observed values.
	studyDat.funchEta_zeroV = funchEta_zeroV;
	studyDat.funchEta_loVar = funchEta_loVar;
	studyDat.funchEta_hiVar = funchEta_hiVar;
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateStudyDat( fModelDat, studyDat, prm );
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
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	strState = fModelDat.strState; % String indicating "state" of fModelDat.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	vecY_ideal = studyDat.vecY_ideal;
	vecY_zeroV = studyDat.vecY_zeroV;
	vecY_loVar = studyDat.vecY_loVar;
	vecY_hiVar = studyDat.vecY_hiVar;
	funchEta_zeroV = studyDat.funchEta_zeroV;
	funchEta_loVar = studyDat.funchEta_loVar;
	funchEta_hiVar = studyDat.funchEta_hiVar;
	%
	%
	%
	omega = sumsq(vecF)/2.0;
	omega_initial = sumsq(vecF)/2.0;
	
	
	hack0529 = true;
	if (hack0529)
		clear vecY_loVar;
		clear vecY_hiVar;
		clear funchEta_loVar;
		clear funchEta_hiVar;
		%
		%
		omegaDynaThresh = 0.1 * omega * ( omega / omega_initial )^0.5;
		if ( funchEta_zeroV(vecY_ideal) > max([ 0.1*omegaDynaThresh, prm.omegaTol ]) )
			vecRhoF = matW(:,end);
			vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
			vecV = __calcOrthonorm( vecU, matV, prm );
			if ( norm(vecV) > sqrt(eps) )
				msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, " Action: Expand subspace." );
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
		%
		%
		%%%if ( funchEta_hiVar(vecY_zeroV) <= max([ omegaDynaThresh, prm.omegaTol ]) )
		if ( funchEta_zeroV(vecY_zeroV) <= max([ omegaDynaThresh, prm.omegaTol ]) ) %%%
		vecYLocal = matVLocal'*(matV*vecY_zeroV); %%%
		if ( norm(vecYLocal) > (1.0-sqrt(eps))*norm(vecY_zeroV) ) %%%
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( " Action: Try step ( approach %9.2e -> %9.2e ).", ...
			 omega, funchEta_zeroV(vecY_zeroV) ) );
			[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY_zeroV, funchF, fModelDat, studyDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif %%%
		endif
		%
		%
		if (1)
			vecV = __calcOrthonorm( matV*vecY_zeroV, matVLocal, prm );
			if ( norm(vecV) > sqrt(eps) )
				vecV = matV*(matV'*vecV); % Force in subspace for numerical stability..
				msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, " Action: Re-evaluate direction." );
				[ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm );
				taFevalCount += fevalIncr; clear fevalIncr;
				if ( 0~= retCode )
					msgretcodeif( true, __FILE__, __LINE__, retCode );
				endif
				return;
			endif
		endif
		%
		%
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "No good action; giving up." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	endif
	
	
	%
	%
	if (1)
	% Initially, be fairly aggressive in expanding subspace.
	if ( strcmpi(strState,"init") ) % CAUTION: MATLAB strmpi is opposite of C/C++.
	%%%if ( funchEta_zeroV(vecY_ideal) > prm.omegaTol * ( ( omega / prm.omegaTol ) ^ (1.0-prm.initExpandExp) ) )
	%%%if ( funchEta_zeroV(vecY_ideal) > max([ prm.omegaTol, prm.expandRelThresh^2 * omega]) )
	if ( funchEta_zeroV(vecY_ideal) > prm.expandRelThresh * prm.omegaTol * ( ( omega / prm.omegaTol ) ^ 0.9 ) )
		vecRhoF = matW(:,end);
		vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
		vecV = __calcOrthonorm( vecU, matV, prm );
		if ( norm(vecV) > sqrt(eps) )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, " Action: Expand subspace (init)." );
			[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif
	endif
	endif
	endif
	%
	%
	%%%omega_initial = sumsq(fModelDat.vecF_initial)/2.0;
	%%%expandOmegaTol = max([ prm.omegaTol, ...
	%%%  prm.expandRelThresh * prm.omegaTol * ( omega / prm.omegaTol )^0.1, ...
	%%%  prm.expandRelThresh * omega * (omega/omega_initial)^0.1 ]);
	if ( funchEta_zeroV(vecY_ideal) > max([ prm.omegaTol, prm.expandRelThresh * omega]) )
	%%%if ( funchEta_zeroV(vecY_ideal) > max([ prm.omegaTol, prm.expandRelThresh * omega * sqrt(omega/omega_initial) ]) )
	%%%if ( funchEta_zeroV(vecY_ideal) > expandOmegaTol )
		vecRhoF = vecF - matW * vecY_ideal;
		vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
		vecV = __calcOrthonorm( vecU, matV, prm );
		if ( norm(vecV) > sqrt(eps) )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, " Action: Expand subspace (ongoing)." );
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
	%%%if ( funchEta_hiVar(vecY_zeroV) <= prm.omegaTol )
	if ( funchEta_zeroV(vecY_zeroV) <= prm.omegaTol )
		%error( "TODO: Try striking at vecY_zeroV." );
		%% Note: __tryStep() may internally update fModelDat and call __studyFModel(),
		%%  making the next itertion's call to __studyFModel() redundant. POITROME.
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( " Action: Try step ( strike %9.2e -> %9.2e / %9.2e / %9.2e ).", ...
		 omega, funchEta_zeroV(vecY_zeroV), funchEta_loVar(vecY_zeroV), funchEta_hiVar(vecY_zeroV) ) );
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY_zeroV, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	if ( funchEta_zeroV(vecY_zeroV) > omega * ( 1.0 - prm.stallRelThresh ) )
		error( "TODO: If zeroV says little decrease is possible, give up or do BLM handling." );
	endif
	%
	%
	if ( funchEta_hiVar(vecY_hiVar) <= omega + prm.tryFallRelThresh * ( funchEta_zeroV(vecY_zeroV) - omega ) )
	%%%if ( funchEta_hiVar(vecY_zeroV) <= omega + prm.tryFallRelThresh * ( funchEta_zeroV(vecY_zeroV) - omega ) )
		%error( "TODO: Something like __tryStep( vecY_hiVar, funchF, fModelDat, prm );" );
		%% Note: __tryStep() may internally update fModelDat and call __studyFModel(),
		%%  making the next itertion's call to __studyFModel() redundant. POITROME.
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( " Action: Try step ( approach %9.2e -> %9.2e / %9.2e / %9.2e ).", ...
		 omega, funchEta_zeroV(vecY_hiVar), funchEta_loVar(vecY_hiVar), funchEta_hiVar(vecY_hiVar) ) );
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY_hiVar, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	%%%if ( funchEta_loVar(vecY_loVar) >= omega + prm.reevalFallRelThresh * ( funchEta_zeroV(vecY_zeroV) - omega ) )
	if (1)
		vecV = __calcOrthonorm( matV*vecY_loVar, matVLocal, prm );
		%%%vecV = __calcOrthonorm( matV*vecY_zeroV, matVLocal, prm );
		if ( norm(vecV) > sqrt(eps) )
			vecV = matV*(matV'*vecV); % Force in subspace for numerical stability..
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, " Action: Re-evaluate direction." );
			[ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm );
			taFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
			endif
			return;
		endif
	endif
	%
	%
	if (1)
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( " Action: Try step, beacuse not sure what else to do ( strike %9.2e -> %9.2e / %9.2e / %9.2e ).", ...
		 omega, funchEta_zeroV(vecY_zeroV), funchEta_loVar(vecY_zeroV), funchEta_hiVar(vecY_zeroV) ) );
		[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY_zeroV, funchF, fModelDat, studyDat, prm );
		taFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
		endif
		return;
	endif
	%
	%
	error( "TODO: Consider more desperate actions, particularly expand or reeval with different vectors." );
	return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "ACTION" FUNCTIONS
%

% TODO: Allow for evaluating a subspace vector that is not entirely outside matV
% if the preconditioner seems highly reliable.
function [ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm )
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
	%strState = fModelDat.strState; % String indicating "state" of fModelDat.
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
	fModelDat.matALo = zeros(sizeV+1,sizeV+1);
	fModelDat.matALo(1:sizeV,1:sizeV) = matALo;
	fModelDat.matAHi = zeros(sizeV+1,sizeV+1);
	fModelDat.matAHi(1:sizeV,1:sizeV) = matAHi;
	fModelDat.matB = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = matB;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	%fModelDat.strState = strState;
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
	retCode = RETCODE__SUCCESS;
	return;
endfunction


%TODO: Implement a non-_crude() __tryStep().
function [ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm )
	[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep_crude( vecY, funchF, fModelDat, studyDat, prm );
	if ( 0~= retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	return;
endfunction
function [ retCode, tsFevalCount, fModelDat, vecX_next, vecF_next ] = __tryStep_crude( vecY, funchF, fModelDat, studyDat, prm )
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
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	strState = fModelDat.strState; % String indicating "state" of fModelDat.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	%
	vecY_ideal = studyDat.vecY_ideal;
	vecY_zeroV = studyDat.vecY_zeroV;
	vecY_loVar = studyDat.vecY_loVar;
	vecY_hiVar = studyDat.vecY_hiVar;
	funchEta_zeroV = studyDat.funchEta_zeroV;
	funchEta_loVar = studyDat.funchEta_loVar;
	funchEta_hiVar = studyDat.funchEta_hiVar;
	%
	vecX_trial = vecX + matV * vecY;
	vecF_trial = funchF( vecX_trial );
	tsFevalCount++;
	%
	omega = sumsq(vecF)/2.0;
	omega_trial = sumsq(vecF_trial)/2.0;


	hack0529 = true;
	if (hack0529)
		clear vecY_loVar;
		clear vecY_hiVar;
		clear funchEta_loVar;
		clear funchEta_hiVar;
		%
		omgaAcceptThresh = 0.5*omega;
		if ( omega_trial < omgaAcceptThresh)
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
			  "  Accepting step: %9.2e -> %9.2e < %9.2e ( down frac %9.2e, remain frac %9.2e ).", ...
			  omega, omega_trial, omgaAcceptThresh, 1.0 - omega_trial/omega, omega_trial/omega ) );
			vecX_next = vecX_trial;
			vecF_next = vecF_trial;
			if ( prm.trExpandCoeff*norm(matB*vecY) >= 1.0 )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Expanding trust region." );
				[ retCode, fevalIncr, fModelDat ] = __expandTR( funchF, prm.trExpandCoeff*vecY, fModelDat, prm );
				tsFevalCount += fevalIncr; clear fevalIncr;
				if ( 0~= retCode )
					msgretcodeif( true, __FILE__, __LINE__, retCode );
					return;
				endif
			endif
			[ retCode, fevalIncr, fModelDat ] = __moveTo( vecY, vecF_next, funchF, fModelDat, prm );
			tsFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
			retCode = RETCODE__SUCCESS;
			return;
		endif
		%
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		  "  Rejecting step: %9.2e -> %9.2e > %9.2e ( up frac %9.2e, scale frac %9.2e ).", ...
		  omega, omega_trial, omgaAcceptThresh, omega_trial/omega - 1.0, omega_trial/omega ) );
		%
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Shrinking trust region." );
		% Add barrier to prevent re-trial.
		[ retCode, fevalIncr, fModelDat ] = __shrinkTR( funchF, prm.trShrinkCoeff*vecY, fModelDat, prm );
		tsFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		return;
	endif
	
	
	%
	omgaAcceptThresh = omega + prm.acceptFallRelThresh * ( funchEta_zeroV(vecY_zeroV) - omega );
	if ( omega_trial < omgaAcceptThresh )
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
		  "  Accepting step: %9.2e -> %9.2e < %9.2e ( down frac %9.2e, remain frac %9.2e ).", ...
		  omega, omega_trial, omgaAcceptThresh, 1.0 - omega_trial/omega, omega_trial/omega ) );
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		if ( prm.trExpandCoeff*norm(matB*vecY) >= 1.0 )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Expanding trust region." );
			[ retCode, fevalIncr, fModelDat ] = __expandTR( funchF, prm.trExpandCoeff*vecY, fModelDat, prm );
			tsFevalCount += fevalIncr; clear fevalIncr;
			if ( 0~= retCode )
				msgretcodeif( true, __FILE__, __LINE__, retCode );
				return;
			endif
		endif
		[ retCode, fevalIncr, fModelDat ] = __moveTo( vecY, vecF_next, funchF, fModelDat, prm );
		tsFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		retCode = RETCODE__SUCCESS;
		return;
	endif
	msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, sprintf( ...
	  "  Rejecting step: %9.2e -> %9.2e > %9.2e ( up frac %9.2e, scale frac %9.2e ).", ...
	  omega, omega_trial, omgaAcceptThresh, omega_trial/omega - 1.0, omega_trial/omega ) );
	%
	didReeval = false;
	assert( abs(norm(matV*vecY) - norm(vecY)) < sqrt(eps) );
	vecU = matV*vecY;
	vecV = __calcOrthonorm( vecU, matVLocal, prm );
	if ( norm(vecV) >= 1.0-sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Re-ealuating (component) of step direction." );
		vecV = matV*(matV'*vecV); % Force in subspace for numerical stability.
		[ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm );
		tsFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
		didReeval = true;
	endif
	%
	%%%vecFModel = vecF + matW * vecY;
	%%%if ( norm(vecFModel) < norm(vecF) )
	if (1)
		msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "  Shrinking trust region." );
		% Add barrier to prevent re-trial.
		[ retCode, fevalIncr, fModelDat ] = __shrinkTR( funchF, prm.trShrinkCoeff*vecY, fModelDat, prm );
		tsFevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			return;
		endif
	elseif (~didReeval)
		assert( "Neither re-eval nor shrink." );
	endif
	%
	retCode = RETCODE__SUCCESS;
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
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%strState = fModelDat.strState; % String indicating "state" of fModelDat.
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
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
	matIV = eye(sizeV,sizeV);
	vecYHat = matV'*vecV;
	assert( abs(norm(vecYHat)-1.0) < sqrt(eps) );
	matW_updated = matW + ( vecW - matW*vecYHat )*(vecYHat');
	matALo_updated = ( matIV - vecYHat*(vecYHat') ) * matALo * ( matIV - vecYHat*(vecYHat') );
	matAHi_updated = ( matIV - vecYHat*(vecYHat') ) * matAHi * ( matIV - vecYHat*(vecYHat') );
	matALo_updated = (matALo_updated'+matALo_updated)/2.0;
	matAHi_updated = (matAHi_updated'+matAHi_updated)/2.0;
	[ matRLo, cholFlag ] = chol(matALo_updated);
	if (0)
	if ( 0 ~= cholFlag || min(diag(matRLo)) <= prm.cholRelTol * max(abs(diag(matRLo))) )
		matALo_updated += max(max(abs(matALo_updated)))*sqrt(eps)*matIV;
	endif
	[ matRHi, cholFlag ] = chol(matAHi_updated);
	if ( 0 ~= cholFlag || min(diag(matRHi)) <= prm.cholRelTol * max(abs(diag(matRHi))) )
		matAHi_updated += max(max(abs(matAHi_updated)))*sqrt(eps)*matIV;
	endif
	elseif (1)
	if ( min(diag(matALo_updated)) < 0.0 )
		matALo_updated += abs(min(diag(matALo_updated)))*(1.0+sqrt(eps))*matIV;
	endif
	if ( min(diag(matAHi_updated)) < 0.0 )
		matAHi_updated += abs(min(diag(matAHi_updated)))*(1.0+sqrt(eps))*matIV;
	endif
	endif
	%
	%
	%
	%fModelDat.vecX = vecX;
	%fModelDat.vecF = vecF;
	%fModelDat.matV = matV;
	fModelDat.matW = matW_updated;
	fModelDat.matALo = matALo_updated;
	fModelDat.matAHi = matAHi_updated;
	%fModelDat.matB = matB;
	fModelDat.matVLocal = [ matVLocal, vecV ];
	%fModelDat.strState = strState;
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		__validateFModelDat( fModelDat, prm );
	endif
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
	%strState = fModelDat.strState; % String indicating "state" of fModelDat.
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
	if (~useQuadUpdate)
		matW_updated = matW + vecRhoF * (vecY')/(yNorm^2);
	else
		rhoSumsq = sumsq(vecRhoF);
		ytay = vecY'*matAHi*vecY;
		if ( rhoSumsq <= ytay )
			s = 1.0;
		else
			s = ytay/rhoSumsq;
		endif
		matW_updated = matW + (2.0-s) * vecRhoF * (vecY')/(yNorm^2);
	endif
	%
	vecDW = sum(matW.^2,1)';
	vecD1 = ones(sizeV,1)/sumsq(vecY);
	vecD2 = vecDW/sumsq(matW*vecY);
	vecD3  = vecDW/(vecY'*diag(vecDW)*vecY);
	vecDLo = min( [ vecD1, vecD2, vecD3 ], [], 2 );
	vecDHi = max( [ vecD1, vecD2, vecD3 ], [], 2 );
	vecYHat = vecY/norm(vecY);
	matELo = matIV - prm.moveToELoCoeff*vecYHat*(vecYHat');
	matEHi = matIV - prm.moveToEHiCoeff*vecYHat*(vecYHat');
	matALo_updated = matELo' * ( matALo + diag(vecDLo) ) * matELo;
	matAHi_updated = matEHi' * ( matAHi + diag(vecDHi) ) * matEHi;
	matALo_updated = (matALo_updated'+matALo_updated)/2.0;
	matAHi_updated = (matAHi_updated'+matAHi_updated)/2.0;
	if ( 0 )
	[ matRLo, cholFlag ] = chol(matALo_updated);
	if ( 0 ~= cholFlag || min(diag(matRLo)) <= prm.cholRelTol * max(abs(diag(matRLo))) )
		matALo_updated += max(max(abs(matALo_updated)))*sqrt(eps)*matIV;
	endif
	[ matRHi, cholFlag ] = chol(matAHi_updated);
	if ( 0 ~= cholFlag || min(diag(matRHi)) <= prm.cholRelTol * max(abs(diag(matRHi))) )
		matAHi_updated += max(max(abs(matAHi_updated)))*sqrt(eps)*matIV;
	endif
	elseif (1)
	if ( min(diag(matALo_updated)) < 0.0 )
		matALo_updated += abs(min(diag(matALo_updated)))*(1.0+sqrt(eps))*matIV;
	endif
	if ( min(diag(matAHi_updated)) < 0.0 )
		matAHi_updated += abs(min(diag(matAHi_updated)))*(1.0+sqrt(eps))*matIV;
	endif
	endif
	%
	%
	%
	fModelDat.vecX = vecX + matV*vecY;
	fModelDat.vecF = vecF_next;
	%fModelDat.matV = matV;
	fModelDat.matW = matW_updated;
	fModelDat.matALo = matALo_updated;
	fModelDat.matAHi = matAHi_updated;
	%fModelDat.matB = matB;
	fModelDat.matVLocal = zeros(sizeX,0);
	fModelDat.strState = "moved";
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

function vecU = __applyPrecon( vecRhoF, prm, vecX, vecF )
	% TODO: Add scaling-matrix preconditioning.
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
	% TODO: Update epsFD based on steps taken.
	mydefs;
	v = norm(vecV);
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < v );
	endif
	vecFP = funchF( vecX + prm.epsFD*vecV );
	vecW = ( vecFP - vecF ) / (prm.epsFD * v);
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
	switch ( tolower(prm.curveType) )
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


function [ retCode, fevalIncr, fModelDat ] = __shrinkTR( funchF, vecY, fModelDat, prm )
	[ retCode, fevalIncr, fModelDat ] = __modifyB( funchF, vecY, fModelDat, prm );
	return;
endfunction
function [ retCode, fevalIncr, fModelDat ] = __expandTR( funchF, vecY, fModelDat, prm )
	[ retCode, fevalIncr, fModelDat ] = __modifyB( funchF, vecY, fModelDat, prm );
	return;
endfunction
function [ retCode, fevalIncr, fModelDat ] = __modifyB( funchF, vecY, fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
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
		assert( isrealscalar(prm.initExpandExp) );
		assert( 0.0 < prm.fTol );
		assert( 0.0 < prm.epsFD );
		assert( 0.0 < prm.orthoTol );
		assert( 0.0 < prm.initExpandExp );
		assert( prm.initExpandExp < 1.0 );
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
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	strState = fModelDat.strState; % String indicating "state" of fModelDat.
	sizeX = size(fModelDat.vecX,1);
	sizeF = size(fModelDat.vecF,1);
	sizeV = size(fModelDat.matV,2);
	sizeB = size(fModelDat.matB,1);
	sizeVLocal = size(fModelDat.matVLocal,2);
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__MEDIUM )
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
		switch ( tolower(strState) )
		case { "init" }
			% Okay.
		case { "moved" }
			% Okay.
		otherwise
			error( "Invalid value of strState." );
		endswitch
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
		assert( reldiff(matVLocal'*matVLocal,eye(sizeVLocal,sizeVLocal)) < sqrt(eps) );
		assert( reldiff(matV*(matV'*matVLocal),matVLocal) < sqrt(eps) );
	endif
	%
	%
	%
	return;
endfunction


function __validateStudyDat( fModelDat, studyDat, prm )
	mydefs;
	if ( prm.valdLev < VALDLEV__LOW )
		return;
	endif
	if ( isempty(fModelDat) )
		return;
	endif
	%
	%vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	%matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	%matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%sizeX = size(fModelDat.vecX,1);
	%sizeF = size(fModelDat.vecF,1);
	sizeV = size(fModelDat.matV,2);
	%sizeB = size(fModelDat.matB,1);
	%sizeVLocal = size(fModelDat.matVLocal,2);
	%
	matWTW = studyDat.matWTW;
	matD = studyDat.matD;
	vecWTF = studyDat.vecWTF;
	matBTB = studyDat.matBTB;
	matC = studyDat.matC;
	vecY_ideal = studyDat.vecY_ideal;
	vecY_zeroV = studyDat.vecY_zeroV;
	vecY_loVar = studyDat.vecY_loVar;
	vecY_hiVar = studyDat.vecY_hiVar;
	funchEta_zeroV = studyDat.funchEta_zeroV;
	funchEta_loVar = studyDat.funchEta_loVar;
	funchEta_hiVar = studyDat.funchEta_hiVar;
	%
	omega = sumsq(vecF)/2.0;
	%
	%
	%
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		assert( isrealarray(matWTW,[sizeV,sizeV]) );
		assert( issymmetric(matWTW) );
		assert( min(diag(matWTW)) >= 0.0 );
		assert( isrealarray(matD,[sizeV,sizeV]) );
		assert( isdiag(matD) );
		assert( min(diag(matD)) >= 0.0 );
		assert( reldiff(diag(matD),diag(matWTW)) < sqrt(eps) );
		assert( isrealarray(vecWTF,[sizeV,1]) );
		assert( isrealarray(matBTB,[sizeV,sizeV]) );
		assert( issymmetric(matBTB) );
		assert( min(diag(matBTB)) >= 0.0 );
		assert( isrealarray(matC,[sizeV,sizeV]) );
		assert( issymmetric(matC) );
		assert( min(diag(matC)) > 0.0 );
		%
		if ( ~(norm(matB*vecY_zeroV) <= 1.0 + prm.candStepRelTol) )
			norm(matB*vecY_zeroV) - 1.0
			prm.candStepRelTol
		endif
		assert( norm(matB*vecY_zeroV) <= 1.0 + min([ 0.5, 100.0*prm.candStepRelTol ]) );
		assert( norm(matB*vecY_loVar) <= 1.0 + min([ 0.5, 100.0*prm.candStepRelTol ]) );
		assert( norm(matB*vecY_hiVar) <= 1.0 + min([ 0.5, 100.0*prm.candStepRelTol ]) );
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
	switch ( tolower(prm.curveScaling) )
	case { "b", "btb", "boundary", "optimal" }
	switch ( tolower(prm.curveType) )
	case { "l", "lev", "levenberg" }
		assert( funchEta_zeroV(vecY_ideal) <= (1.0+100.0*sqrt(eps))*abs(funchEta_zeroV(vecY_zeroV)) + eps*100.0*sumsq(vecF) );
		assert( funchEta_zeroV(vecY_zeroV) <= (1.0+100.0*sqrt(eps))*abs(funchEta_zeroV(vecY_loVar)) + eps*100.0*sumsq(vecF) );
		assert( funchEta_zeroV(vecY_loVar) <= (1.0+100.0*sqrt(eps))*abs(funchEta_zeroV(vecY_hiVar)) + eps*100.0*sumsq(vecF) );
		assert( funchEta_zeroV(vecY_hiVar) <= (1.0+100.0*sqrt(eps))*sumsq(vecF)/2.0 );
		assert( funchEta_loVar(vecY_loVar) <= (1.0+100.0*sqrt(eps))*abs(funchEta_loVar(vecY_zeroV)) + eps*100.0*sumsq(vecF) );
		assert( funchEta_loVar(vecY_loVar) <= (1.0+100.0*sqrt(eps))*abs(funchEta_loVar(vecY_hiVar)) + eps*100.0*sumsq(vecF) );
		assert( funchEta_hiVar(vecY_hiVar) <= (1.0+100.0*sqrt(eps))*abs(funchEta_hiVar(vecY_zeroV)) + eps*100.0*sumsq(vecF) );
		assert( funchEta_hiVar(vecY_hiVar) <= (1.0+100.0*sqrt(eps))*abs(funchEta_hiVar(vecY_loVar)) + eps*100.0*sumsq(vecF) );
		%
		% Note that this test is a property of the "optim" curve,
		% not a propery of hitting the boundary exactly.
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
			assert( eta_temp + 1000.0*eps*omega >= eta*(1.0-100.0*sqrt(eps)) || b_temp + 1000.0*eps >= b*(1.0-100.0*sqrt(eps)) );
			eta_temp_etaTempMin = eta_temp;
			b_temp_etaTempMin = b_temp;
			eta_temp_bTempMin = eta_temp;
			b_temp_bTempMin = b_temp;
			for n=1:100
				vecY_temp = vecY + yNorm * sqrt(eps)*(2.0*rand(sizeV,1)-1.0);
				b_temp = norm(matB*vecY_temp);
				eta_temp = funchEta(vecY_temp);
				assert( eta_temp + 1000.0*eps*omega >= eta*(1.0-100.0*sqrt(eps)) || b_temp + 1000.0*eps >= b*(1.0-100.0*sqrt(eps)) );
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
	%
	%
	return;
endfunction
