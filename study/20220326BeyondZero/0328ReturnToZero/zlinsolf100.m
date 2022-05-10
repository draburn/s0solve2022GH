% DRaburn 2022.05.07
%  zlinsolf100
%  First concrete attempt at MVP.
%
% DRaburn 2022.05.08...
%   - Overhaulled trust region / boundary stuff.
%
% DRaburn 2022.05.07...
%   Todo:
%    - Implement partial-quadratic update to A when taking a step.
%    - Refactor trial step handling.
%   Soon:
%    - Implement halding of BLM / overly small TR as momentum / restart jump, consider quadratic terms?
%    - Allow for a constant matrix preconditioner.
%    - Allow dog-leg intstead of Levenberg curve.
%    - Refactor from scratch, test, compare, etc.
%   Much later:
%    - Optimize Octave code.
%    - Code in C/C++.
%    - Try with PIES.
%    - Test against available solvers.
%   Post-MVP:
%    ~ Reconsider anything and everything.
%    - Sepratrix domain enumeration for "basin" BLM handling?
%    - Non-smooth functions: automatically adjust epsFD for differentiation, what else?
%    - Automatic scaling in X, (makes subspace basis non-orthogonal)?

function [ vecX_best, vecF_best, datOut ] = zlinsolf100( funchF, vecX_initial, vecF_initial=[], prm=[] )
	% Init...
	datOut = [];
	fevalCount = 0;
	if (isempty(vecF_initial))
		vecF_initial = funchF( vecX_initial );
		fevalCount++;
	endif
	prm = __initPrm( vecX_initial, vecF_initial, prm );
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, "Initializing subspace." );
	[ fModelDat, datOut_initModel ] = __initModel( funchF, vecX_initial, vecF_initial, prm );
	fevalCount += datOut_initModel.fevalCount;
	%
	% Current local vecX and vecF are stored only in fModelDat to avoid redundancy.
	iterCount = 0;
	vecX_cand = []; % Candidate for next vecX.
	vecF_cand = [];
	vecX_best = vecX_initial;
	vecF_best = vecF_initial;
	%
	stepCount = 0;
	datOut.iterCountOfSteps(stepCount+1) = 1;
	%
	while (1)
		%
		datOut.iterCountVals(iterCount+1) = iterCount;
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.stepCountVals(iterCount+1) = stepCount;
		datOut.fNormVals(iterCount+1) = norm(vecF_best);
		datOut.vecXVals(:,iterCount+1) = fModelDat.vecX;
		datOut.vecFVals(:,iterCount+1) = fModelDat.vecF;
		iterCount++;
		%
		%
		fModelDat = __analyzeModel( fModelDat, prm );
		%
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( ...
		  "  iter:  %5d / %5d;  omega: %10.3e / %10.3e;  sizeV: %d;  sizeVLoc: %d;  sizeB: %d.", ...
		  iterCount, prm.iterMax, sumsq(fModelDat.vecF)/2.0, prm.omegaTol, ...
		  size(fModelDat.matV,2), size(fModelDat.matVLocal,2), size(fModelDat.matB,2) ) );
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( ...
		  "  omega:  %10.3e, %10.3e IU;  %10.3e, %10.3e IB;  %10.3e, %10.3e PB.", ...
		  fModelDat.omegaModelAvgIU, fModelDat.omegaModelVarIU, ...
		  fModelDat.omegaModelAvgIB, fModelDat.omegaModelVarIB, ...
		  fModelDat.omegaModelAvgPB, fModelDat.omegaModelVarPB ) );
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( ...
		  "  steps:  %10.3e, %10.3e IU;  %10.3e, %10.3e IB;  %10.3e, %10.3e PB.", ...
		  norm(fModelDat.vecYIU), fModelDat.bIU, ...
		  norm(fModelDat.vecYIB), fModelDat.bIB, ...
		  norm(fModelDat.vecYPB), fModelDat.bPB ) );
		%
		if ( fModelDat.bPB > 1.0 + sqrt(eps) || fModelDat.bIB > 1.0 + sqrt(eps) )
			msg( __FILE__, __LINE__, "ERROR: fModelDat.bPB or bIB > 1.0 + sqrt(eps)" );
			break;
		endif
		if ( fModelDat.omegaModelAvgIU > fModelDat.omega ...
		  || fModelDat.omegaModelAvgIB > fModelDat.omega ...
		  || fModelDat.omegaModelAvgPB > fModelDat.omega )
			msg( __FILE__, __LINE__, "ERROR: omegaModelAvg (of at least one flavor) > omega." );
			break;
		endif
		%
		% Simple stoping criteria.
		if ( norm(vecF_best) <= prm.fTol )
			msgif( prm.msgMain, __FILE__, __LINE__, "SUCCESS: sumsq(vecF_best) <= prm.fTol^2." );
			break;
		endif
		%
		if ( iterCount > prm.iterMax )
			msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterMax." );
			if (1)
				msgif( prm.msgCopious, __FILE__, __LINE__, "vvvvv Data dump..." );
				__dumpModel( fModelDat, prm );
				%vecX_cand
				%vecF_cand
				msgif( prm.msgCopious, __FILE__, __LINE__, "^^^^^ End data dump." );
				msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterMax." );
			endif
			break;
		endif
		%
		if ( stopsignalpresent() )
			msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		%
		%if ( fModelDat.omegaModelAvgIU > prm.omegaTol )
		if ( fModelDat.omegaModelAvgIU > fModelDat.omega/100.0 )
		if ( size(fModelDat.matV,2) < size(vecX_initial,1) )
			% Conceptually, we could also consider "refreshing" something already in our space;
			%  alternatively, we could be more conservative about expanding the subspace.
			% This is an area for future analysis.
			msgif( prm.msgCopious, __FILE__, __LINE__, "Expanding subspace." );
			[ fModelDat, datOut_expandModel ] = __expandModel( fModelDat.vecFModelIU, funchF, fModelDat, prm );
			fevalCount += datOut_expandModel.fevalCount;
			continue;
		endif
		endif
		%
		if ( fModelDat.bIU <= 1.0 && fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU <= prm.omegaTol )
			msgif( prm.msgCopious, __FILE__, __LINE__, "Trying ideal unbound step." );
			vecX_trial = fModelDat.vecXIU;
			vecF_trial = funchF( vecX_trial );
			omega_trial = sumsq(vecF_trial)/2.0;
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  omega_trial = %g.", omega_trial ) );
			fevalCount++;
			if ( norm(vecF_trial) < norm(vecF_best) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Step is new best." );
				vecX_best = vecX_trial;
				vecF_best = vecF_trial;
			endif
			if ( ~isempty(vecF_cand) )
			if ( norm(vecF_trial) >= norm(vecF_cand) )
				msgif( prm.msgNotice, __FILE__, __LINE__, "Current trial is worse than earlier candidate; forcing acceptance of earlier candidate." );
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_cand)/2.0, fModelDat.omega-sumsq(vecF_cand)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_cand, vecF_cand, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			endif
			%
			avefaThresh = mygetfield( prm, "avefaThresh", 0.5 ); % Actual vs expect fall acceptace threshold
			assert( 0.0 < avefaThresh );
			assert( avefaThresh < 1.0 );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  expected fall = %g.", fModelDat.omega - (fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU) ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  actual fall = %g.", fModelDat.omega - omega_trial ) );
			if ( omega_trial <= fModelDat.omega - avefaThresh*( fModelDat.omega - (fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU) ) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Accepting step." );
				excellentThresh = mygetfield( prm, "excellentThresh", 0.1 );
				if ( norm(vecF_trial-fModelDat.vecFModelIU) <= excellentThresh*norm(fModelDat.vecFModelIU) )
					msg( __FILE__, __LINE__, "  Model was very accurate; removing wall(s)." );
					bRemoveFactor = mygetfield( prm, "bRemoveFactor", 1.5 );
					fModelDat = __removeB( bRemoveFactor*fModelDat.vecYIU, fModelDat, prm );
				endif
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_trial)/2.0, fModelDat.omega-sumsq(vecF_trial)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			msgif( prm.msgCopious, __FILE__, __LINE__, "  Rejecting step." );
			%
			if ( norm(vecF_trial) < norm(fModelDat.vecF) )
				% Trial is better than current, at least, so it's a candidate.
				vecX_cand = vecX_trial;
				vecF_cand = vecF_trial;
			endif
			%
			% Consider adding a new barrier/wall.
			% But, first, make sure uncertainty in W was not the issue.
			vecY = fModelDat.vecYIU;
			vecU = fModelDat.matV*vecY;
			vecV = __calcOrthonorm( vecU, fModelDat.matVLocal, prm );
			if ( 0.0 ~= norm(vecV) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Refreshing trial step before adding wall." );
				[ fModelDat, datOut_refresh ] = __refresh( vecY, funchF, fModelDat, prm );
				fevalCount += datOut_refresh.fevalCount;
			endif
			vecFModel = fModelDat.vecF + fModelDat.matW*vecY;
			badThresh = mygetfield( prm, "badThresh", 0.5 );
			if ( norm( vecF_trial - vecFModel ) > badThresh * norm(vecFModel) )
				msg( __FILE__, __LINE__, "  Model was very bad; adding wall." );
				bAddFactor = mygetfield( prm, "bAddFactor", 0.5 );
				fModelDat = __addB( bAddFactor*vecY, fModelDat, prm );
			endif
			%
			clear vecX_trial;
			clear vecF_trial;
			clear omega_trial;
			continue;
		endif
		%
		minRelFallThresh = mygetfield( prm, "minRelFallThresh", 1.0E-4 );
		if ( fModelDat.omegaModelAvgIB > fModelDat.omega*(1.0-minRelFallThresh) )
			msgif( prm.msgMain, __FILE__, __LINE__, "We seem to have no way to reduce omega much." );
			msgif( prm.msgMain, __FILE__, __LINE__, "  This is expected to happen near a bad local minimum." );
			msgif( prm.msgMain, __FILE__, __LINE__, "  Todo: add handling for this case." );
			break;
		endif
		%
		%
		%
		practicalRelFallThresh = mygetfield( prm, "practicalRelFallThresh", 0.5 );
		if ( fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB <= fModelDat.omega - practicalRelFallThresh*(fModelDat.omega-fModelDat.omegaModelAvgIB) )
			msgif( prm.msgCopious, __FILE__, __LINE__, "Trying practical bound step." );
			vecX_trial = fModelDat.vecXPB;
			vecF_trial = funchF( vecX_trial );
			omega_trial = sumsq(vecF_trial)/2.0;
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  omega_trial = %g.", omega_trial ) );
			fevalCount++;
			if ( norm(vecF_trial) < norm(vecF_best) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Step is new best." );
				vecX_best = vecX_trial;
				vecF_best = vecF_trial;
			endif
			if ( ~isempty(vecF_cand) )
			if ( norm(vecF_trial) >= norm(vecF_cand) )
				msgif( prm.msgNotice, __FILE__, __LINE__, "Current trial is worse than earlier candidate; forcing acceptance of earlier candidate." );
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_cand)/2.0, fModelDat.omega-sumsq(vecF_cand)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_cand, vecF_cand, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			endif
			%
			avefaThresh = mygetfield( prm, "avefaThresh", 0.5 ); % Actual vs expect fall acceptace threshold
			assert( 0.0 < avefaThresh );
			assert( avefaThresh < 1.0 );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  expected fall = %g.", fModelDat.omega - (fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB) ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  actual fall = %g.", fModelDat.omega - omega_trial ) );
			if ( omega_trial <= fModelDat.omega - avefaThresh*( fModelDat.omega - (fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB) ) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Accepting step." );
				excellentThresh = mygetfield( prm, "excellentThresh", 0.1 );
				if ( norm(vecF_trial-fModelDat.vecFModelPB) <= excellentThresh*norm(fModelDat.vecFModelPB) )
					msg( __FILE__, __LINE__, "  Model was very accurate; removing wall(s)." );
					bRemoveFactor = mygetfield( prm, "bRemoveFactor", 1.5 );
					fModelDat = __removeB( bRemoveFactor*fModelDat.vecYPB, fModelDat, prm );
				endif
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_trial)/2.0, fModelDat.omega-sumsq(vecF_trial)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			msgif( prm.msgCopious, __FILE__, __LINE__, "  Rejecting step." );
			%
			if ( norm(vecF_trial) < norm(fModelDat.vecF) )
				% Trial is better than current, at least, so it's a candidate.
				vecX_cand = vecX_trial;
				vecF_cand = vecF_trial;
			endif
			%
			%
			% Consider adding a new barrier/wall.
			% But, first, make sure uncertainty in W was not the issue.
			vecY = fModelDat.vecYPB;
			vecU = fModelDat.matV*vecY;
			vecV = __calcOrthonorm( vecU, fModelDat.matVLocal, prm );
			if ( 0.0 ~= norm(vecV) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Refreshing trial step before adding wall." );
				[ fModelDat, datOut_refresh ] = __refresh( vecY, funchF, fModelDat, prm );
				fevalCount += datOut_refresh.fevalCount;
			endif
			vecFModel = fModelDat.vecF + fModelDat.matW*vecY;
			badThresh = mygetfield( prm, "badThresh", 0.5 );
			if ( norm( vecF_trial - vecFModel ) > badThresh * norm(vecFModel) )
				msg( __FILE__, __LINE__, "  Model was very bad; adding wall." );
				bAddFactor = mygetfield( prm, "bAddFactor", 0.5 );
				fModelDat = __addB( bAddFactor*vecY, fModelDat, prm );
			endif
			%
			clear vecX_trial;
			clear vecF_trial;
			clear omega_trial;
			continue;
		endif
		%
		msgif( prm.msgCopious, __FILE__, __LINE__, "Refreshing subspace." );
		[ fModelDat, datOut_refresh ] = __refresh( fModelDat.vecYPB, funchF, fModelDat, prm );
		fevalCount += datOut_refresh.fevalCount;
		continue;
	endwhile
	%
	%
	datOut.iterCountVals(iterCount+1) = iterCount;
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.fNormVals(iterCount+1) = norm(vecF_best);
	datOut.vecXVals(:,iterCount+1) = fModelDat.vecX;
	datOut.vecFVals(:,iterCount+1) = fModelDat.vecF;
	datOut.iterCountOfSteps(stepCount+1) = iterCount;
	%
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
	datOut.stepCount = stepCount;
return;
endfunction


function prm = __initPrm( vecX, vecF, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	prm.msgCopious = mygetfield( prm, "msgCopious", verbLev >= VERBLEV__COPIOUS );
	prm.msgProgress = mygetfield( prm, "msgProgress", verbLev >= VERBLEV__PROGRESS );
	prm.msgMain = mygetfield( prm, "msgMain", verbLev >= VERBLEV__MAIN );
	prm.msgNotice = mygetfield( prm, "msgNotice", verbLev >= VERBLEV__NOTICE );
	prm.debugMode = mygetfield( prm, "debugMode", true );
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	%fTol = max([ sqrt(eps)*norm(vecF), eps ]);
	fTol = eps;
	prm.fTol = mygetfield( prm, "fTol", fTol );
	prm.iterMax = mygetfield( prm, "iterMax", ceil(20 + 10*sqrt(sizeX) + 0.01*sizeX) );
	%
	prm.omegaTol = fTol^2/2.0;
return;
endfunction


function [ fModelDat, datOut ] = __initModel( funchF, vecX, vecF, prm )
	datOut = [];
	fevalCount = 0;
	%
	fModelDat = mygetfield( prm, "fModelDat_initial", [] );
	if (~isempty(fModelDat))
		assert( ~isempty(fModelDat,"vecX",[]) );
		assert( ~isempty(fModelDat,"vecF",[]) );
		assert( reldiff(fModelDat.vecX,vecX,eps^2) <= eps );
		if ( isempty(vecF) )
			vecF = fModelDat.vecF;
		else
			assert( reldiff(fModelDat.vecF,vecF,eps^2) <= eps );
		endif
		%
		% Do more checks here?
		%
		return;
	endif
	%
	if (isempty(vecF))
		vecF = funchF(vecX);
		fevalCount++;
	endif
	if (prm.debugMode)
		sizeX = size(vecX,1);
		sizeF = size(vecF,1);
		assert( sizeX >= 1 );
		assert( sizeF >= 1 );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
	endif
	%
	%
	fNorm = norm(vecF);
	if ( 0.0 == fNorm )
		error( "Initial vecF is zero." );
	endif
	vecV = vecF/fNorm;
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecW) )
		error( "Initial vecW is zero." );
	endif
	%
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matVLocal = [ vecV ];
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matA = [ 0.0 ];
	fModelDat.matB = [];
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function vecW = __calcJV( vecV, funchF, vecX, vecF, prm )
	v = norm(vecV);
	assert( 0.0 < v );
	epsFD = mygetfield( prm, "epsFD", eps^0.3 );
	assert( 0.0 < epsFD );
	vecFP = funchF( vecX + epsFD*vecV );
	vecW = ( vecFP - vecF ) / (epsFD*norm(vecV));
return;
endfunction


function fModelDat = __analyzeModel( fModelDat, prm )
	% Unpack.
	%matVLocal = fModelDat.matVLocal;
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	% Basics.
	matWTW = matW'*matW;
	vecMG = -(matW'*vecF);
	matSCurve = eye(sizeV,sizeV);
	%
	vecYIU = __calcStep( matWTW, vecMG, prm );
	vecYIB = __calcBoundStep( matWTW, vecMG, matB, matSCurve, prm );
	vecYPB = __calcBoundStep( matWTW + matA, vecMG, matB, matSCurve, prm );
	%
	vecFModelIU = vecF + (matW*vecYIU);
	vecFModelIB = vecF + (matW*vecYIB);
	vecFModelPB = vecF + (matW*vecYPB);
	%
	fModelDat.vecYIU = vecYIU;
	fModelDat.vecYIB = vecYIB;
	fModelDat.vecYPB = vecYPB;
	%
	fModelDat.vecXIU = vecX + (matV*vecYIU);
	fModelDat.vecXIB = vecX + (matV*vecYIB);
	fModelDat.vecXPB = vecX + (matV*vecYPB);
	%
	fModelDat.vecFModelIU = vecFModelIU;
	fModelDat.vecFModelIB = vecFModelIB;
	fModelDat.vecFModelPB = vecFModelPB;
	%
	fModelDat.omegaModelAvgIU = sumsq(vecFModelIU)/2.0;
	fModelDat.omegaModelAvgIB = sumsq(vecFModelIB)/2.0;
	fModelDat.omegaModelAvgPB = sumsq(vecFModelPB)/2.0;
	%
	fModelDat.omegaModelVarIU = vecYIU'*matA*vecYIU;
	fModelDat.omegaModelVarIB = vecYIB'*matA*vecYIB;
	fModelDat.omegaModelVarPB = vecYPB'*matA*vecYPB;
	%
	if ( isempty(matB) )
		fModelDat.bIU = 0.0;
		fModelDat.bIB = 0.0;
		fModelDat.bPB = 0.0;
	else
		fModelDat.bIU = max(abs(vecYIU'*matB));
		fModelDat.bIB = max(abs(vecYIB'*matB));
		fModelDat.bPB = max(abs(vecYPB'*matB));
	endif
	%
	fModelDat.omega = sumsq(vecF)/2.0;
	%
	assert( fModelDat.omegaModelAvgIU <= fModelDat.omega );
	assert( fModelDat.omegaModelAvgIB <= fModelDat.omega );
	assert( fModelDat.omegaModelAvgPB <= fModelDat.omega );
return;
endfunction


function vecY = __calcStep( matH, vecMG, prm )
	[ matR, cholFlag ] = chol( matH );
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	if ( 0 == cholFlag && min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		vecY = matR \ (matR'\vecMG);
	else
		msgif( prm.msgNotice, __FILE__, __LINE__, "Extrapolating step to singular point. (Perhaps should bail instead?)" );
		hScale = max(max(abs(matH)));
		assert( 0.0 < hScale );
		sizeV = size(matH,1);
		matIV = eye(sizeV,sizeV);
		epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
		matH1 = matH + ( 1.0 * epsRelRegu * hScale ) * matIV;
		matH2 = matH + ( 2.0 * epsRelRegu * hScale ) * matIV;
		matR1 = chol( matH1 );
		matR2 = chol( matH2 );
		vecY1 = matR1 \ (matR1'\vecMG);
		vecY2 = matR2 \ (matR2'\vecMG);
		vecY = (2.0*vecY1) - vecY2;
	endif
return;
endfunction


function vecY = __calcBoundStep( matH, vecMG, matB, matSCurve, prm );
	if ( isempty(matB) || 0.0==max(max(abs(matB))) )
		vecY = __calcStep( matH, vecMG, prm );
		return;
	endif
	%
	%cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	%
	%[ matR, cholFlag ] = chol( matH );
	%if ( 0 == cholFlag && min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
	% Looks like chol isn't accurate enough, and/or "\" triggers a check based on rcond()?
	% So, we'l use rcond too.
	rc = rcond( matH );
	if ( rcond(matH) > sqrt(eps) )
		s1 = 1.0;
	else
		epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
		hScale = max(max(abs(matH)));
		sScale = max(max(abs(matSCurve)));
		assert( 0.0 < hScale );
		assert( 0.0 < sScale );
		s1 = 1.0 - (epsRelRegu*hScale/(epsRelRegu*hScale+sScale));
	endif
	assert( s1 >= 0.0 );
	%
	funchYOfS = @(s)( ( s*matH + (1.0-s)*matSCurve ) \ (s*vecMG) );
	funchBOfY = @(y)( max(abs(y'*matB)) );
	vecY1 = funchYOfS(s1);
	b1 = funchBOfY(vecY1);
	if ( b1 <= 1.0 )
		vecY = vecY1;
		return;
	endif
	%
	funchBM1OfS = @(s)( funchBOfY(funchYOfS(s)) - 1.0 );
	%
	s = fzero( funchBM1OfS, [0.0, s1] );
	vecY = funchYOfS(s);
return;
endfunction


function [ fModelDat, datOut ] = __expandModel( vecU, funchF, fModelDat, prm )
	fevalCount = 0;
	% Unpack.
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,2);
	%
	%
	uNorm = norm(vecU);
	assert( 0.0 < uNorm );
	vecV = __calcOrthonorm( vecU, matV, prm );
	if ( 0.0 == norm(vecV) )
		error( "Failed to generate new subspace basis vector." );
	endif
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecV) )
		error( "Jacobian along new subspace basis vector is zero." );
	endif
	%
	%
	fModelDat.matVLocal = [ matVLocal, vecV ];
	fModelDat.matV = [ matV, vecV ];
	fModelDat.matW = [ matW, vecW ];
	fModelDat.matA = zeros(sizeV+1,sizeV+1);
	fModelDat.matA(1:sizeV,1:sizeV) = matA;
	fModelDat.matB = [ matB; zeros(1,sizeB) ];
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function vecV = __calcOrthonorm( vecU, matV, prm )
	numPasses = 2;
	u0 = norm(vecU);
	if (0.0==u0)
		return;
	endif
	if (isempty(matV))
		vecV = vecU/u0;
		return;
	endif
	orthoTol = mygetfield( prm, "orthoTol", 1.0e-10 );
	vecV = vecU;
	for n=1:numPasses
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( v <= orthoTol*u0 )
			vecV(:) = 0.0;
			return;
		else
			vecV /= v;
		endif
	endfor
return;
endfunction


function fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm )
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	%matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( " ( ||F||: %10.3e -> %10.3e. )", norm(vecF), norm(vecF_trial) ) );
	%
	%
	vecDeltaX = vecX_trial - vecX;
	vecY = matV'*vecDeltaX;
	yNorm = norm(vecY);
	assert( 0.0 < yNorm );
	%
	%
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, "TODO: consider partial quadratic update (OSQU)." );
	msgif( prm.msgCopious, __FILE__, __LINE__, "  Here, only linear (Broyden) update is applied." );
	vecFModel_trial = vecF + matW*vecY;
	vecRho = vecF_trial - vecFModel_trial;
	matW_plus = matW + vecRho * (vecY')/(yNorm^2);
	%
	%
	%
	vecYHat = vecY/yNorm;
	stepUpdateAccuracyCoeff = mygetfield( prm, "stepUpdateAccuracyCoeff", 0.0 );
	assert( 0.0 <= stepUpdateAccuracyCoeff );
	assert( stepUpdateAccuracyCoeff <= 1.0 );
	matE = eye(sizeV,sizeV) - (stepUpdateAccuracyCoeff*vecYHat)*(vecYHat');
	%
	matD = diag(max([ abs(diag(matW'*matW)), abs(diag(matW_plus'*matW_plus)) ]'));
	foo1 = sumsq( matW_plus*vecY - matW*vecY ) - vecY'*matA*vecY;
	foo2 = vecY'*matD*vecY;
	if ( foo1 <= 0.0 )
		s = 0.0;
	elseif ( foo2 <= foo1 )
		s = 1.0;
	else
		s = foo1 / foo2;
		assert( 0.0 <= s );
		assert( s <= 1.0 );
	endif
	coeffD = mygetfield( prm, "coeffD", 1.0 );
	s*=coeffD;
	matA = matE*( matA + s*matD )*matE;
	%
	%
	%
	fModelDat.matVLocal = [];
	fModelDat.vecX = vecX_trial;
	fModelDat.vecF = vecF_trial;
	fModelDat.matW = matW;
	fModelDat.matA = matA;
return;
endfunction


function [ fModelDat, datOut ] = __refresh( vecY, funchF, fModelDat, prm )
	fevalCount = 0;
	% Unpack.
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	%matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	%
	yNorm = norm(vecY);
	assert( 0.0 < yNorm );
	vecYHat = vecY/yNorm;
	vecU = matV*vecYHat;
	vecV = __calcOrthonorm( vecU, matVLocal, prm );
	if ( 0.0 == norm(vecV) )
		error( "Failed to generate local subspace basis vector." );
	endif
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecW) )
		error( "Jacobian along local subspace basis vector is zero." );
	endif
	%
	matVLocal = [ matVLocal, vecV ];
	%
	matE = eye(sizeV,sizeV) - matV'*matVLocal*(matVLocal')*matV;
	%
	fModelDat.matVLocal = matVLocal;
	fModelDat.matW = matW + (matW*vecYHat-vecW)*(vecYHat');
	fModelDat.matA = matE * matA * matE;
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function fModelDat = __addB( vecY, fModelDat, prm )
	% Unpack.
	%vecX = fModelDat.vecX;
	%vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%matV = fModelDat.matV; % Subspace basis matrix.
	%matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	%matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	%sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	ysumsq = sumsq(vecY);
	assert( 0.0 < ysumsq );
	if (isempty(matB))
		fModelDat.matB = vecY/ysumsq;
	else
		fModelDat.matB = [ matB, vecY/ysumsq ];
	endif
return;
endfunction



function fModelDat = __removeB( vecY, fModelDat, prm )
	% Unpack.
	%vecX = fModelDat.vecX;
	%vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%matV = fModelDat.matV; % Subspace basis matrix.
	%matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	%matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	%sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	if ( isempty(matB) )
		return;
	endif
	msk = logical( abs(vecY'*matB) >= 1.0 );
	fModelDat.matB = matB( :, msk );
return;
endfunction

function __dumpModel( fModelDat, prm )
	msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvv Begin __dumpModel()..." );
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,2);
	%
	echo__omega = omega
	sizeVLocal = size(matVLocal,2)
	sizeV = size(matV,2)
	%
	vecYIU = fModelDat.vecYIU;
	vecYIB = fModelDat.vecYIB;
	vecYPB = fModelDat.vecYPB;
	%
	vecFModelIU = vecF + (matW*vecYIU);
	vecFModelIB = vecF + (matW*vecYIB);
	vecFModelPB = vecF + (matW*vecYPB);
	%
	omegaModelAvgIU = sumsq(vecFModelIU)/2.0
	omegaModelAvgIB = sumsq(vecFModelIB)/2.0
	omegaModelAvgPB = sumsq(vecFModelPB)/2.0
	%
	omegaModelVarIU = vecYIU'*matA*vecYIU
	omegaModelVarIB = vecYIB'*matA*vecYIB
	omegaModelVarPB = vecYPB'*matA*vecYPB
	%
	if (isempty(matB))
		bIU = 0.0
		bIB = 0.0
		bPB = 0.0
	else
		bIU = max(abs(vecYIU'*matB))
		bIB = max(abs(vecYIB'*matB))
		bPB = max(abs(vecYPB'*matB))
	endif
	%
	msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^ End __dumpModel()." );
return;
endfunction
