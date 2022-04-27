function [ vecXF, vecFF, datOut ] = slinsolf( funchF, vecX0, vecF0, prm, datIn )
	% Parse input.
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	fevalCount = 0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	stepConstraintDat = mygetfield( datIn, "stepConstraintDat", [] );
	preconDat = mygetfield( datIn, "preconDat", [] );
	trialActionDat = mygetfield( datIn, "trialActionDat", [] );
	if (~isempty(trialActionDat))
		msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Input trialActionDat is not empty. Are you sure you know what you're doing?" );
	endif
	%
	localModelDat = mygetfield( datIn, "localModelDat", [] );
	if (isempty(localModelDat))
		localModelDat.vecX0 = vecX0;
		localModelDat.vecF0 = vecF0;
	else
		msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Input localModelDat is not empty. Are you sure you know what you're doing?" );
		assert( reldiff(vecX0,localModelDat.vecX0,realmin) <= eps )
		assert( reldiff(vecF0,localModelDat.vecF0,realmin) <= eps )
	endif
	%
	%
	% Prep result.
	vecX_best = [];
	vecF_best = [];
	stepConstraintDat_best = []; % Essentially, SCD *implied by* best.
	%
	% Prep loop.
	vecX = vecX0;
	vecF = vecF0;
	iterCount = 0;
	while (1)
		iterMax = mygetfield( prm, "iterMax", 10+sizeX );
		if ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		iterCount++;
		%
		%
		% Calculate trial step.
		vecX_trial = __calcDelta( localModelDat, stepConstraintDat, prm );
		vecFModel_trial = __calcFModel( vecX_trial, localModelDat, prm );
		%
		%
		% Decide what to do.
		trialActionDat.vecX_best = vecX_best;
		trialActionDat.vecF_best = vecF_best;
		trialAction = __determineTrialAction( vecX_trial, vecFModel_trial, localModelDat, trialActionDat, prm );
		msg( __FILE__, __LINE__, sprintf( "trialAction = '%s'.", trialAction ) );
		switch (trialAction)
		case "giveUp"
			break;
		case "expandSubspace"
			[ localModelDat_new, ess_datOut ]  = __expandSubspace( funchF, localModelDat, preconDat, prm );
			fevalCount += ess_datOut.fevalCount;
			if ( ~isempty(localModelDat_new) )
				% Retry, from the top, with new model.
				localModelDat = localModelDat_new;
				continue;
			else
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: __expandSubspace() failed." );
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "  Forcing findGoodStep." );
				trialAction = "findGoodStep";
				% Go outside switch.
			endif
		case "findGoodStep"
			% Go below.
		case "tryStep"
			% Go below.
		otherwise
			error( "Invalid value of trialAction." );
		endswitch
		%
		%
		%
		if ( strcmp(trialAction,"findGoodStep") )
			% This is hack-ish, to allow this code to match conventional behavior, without in-loop step search.
			[ vecX_trial, vecF_trial, stepConstraintDat_trial, fgs_datOut ] = __findGoodStep( funchF, vecX0, vecF0, localModelDat, stepConstraintDat, prm );
			fevalCount += fgs_datOut.fevalCount;
			if (isempty(vecX_trial))
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: __findGoodStep() returned []." );
				break;
			endif
			%
			stepIsBest = false;
			if ( isempty(vecX_best) )
				stepIsBest = true;
			elseif ( norm(vecF_step) < norm(vecF_best) )
				stepIsBest = true;
			endif
			if ( stepIsBest )
				vecX_best = vecX_trial;
				vecF_best = vecF_trial;
				stepConstraintDat_best = stepConstraintDat_trial;
			endif
			stepConstrantDat = stepConstraintDat_trial; % Not that this assignment matters.
			%
			break;
		endif
		%
		%
		%
		vecF_trial = funchF( vecX_trial ); fevalCount++;
		trialResult = __determineTrialResult( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, prm );
		msg( __FILE__, __LINE__, sprintf( "||F||: %g, %g, %g.", norm(vecF0), norm(vecFModel_trial), norm(vecF_trial) ) );
		msg( __FILE__, __LINE__, sprintf( "trialResult = '%s'.", trialResult ) );
		%
		switch (trialResult)
		case "excellent"
			% Update (loosen) SCD, then accept step.
			stepConstraintDat = __updateSCD_excellent( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			stepConstraintDat_best = stepConstraintDat;
			break;
		case "good"
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			stepConstraintDat_best = stepConstraintDat;
			break;
		case "okay"
			if ( ~isempty(vecX_best) )
			if ( norm(vecF_best) < norm(vecF_trial) )
				% Previous was better.
				% Update SCD and return pre-existing best.
				stepConstraintDat = __updateSCD_okay( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
				stepConstraintDat_best = stepConstraintDat;
				break;
			endif
			endif
			% Update best *then* update SCD.
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			stepConstraintDat_best = stepConstraintDat;
			stepConstraintDat = __updateSCD_okay( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			continue;
		case "bad"
			if ( ~isempty(vecX_best) )
				msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Trial result was 'bad' even though a previous result was at least 'okay'." );
			endif
			stepConstraintDat = __updateSCD_bad( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			continue;
		case "horrid"
			if ( ~isempty(vecX_best) )
				msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Trial result was 'horrid' even though a previous result was at least 'okay'." );
			endif
			stepConstraintDat = __updateSCD_horrid( vecX_trial, localModelDat, stepConstraintDat, prm );
			continue;
		otherwise
			error( "Invalid value of trialResult." );
		endswitch
		error( "Inaccessible code." );
	endwhile
	%
	vecXF = vecX_best;
	vecFF = vecF_best;
	vecFModelF = __calcFModel( vecX_best, localModelDat, prm );
	datOut.stepConstraintDat = stepConstraintDat;
	datOut.preconDat = preconDat;
	datOut.localModelDat = localModelDat; % Let someone externally updated the data to new point.
	datOut.vecFModelF = vecFModelF;
	datOut.fevalCount = fevalCount;
return;
endfunction



function vecX_trial = __calcDelta( localModelDat, stepConstraintDat, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		% We don't have a valid model (yet).
		vecX_trial = [];
		return;
	endif
	matW = localModelDat.matW;
	sizeV = size(matV,2);
	assert( 0 < sizeV );
	assert( isrealarray(matV,[sizeX,sizeV]) );
	assert( isrealarray(matW,[sizeF,sizeV]) );
	assert( reldiff(matV'*matV,eye(size(matV,2)),eps) <= sqrt(eps) );
	%
	vecG = localModelDat.vecG;
	matH = localModelDat.matH;
	%
	%
	if ( mygetfield( prm, "useMarquardtScaling", false ) )
		hScale = sqrt(sum(sum(matH.^2))/sizeV);
		matS_curve = diag( 1.0 ./ ( sqrt(abs(diag(matH))) + eps*hScale ) ); % Needs testing.
		vecG_curve = matS_curve'*vecG;
		matH_curve = calcHRegu( matS_curve'*matH*matS_curve );
		%
		pCauchy = calcLinishRootOfQuad( 0.5*(vecG_curve'*matH_curve*vecG_curve), -sumsq(vecG_curve), sumsq(vecF0)/2.0 );
		assert( pCauchy >= 0.0 );
		vecDeltaCauchy = pCauchy*matS_curve*(-vecG_curve); % Needs testing.
		vecDeltaNewton = matS_curve*(matH_curve\(-vecG_curve));
	else
		vecG_curve = vecG;
		matH_curve = calcHRegu( matH );
		pCauchy = calcLinishRootOfQuad( 0.5*(vecG_curve'*matH_curve*vecG_curve), -sumsq(vecG_curve), sumsq(vecF0)/2.0 );
		assert( pCauchy >= 0.0 );
		vecDeltaCauchy = pCauchy*(-vecG_curve); % Needs testing.
		vecDeltaNewton = matH_curve\(-vecG_curve);
	endif
	%
	%
	% Crude placeholder...
	if (~isempty(stepConstraintDat))
		error( "Step constraints are not yet supported!" );
	endif
	vecX_trial = vecX0 + matV*vecDeltaNewton;
	%vecX_trial = vecX0 + 0.01*matV*vecDeltaCauchy; % HACK.
return;
endfunction



function trialAction = __determineTrialAction( vecX_trial, vecFModel_trial, localModelDat, trialActionDat, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%
	if (isempty(vecX_trial))
		% We don't have a valid trial.
		trialAction = "expandSubspace";
		return;
	endif
	%
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( isrealarray(vecX_trial,[sizeX,1]) );
	assert( isrealarray(vecFModel_trial,[sizeF,1]) );
	%
	%vecX_best = mygetfield( trialActionDat, "vecX_best", [] );
	%if (~isempty(vecX_best))
	%	vecF_best = trialActionDat.vecF_best;
	%	assert( isrealarray(vecX_best,[sizeX,1]) );
	%	assert( isrealarray(vecF_best,[sizeF,1]) );
	%endif
	%
	%
	% Crude placeholder...
	dta_c0 = mygetfield( prm, "dta_c0", 0.1 );
	if ( norm(vecFModel_trial) < dta_c0 * norm(vecF0) )
		%trialAction = "tryStep";
		trialAction = "findGoodStep"; % This makes us stop expanding subspace, per traditional codes.
		return;
	endif
	trialAction = "expandSubspace";
	return;
return;
endfunction



function [ localModelDat, ess_datOut ]  = __expandSubspace( funchF, localModelDat, preconDat, prm )
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%
	sizeX = size(localModelDat.vecX0,1);
	sizeF = size(localModelDat.vecF0,1);
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	matV = mygetfield( localModelDat, "matV", zeros(sizeX,0) );
	matW = mygetfield( localModelDat, "matW", zeros(sizeF,0) );
	matPhi = mygetfield( localModelDat, "matPhi", zeros(sizeX,0) ); % Make sizeV x sizeV?
	matGamma = mygetfield( localModelDat, "matGamma", zeros(sizeF,0) );
	%
	if ( 0 == size(matW,2) )
		vecU = localModelDat.vecF0;
	else
		vecU = localModelDat.matW(:,end);
	endif
	%
	vecV = __calcGuessV( vecU, preconDat, prm );
	vecV = __calcOrthonorm( vecV, matV, prm );
	if ( norm(vecV) < eps )
		msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Calculation of new basis vector failed." );
		localModelDat = [];
		ess_datOut.fevalCount = fevalCount;
		return;
	endif
	epsFD = mygetfield( prm, "epsFD", eps^0.3 );
	vecFP = funchF( vecX0 + epsFD*vecV ); fevalCount++;
	vecW = ( vecFP - vecF0 ) / (epsFD*norm(vecV));
	matV = [ matV, vecV ];
	matW = [ matW, vecW ];
	sizeV = size(matV,2);
	%
	assert( reldiff(matV'*matV,eye(size(matV,2)),eps) <= sqrt(eps) );
	%
	[ matPhi, matGamma, pp_datOut ] = __phiPatch( funchF, vecX0, vecF0, matV, matW, matPhi, matGamma, prm );
	fevalCount += pp_datOut.fevalCount;
	sizePhi = size(matPhi,2);
	if ( 0 < sizePhi )
		assert( isrealarray(matPhi,[sizeX,sizePhi]) );
		assert( isrealarray(matGamma,[sizeF,sizePhi]) );
		assert( reldiff(matPhi'*matPhi,eye(size(matPhi,2)),eps) <= sqrt(eps) );
	endif
	matVTPhi = matV'*matPhi;
	%
	vecG = matW'*vecF0;
	matWTW = matW'*matW;
	matH = matWTW;
	for n=1:sizePhi
		matH += (vecF0'*matGamma(:,n)) * matVTPhi(:,n) * (matVTPhi(:,n)');
	endfor
	%
	localModelDat.matV = matV;
	localModelDat.matW = matW;
	localModelDat.matWTW = matWTW;
	localModelDat.matPhi = matPhi;
	localModelDat.matGamma = matGamma;
	localModelDat.matVTPhi = matVTPhi;
	localModelDat.vecG = vecG;
	localModelDat.matH = matH;
	ess_datOut.fevalCount = fevalCount;
return;
endfunction
%
function vecV = __calcGuessV( vecU, preconDat, prm )
	% Support conventional AP.
	% Use LU instead?
	matA = mygetfield( preconDat, "matA", [] );
	if (isempty(matA))
		vecV = vecU;
	else
		vecV = matA\vecU;
	endif
return;
endfunction
%
function vecV = __calcOrthonorm( vecV, matV, prm )
	numPasses = 2;
	v = norm(vecV);
	if (0.0==v)
		return;
	endif
	orthoTol = mygetfield( prm, "orthoTol", 1.0e-10 );
	for n=1:numPasses
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( v <= orthoTol )
			vecV(:) = 0.0;
			return;
		else
			vecV /= v;
		endif
	endfor
return;
endfunction
%
function [ matPhi, matGamma, pp_datOut ] = __phiPatch( funchF, vecX0, vecF0, matV, matW, matPhi, matGamma, prm )
	fevalCount = 0;
	sizeF = size(vecF0,1);
	sizeV = size(matV,2);
	matWTW = matW'*matW;
	[ matPsi, matLambda ] = eig( matWTW );
	assert( isrealarray(matPsi,[sizeV,sizeV]) );
	%
	% Re-patch everything for now.
	%
	% Patch every eigenvalue for now.
	matPhi = matV*matPsi;
	sizePhi = sizeV;
	%
	% Patch using second order for now.
	epsFD = mygetfield( prm, "epsFD", eps^0.3 );
	matGamma = zeros(sizeF,sizePhi);
	for n=1:sizePhi
		vecFP = funchF( vecX0 + epsFD * matPhi(:,n) ); fevalCount++;
		vecFM = funchF( vecX0 - epsFD * matPhi(:,n) ); fevalCount++;
		matGamma(:,n) = ( vecFP + vecFM - 2.0*vecF0 ) / (epsFD^2);
	endfor
	%
	pp_datOut.fevalCount = fevalCount;
return;
endfunction



function [ vecX_next, vecF_next, stepConstraintDat_next, fgs_datOut ] = __findGoodStep( funchF, vecX0, vecF0, localModelDat, stepConstraintDat, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	fevalCount = 0;
	msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Performing __findGoodStep()." );
	if (~isempty(stepConstraintDat))
		error( "Step constraints are not yet supported!" );
	endif
	%
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	matV = localModelDat.matV;
	matH = localModelDat.matH;
	vecG = localModelDat.vecG;
	sizeV = size(matV,2);
	assert( 0 < sizeV );
	%
	%
	matIV = eye(sizeV,sizeV);
	if ( mygetfield( prm, "useMarquardtScaling", false ) )
		hScale = sqrt(sum(sum(matH.^2))/sizeV);
		matS_curve = diag( 1.0 ./ ( sqrt(abs(diag(matH))) + eps*hScale ) ); % Needs testing.
		vecG_curve = matS_curve'*vecG;
		matH_curve = calcHRegu( matS_curve'*matH*matS_curve );
		funchYOfP = @(p)(matS_curve*( (p*matH_curve + (1.0-p)*matIV) \ (-p*vecG_curve) ));
	else
		matS_curve = eye(sizeV,sizeV);
		vecG_curve = vecG;
		matH_curve = calcHRegu( matH );
		funchYOfP = @(p)(( (p*matH_curve + (1.0-p)*matIV) \ (-p*vecG_curve) ));
	endif
	%
	dTreg = mygetfield( stepConstraintDat, "dTreg", [] ); % Is this right?
	if (isempty(dTreg))
		pMax = 1.0;
	else
		pMax = __findPOfDeltaNorm( dTreg, funchDeltaOfP  );
	endif
	vecY_pMax = funchYOfP( pMax );
	vecX_pMax = vecX0 + matV*vecY_pMax;
	vecFModel_pMax = __calcFModel( vecX_pMax, localModelDat, prm );
	%
	btMax = mygetfield( prm, "btMax", 30 );
	deltaNormTol = mygetfield( prm, "deltaNormTol", 100.0*eps );
	btCount = 0;
	vecDelta_rejected = [];
	p = pMax;
	while (1)
		vecY = funchYOfP( p );
		vecDelta = matV*vecY;
		vecX_next = vecX0 + vecDelta;
		vecF_next = funchF( vecX_next );
		fevalCount++;
		%
		if ( norm(vecF_next) < 0.5*norm(vecF0) + 0.5*norm(vecFModel_pMax) )
			break;
		endif
		if ( norm(vecDelta) < deltaNormTol )
			msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, "  step: IMPOSED STOP: norm(vecDelta) < deltaNormTol." );
			break;
		endif
		btCount++;
		if ( btCount > btMax )
			msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, "  step: IMPOSED STOP: btCount > btMax." );
			break;
		endif
		p /= 2.0;
		vecDelta_rejected = vecDelta;
		continue;
	endwhile
	if ( ~isempty(vecDelta_rejected) )
		dTreg = min([ dTreg, norm(vecDelta_rejected) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  step: Have a rejected step. Set dTreg = %0.3e.", dTreg ) )
	endif
	%
	%
	vecFModel_next = __calcFModel( vecX_next, localModelDat, prm );
	rhoThresh0 = mygetfield( prm, "rhoThresh0", 0.05 );
	rhoThresh1 = mygetfield( prm, "rhoThresh1", 0.30 );
	rho = norm(vecF_next-vecFModel_next)/norm(vecF0);
	if ( rho < rhoThresh0 )
		% Model is very accurate at the point.
		dTreg = max([ dTreg, 2.0*norm(vecDelta) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  step: Model was very accurate. Set dTreg = %0.3e.", dTreg ) )
	elseif ( rho > rhoThresh1 )
		% Model was inaccurate at the point.
		dTreg = min([ dTreg, norm(vecDelta) ]); % "min([ dTreg," should be superfluous.
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  step: Model was very inaccurate. Set dTreg = %0.3e.", dTreg ) )
	endif
	%
	stepConstraintDat_next.dTreg = dTreg;
	fgs_datOut.fevalCount = fevalCount;
return;
endfunction



function trialResult = __determineTrialResult( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, prm )
	% Crude placeholder...
	vecF0 = localModelDat.vecF0;
	dtr_c0 = 0.9;
	if ( norm(vecF_trial) < dtr_c0*norm(vecF0) + (1.0-dtr_c0)*norm(vecFModel_trial) )
		trialResult = "good";
		return;
	endif
	trialResult = "bad";
	return;
return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM ACCESS FUNCTIONS
%


function vecFModel_trial = __calcFModel( vecX_trial, localModelDat, prm )
	if (isempty(vecX_trial))
		vecFModel_trial = [];
		return;
	endif
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		% We don't have a valid model (yet).
		vecFModel_trial = [];
		return;
	endif
	matW = localModelDat.matW;
	sizeV = size(matV,2);
	assert( 0 < sizeV );
	assert( isrealarray(matV,[sizeX,sizeV]) );
	assert( isrealarray(matW,[sizeF,sizeV]) );
	assert( reldiff(matV'*matV,eye(size(matV,2)),eps) <= sqrt(eps) );
	%
	assert( isrealarray(vecX_trial,[sizeX,1]) );
	vecD = vecX_trial - vecX0;
	vecFModel_trial = vecF0 + matW*(matV'*vecD);
	%
	matPhi = localModelDat.matPhi;
	matGamma = localModelDat.matGamma;
	sizePhi = size(matPhi,2);
	if ( 0 < sizePhi )
		assert( isrealarray(matPhi,[sizeX,sizePhi]) );
		assert( isrealarray(matGamma,[sizeF,sizePhi]) );
		assert( reldiff(matPhi'*matPhi,eye(sizePhi),eps) <= sqrt(eps) );
		assert( reldiff(matPhi'*matPhi,eye(size(matPhi,2)),eps) <= sqrt(eps) );
		%
		for n=1:sizePhi
			vecFModel_trial += 0.5 * matGamma(:,n) * (matPhi(:,n)'*vecD)^2;
		endfor
	endif
return;
endfunction


function p = __findPOfDeltaNorm( deltaNormMax, funchDeltaOfP )
	if (isempty(deltaNormMax))
		p = 1.0;
		return;
	endif
	fzero_fun = @(p_dummy)( norm(funchDeltaOfP(p_dummy)) - deltaNormMax );
	if ( fzero_fun(1.0) <= 0.0 )
		p = 1.0;
		return;
	endif
	p = fzero( fzero_fun, [0.0, 1.0] );
return;
endfunction
