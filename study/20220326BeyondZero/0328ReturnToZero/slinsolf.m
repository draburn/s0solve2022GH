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
		switch (trialAction)
		case "giveUp"
			break;
		case "expandSubspace"
			[ localModelDat, ess_datOut ]  = __expandSubspace( funchF, localModelDat, preconDat, prm );
			fevalCount += ess_datOut.fevalCount;
			if ( isempty(localModelDat) )
				msgif( verbLev >= VERBLEV_MAIN, __FILE__, __LINE__, sprintf( "ALGORITHM BREAKDOWN: __expandSubspace() returned %d.", retCode ) );
				break;
			endif
			continue;
		case "tryStep"
			% Go below.
		otherwise
			error( "Invalid value of trialAction." );
		endswitch
		%
		%
		%
		vecF_trial = funchF( vecX_trial ); fevalCount++;
		trialResult = __determineTrialResult( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, prm );
		error( "HALT!" );
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
	datOut.stepConstraintDat = stepConstraintDat;
	datOut.preconDat = preconDat;
	datOut.localModelDat = localModelDat; % Let someone externally updated the data to new point.
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
	assert( reldiff(matV'*matV,ones(sizeV),eps) <= eps );
	%
	vecG = matW'*vecF0;
	matH = matW'*matW;
	%
	matVPhi = localModelDat.matVPhi;
	matGamma = localModelDat.matGamma;
	sizePhi = size(matVPhi,2);
	if ( 0 < sizePhi )
		assert( isrealarray(matVPhi,[sizeV,sizePhi]) );
		assert( isrealarray(matGamma,[sizeF,sizePhi]) );
		assert( reldiff(matVPhi'*matVPhi,ones(sizePhi),eps) <= sqrt(eps) );
		%
		for n=1:sizePhi
			msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Negative quad term." );
			matH += (vecF0'*matGamma(:,n)) * matVPhi(:,n) * (matVPhi(:,n)');
		endfor
	endif
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
	vecX_trial = vecX0 + vecDeltaNewton;
return;
endfunction



function vecFModel_trial = __calcFModel( vecX_trial, localModelDat, prm )
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
		vecFModel_trial = [];
		return;
	endif
	matW = localModelDat.matW;
	sizeV = size(matV,2);
	assert( 0 < sizeV );
	assert( isrealarray(matV,[sizeX,sizeV]) );
	assert( isrealarray(matW,[sizeF,sizeV]) );
	assert( reldiff(matV'*matV,ones(sizeV),eps) <= eps );
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
		assert( reldiff(matPhi'*matPhi,ones(sizePhi),eps) <= eps );
		%
		for n=1:sizePhi
			vecFModel_trial += 0.5 * matGamma(:,n) * (matPhi(:,n)'*vecD)^2;
		endfor
	endif
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
	dta_c0 = 0.9;
	if ( norm(vecFModel_trial) < dta_c0 * norm(vecF0) )
		trialAction = "tryStep";
		return;
	endif
	trialAction = "expandSubspace";
	return;
return;
endfunction



function [ localModelDat, datOut ]  = __expandSubspace( funchF, localModelDat, preconDat, prm )
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
		return;
	endif
	epsFD = mygetfield( prm, "epsFD", eps^0.3 );
	vecFP = funchF( vecX0 + epsFD*vecV ); fevalCount++;
	vecW = ( vecFP - vecF0 ) / (epsFD*norm(vecV));
	matV = [ matV, vecV ];
	matW = [ matW, vecW ];
	%
	[ matPhi, matGamma, pp_datOut ] = __phiPatch( funchF, vecX0, vecF0, matV, matW, matPhi, matGamma, prm );
	fevalCount += pp_datOut.fevalCount;
	%
	localModelDat.matV = matV;
	localModelDat.matW = matW;
	localModelDat.matPhi = matPhi;
	localModelDat.matVPhi = matV'*matPhi;
	localModelDat.matGamma = matGamma;
	datOut.fevalCount = fevalCount;
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
	for n=1:numPasses
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( 0.0==v )
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
