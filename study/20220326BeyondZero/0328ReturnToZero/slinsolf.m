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
	%
	localModelDat = mygetfield( datIn, "localModelDat", [] );
	if (isempty(localModelDat))
		localModelDat.vecX0 = vecX0;
		localModelDat.vecF0 = vecF0;
	else
		assert( reldiff(vecX0,localModelDat.vecX0,realmin) <= eps )
		assert( reldiff(vecF0,localModelDat.vecF0,realmin) <= eps )
	endif
	%
	%
	% Prep loop.
	vecX_best = vecX0;
	vecF_best = vecF0;
	vecFModel_best = vecF0;
	stepConstraitDat_best = []; % Essentially, SCD *implied by* best.
	bestIsNot0 = false;
	%
	vecX = vecX0;
	vecF = vecF0;
	iterCount = 0;
	while (1)
		iterMax = mygetfield( prm, "iterMax", 10+sizeX );
		if ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV_MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
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
		error( "HALT!" );
		trialAction = __determineTrialAction( vecX_trial, vecFModel_trial, localModelDat, prm );
		switch (trialAction)
		case "giveUp"
			break;
		case "expandSubspace"
			[ localModelDat, fevalIncrement ]  = __expandSubspace( localModelDat, preconDat, prm );
			fevalCount += fevalIncrement;
			if ( 0~=retCode )
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
		%
		switch (trialResult)
		case "excellent"
			% Update (loosen) SCD, then accept step.
			stepConstraintDat = __updateSCD_excellent( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			vecFModel_best = vecFModel_trial;
			stepConstraintDat_best = stepConstraintDat;
			bestIsNot0 = true;
			break;
		case "good"
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			vecFModel_best = vecFModel_trial;
			stepConstraintDat_best = stepConstraintDat;
			bestIsNot0 = true;
			break;
		case "okay"
			if ( bestIsNot0 )
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
			vecFModel_best = vecFModel_trial;
			stepConstraintDat_best = stepConstraintDat;
			bestIsNot0 = true;
			stepConstraintDat = __updateSCD_okay( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			continue;
		case "bad"
			if ( bestIsNot0 )
				msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Trial result was 'bad' even though a previous result was at least 'okay'." );
			endif
			stepConstraintDat = __updateSCD_bad( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			continue;
		case "horrid"
			if ( bestIsNot0 )
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
	datOut.localModelDat = localModelDat; % Let this be updated to new point externally.
	if ( ~bestIsNot0 )
		msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Failed to find a good new point." );
	endif
	%
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
	matPhi = localModelDat.matPhi;
	matGamma = localModelDat.matGamma;
	sizePhi = size(matPhi,2);
	if ( 0 < sizePhi )
		assert( isrealarray(matPhi,[sizeX,sizePhi]) );
		assert( isrealarray(matGamma,[sizeF,sizePhi]) );
		assert( reldiff(matPhi'*matPhi,ones(sizePhi),eps) <= eps );
		%
		for n=1:sizePhi
			msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Negative quad term." );
			matH += (vecF0'*matGamma(:,n)) * matPhi(:,n) * (matPhi(:,n)');
		endfor
	endif
	%
	%
	if ( mygetfield( prm, "useMarquardtScaling", false ) )
		hScale = sqrt(sum(sum(matH.^2))/sizeV);
		matS_curve = diag( 1.0 ./ ( sqrt(abs(diag(matH))) + eps*hScale ) ); % Needs testing.
		vecG_curve = matS_curve'*vecG;
		matH_curve = calcHRegu( matS_curve'*matH*matS_curve );
	else
		vecG_curve = vecG;
		matH_curve = calcHRegu( matH );
	endif
	%
	pCauchy = calcLinishRootOfQuad( 0.5*(vecG_curve'*matH_curve*vecG_curve), -sumsq(vecG_curve), omega );
	assert( pCauchy >= 0.0 );
	vecDeltaCauchy = pCauchy*matS_curve*(-vecG_curve); % Needs testing.
	vecDeltaNewton = matS_curve*(matH_curve\(-vecG_curve));
	%
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



function [ localModelDat, fevalIncrement ]  = __expandSubspace( localModelDat, preconDat, prm )
	fevalIncrement = 0;
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
	if ( 1 == size(matW,2) )
		matU = localModelDat.vecF;
	else
		matU = localModelDat.matW(:,end);
	endif
	%
	vecV = __calcGuessV( vecU, preconDat, prm );
	vecV = __calcOrthog( vecV, matV, prm );
	if ( norm(vecV) < eps )
		msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: norm(vecV) < eps." );
		retCode = 100;
		return;
	endif
	vecW = __calcJProd( vecV, prm ); fevalIncrement++;
	matV = [ matV, vecV ];
	matW = [ matW, vecW ];
	%
	[ matPhi, matGamma ] = __phiPatch( funchF, vecX0, vecF0, matV, matW, matPhi, matGamma, prm );
	fevalCount += phiPatch_datOut.fevalCount;
return;
endfunction
