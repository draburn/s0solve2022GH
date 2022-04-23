function [ vecXF, vecFF, datOut ] = slinsolf( funchF, vecX0, vecF0, prm, datIn )
	% Parse input.
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	fevalCount = 0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	%
	stepConstraintDat = mygetfield( datIn, "stepConstraintDat", [] );
	preconDat = mygetfield( datIn, "preconDat", [] );
	%
	localModelDat = mygetfield( datIn, "localModelDat", [] );
	localModelDat.vecX0 = vecX0;
	localModelDat.vecF0 = vecF0;
	%
	vecX_best = vecX0;
	vecF_best = vecF0;
	vecFModel_best = vecF0;
	stepConstraitDat_best = []; % Essentially SCD *implied by* best.
	bestIsNot0 = false;
	vecX = vecX0;
	vecF = vecF0;
	iterCount = 0;
	while (1)
		iterMax = mygetfield( prm, "iterMax", 10+sizeX );
		if ( iterCount >= iterMax )
			msg( __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		iterCount++;
		%
		%
		% Calculate trial step.
		error( "HALT!" );
		vecX_trial = __calcDelta( localModelDat, stepConstraintDat, prm );
		vecFModel_trial = __calcFModel( vecX_trial, localModelDat, prm );
		%
		%
		% Decide what to do.
		trialAction = __determineTrialAction( vecX_trial, vecFModel_trial, localModelDat, prm );
		switch (trialAction)
		case "giveUp"
			break;
		case "expandSubspace"
			[ retCode, localModelDat, fevalIncrement ]  = __expandSubspace( localModelDat, preconDat, prm );
			fevalCount += fevalIncrement;
			if ( 0~=retCode )
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
				msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Trial result was 'bad' even though a previous result was 'okay'." );
			endif
			stepConstraintDat = __updateSCD_nad( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			continue;
		case "horrid"
			if ( bestIsNot0 )
				msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Trial result was 'horrid' even though a previous result was 'okay'." );
			endif
			stepConstraintDat = __updateSCD_horrid( vecX_trial, localModelDat, stepConstraintDat, prm );
			continue;
		otherwise
			error( "Invalid value of trialResult." );
		endswitch
		error( "Inaccessible code." );
	endwhile
	%
endfunction



function [ retCode, localModelDat, fevalIncrement ]  = __expandSubspace( localModelDat, preconDat, prm )
	fevalIncrement = 0;
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
endfunction
