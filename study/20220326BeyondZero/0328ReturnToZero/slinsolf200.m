% slinsolf200: let's focus on trial steps & walls;
%  de-prioritizing quadratic terms.

function [ vecXF, vecFF, datOut ] = slinsolf200( funchF, vecX0, vecF0, prm, datIn )
	% Parse input.
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	debugMode = mygetfield( prm, "debugMode", false );
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
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "trialAction = '%s'.", trialAction ) );
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
				%msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: __expandSubspace() failed." );
				if (isempty(localModelDat.matV))
					error( "__expandSubspate failed for first vector!" );
				endif
				%trialAction = "findGoodStep";
				trialAction = "findGoodStep";
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "  Forcing '%s'.", trialAction ) );
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
			fgsIsBest = false;
			if ( isempty(vecX_best) )
				fgsIsBest = true;
			elseif ( norm(vecF_trial) < norm(vecF_best) )
				fgsIsBest = true;
			endif
			if ( fgsIsBest )
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
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "||F||: %g, %g, %g.", norm(vecF0), norm(vecFModel_trial), norm(vecF_trial) ) );
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "trialResult = '%s'.", trialResult ) );
		%
		trialIsNewBest = false;
		if ( isempty(vecX_best) )
			trialIsNewBest = true;
		elseif ( norm(vecF_trial) < norm(vecF_best) )
			trialIsNewBest = true;
		endif
		if (trialIsNewBest)
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			stepConstraintDat_best = stepConstraintDat;
		endif
		%
		switch (trialResult)
		case "accept"
			break;
		case "reject"
			stepConstraintDat = __updateSCD_addWall( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm );
			continue;
		endswitch
		error( "Inaccessible code." );
	endwhile
	%
	vecXF = vecX_best;
	vecFF = vecF_best;
	vecFModelF = __calcFModel( vecX_best, localModelDat, prm ); % Recalc to ensure in-sync with _best.
	stepConstraintDat = stepConstraintDat_best;
	%
	datOut.stepConstraintDat = stepConstraintDat;
	datOut.preconDat = preconDat;
	datOut.localModelDat = localModelDat; % Let someone externally updated the data to new point.
	datOut.vecFModelF = vecFModelF;
	datOut.fevalCount = fevalCount;
return;
endfunction



function vecX_trial = __calcDelta( localModelDat, stepConstraintDat, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	debugMode = mygetfield( prm, "debugMode", false );
	%
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
	endif
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		% We don't have a valid model (yet).
		vecX_trial = [];
		return;
	endif
	matW = localModelDat.matW;
	sizeV = size(matV,2);
	if (debugMode)
		assert( 0 < sizeV );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( reldiff(matV'*matV,eye(size(matV,2)),eps) <= sqrt(eps) );
	endif
	%
	vecG = localModelDat.vecG;
	matH = localModelDat.matH;
	hScale = sqrt(sum(sum(matH.^2))/sizeV);
	matIV = eye(sizeV,sizeV);
	%
	%
	if ( mygetfield( prm, "useMarquardtScaling", false ) )
		matS_curve = diag( 1.0 ./ ( sqrt(abs(diag(matH))) + eps*hScale ) ); % Needs testing.
	else
		matS_curve = matIV;
	endif
	%
	%
	vecG_curve = matS_curve'*vecG;
	matH_curve = calcHRegu( matS_curve'*matH*matS_curve );
	%
	pCauchy = calcLinishRootOfQuad( 0.5*(vecG_curve'*matH_curve*vecG_curve), -sumsq(vecG_curve), sumsq(vecF0)/2.0 );
	assert( pCauchy >= 0.0 );
	vecDeltaCauchy = pCauchy*matS_curve*(-vecG_curve); % Needs testing.
	vecDeltaNewton = matS_curve*(matH_curve\(-vecG_curve));
	%
	if ( isempty(stepConstraintDat) )
		vecX_trial = vecX0 + matV * vecDeltaNewton;
		return;
	endif
	%
	%
	numWalls = mygetfield( stepConstraintDat, "numWalls", 0 );
	funchYOfP = @(p)(matS_curve*( (p*matH_curve + (1.0-p)*matIV) \ (-p*vecG_curve) ));
	%
	switch ( mygetfield(prm,"calcDeltaMethod","levPull_v100") )
	case "levPull_v000"
		% Generate points along curve, pull inwards to satisfy constraints, find best according to model.
		% (An alteriative would be to use fminunc (or the constrained version, if available)
		%  to find a/the min of the full quadratic F, but, this is good enough for now.
		calcDeltaCurvePrm = mygetfield( prm, "calcDeltaCurvePrm", [] );
		numPts = mygetfield( calcDeltaCurvePrm, "numPts", 101 );
		pMax = mygetfield( calcDeltaCurvePrm, "pMax", 1.0 );
		a0 = mygetfield( calcDeltaCurvePrm, "a0", 2.0 );
		a1 = mygetfield( calcDeltaCurvePrm, "a1", 2.0 );
		assert( isrealscalar(numPts) );
		assert( isrealscalar(pMax) );
		assert( isrealscalar(a0) );
		assert( isrealscalar(a1) );
		assert( 2 <= numPts );
		assert( 0.0 < pMax );
		%assert( pMax <= 1.0 );
		assert( 0.0 < a0 );
		assert( 0.0 < a1 );
		%
		foo = linspace( pMax, 0.0, numPts );
		pOfPts = ( 1.0 - (foo.^a0) ).^a1;
		%
		p_best = 0.0;
		vecY_best = zeros(sizeV,1);
		vecFModel_best = vecF0;
		for indexP = 2 : numPts
			p = pOfPts(indexP);
			vecY_raw = funchYOfP(p);
			vecDelta_raw = matV*vecY_raw;
			%
			vecDelta = vecDelta_raw;
			for indexW = 1 : numWalls
				vecDeltaOfWall = stepConstraintDat.vecDeltaOfWall(:,indexW);
				s0 = vecDeltaOfWall'*vecDeltaOfWall;
				assert( s0 > 0.0 );
				s = vecDelta'*vecDeltaOfWall;
				if ( s > s0 )
					vecDelta *= s0/s;
				endif
			endfor
			%
			vecX = vecX0 + vecDelta;
			vecFModel = __calcFModel( vecX, localModelDat, prm );
			%
			if ( norm(vecFModel) < norm(vecFModel_best) )
				p_best = p;
				vecX_best = vecX;
				vecFModel_best = vecFModel;
			endif
		endfor
		%
		vecX_trial = vecX_best;
		return;
	case "levPull_v100"
		% Generate points along curve, pull inwards to satisfy constraints, find best according to model.
		% (An alteriative would be to use fminunc (or the constrained version, if available)
		%  to find a/the min of the full quadratic F, but, this is good enough for now.
		% But, be faster!
		calcDeltaCurvePrm = mygetfield( prm, "calcDeltaCurvePrm", [] );
		numPts = mygetfield( calcDeltaCurvePrm, "numPts", 101 );
		pMax = mygetfield( calcDeltaCurvePrm, "pMax", 1.0 );
		a0 = mygetfield( calcDeltaCurvePrm, "a0", 2.0 );
		a1 = mygetfield( calcDeltaCurvePrm, "a1", 2.0 );
		assert( isrealscalar(numPts) );
		assert( isrealscalar(pMax) );
		assert( isrealscalar(a0) );
		assert( isrealscalar(a1) );
		assert( 2 <= numPts );
		assert( 0.0 < pMax );
		%assert( pMax <= 1.0 );
		assert( 0.0 < a0 );
		assert( 0.0 < a1 );
		%
		foo = linspace( pMax, 0.0, numPts );
		pOfPts = ( 1.0 - (foo.^a0) ).^a1;
		pOfPts(1) = pOfPts(2)/numPts; % point is unused, but hack to avoid db zero.
		muOfPts = (1.0-pOfPts)./pOfPts;
		%
		vecZOfPts = zeros(sizeV,numPts);
		for indexP = 2 : numPts
			matR = chol( matH_curve + muOfPts(indexP)*matIV ); % Two steps with chol() is slightly faster, AFAIK.
			vecZOfPts =  matR \ ( matR' \ vecG ); % Note minus sign is missing.
		endfor
		vecYOfPts = -matS_curve*vecZOfPts;
		% We could try to do constraints in Y space.
		% But, POITROME.
		%
		vecDeltaOfPts = matV*vecYOfPts; % Needs constraints
		if ( numWalls > 0 )
			vecDeltaOfWalls = stepConstraintDat.vecDeltaOfWall;
			s0OfWalls = sum(vecDeltaOfWalls.*vecDeltaOfWalls,1);
			%
			% Could probably do this witout the numWalls loop, but POITROME.
			for indexW=1:numWalls
				s1OfPts = vecDeltaOfWalls(:,indexW)'*vecDeltaOfPts;
				vecDeltaOfPts += ( s1OfPts > s0OfWalls(indexW) ) .* vecDeltaOfPts .* ( s0OfWalls(indexW)./s1OfPts - 1.0 );
			endfor
		endif
		%
		%
		vecFModelOfPts = __calcFModelOfDelta( vecDeltaOfPts, localModelDat, prm );
		omegaModelOfPts = sum( vecFModelOfPts.*vecFModelOfPts, 1 );
		[ omegaModel_best, indexP_best ] = min(omegaModelOfPts);
		vecX_trial = vecX0 + vecDeltaOfPts(:,indexP_best);
		return;
	otherwise
		error( "Invalid value of calcDeltaMethod." );
	endswitch
	error( "Unreachable code." );
return;
endfunction



function trialAction = __determineTrialAction( vecX_trial, vecFModel_trial, localModelDat, trialActionDat, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	debugMode = mygetfield( prm, "debugMode", false );
	%
	if (isempty(vecX_trial))
		% We don't have a valid trial.
		trialAction = "expandSubspace";
		return;
	elseif (stopsignalpresent())
		msgif( verbLev >= VERBLEV__NOTIFY, __FILE__, __LINE__, "Received stop signal." );
		trialAction = "giveUp";
		return;
	endif
	%
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(vecX_trial,[sizeX,1]) );
		assert( isrealarray(vecFModel_trial,[sizeF,1]) );
	endif
	%
	%vecX_best = mygetfield( trialActionDat, "vecX_best", [] );
	%if (~isempty(vecX_best))
	%	vecF_best = trialActionDat.vecF_best;
	%	assert( isrealarray(vecX_best,[sizeX,1]) );
	%	assert( isrealarray(vecF_best,[sizeF,1]) );
	%endif
	%
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (~isempty(matV))
	if (size(matV,2)>=sizeX)
		% More quad terms may be possible, but no more expanding to do.
		%trialAction = "tryStep";
		trialAction = "findGoodStep";
		return;
	endif
	endif
	%
	% Okay-ish placeholder...
	vecFModel_newton = mygetfield( localModelDat, "vecFModel_newton", 0.0 );
	dta_c0 = mygetfield( prm, "dta_c0", 0.1 );
	dta_c1 = mygetfield( prm, "dta_c1", 0.5 );
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, ...
	  sprintf( "dta_c0 = %10.3e, dta_c1 = %10.3e;  ||F0|| = %10.3e,  ||FM_n|| = %10.3e,  ||FM_t|| = %10.3e.", ...
	  dta_c0, dta_c1, norm(vecF0), norm(vecFModel_newton), norm(vecFModel_trial) ) );
	if ( norm(vecFModel_newton) > dta_c0 * norm(vecF0) )
		trialAction = "expandSubspace";
		return;
	elseif ( norm(vecFModel_trial) > dta_c1 * norm(vecF0) )
		trialAction = "expandSubspace";
		return;
	else
		trialAction = "tryStep";
		return;
	endif
	error( "Unreachable code." );
return;
endfunction



function [ localModelDat, ess_datOut ]  = __expandSubspace( funchF, localModelDat, preconDat, prm )
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	debugMode = mygetfield( prm, "debugMode", false );
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
	if (debugMode)
		assert( reldiff(matV'*matV,eye(size(matV,2)),eps) <= sqrt(eps) );
	endif
	%
	if (1)
	matPhi = [];
	matGamma = [];
	sizePhi = 0;
	matVTPhi = [];
	else
	[ matPhi, matGamma, pp_datOut ] = __phiPatch( funchF, vecX0, vecF0, matV, matW, matPhi, matGamma, prm );
	fevalCount += pp_datOut.fevalCount;
	sizePhi = size(matPhi,2);
	if ( 0 < sizePhi )
		if (debugMode)
		assert( isrealarray(matPhi,[sizeX,sizePhi]) );
		assert( isrealarray(matGamma,[sizeF,sizePhi]) );
		assert( reldiff(matPhi'*matPhi,eye(size(matPhi,2)),eps) <= sqrt(eps) );
		endif
	endif
	matVTPhi = matV'*matPhi;
	endif
	%
	vecG = matW'*vecF0;
	matWTW = matW'*matW;
	matH = matWTW;
	for n=1:sizePhi
		matH += (vecF0'*matGamma(:,n)) * matVTPhi(:,n) * (matVTPhi(:,n)');
	endfor
	%
	vecFModel_newton = vecF0-matW*(matH\vecG);
	localModelDat.matV = matV;
	localModelDat.matW = matW;
	localModelDat.matWTW = matWTW;
	localModelDat.matPhi = matPhi;
	localModelDat.matGamma = matGamma;
	localModelDat.matVTPhi = matVTPhi;
	localModelDat.vecG = vecG;
	localModelDat.matH = matH;
	localModelDat.vecFModel_newton = vecFModel_newton;
	ess_datOut.fevalCount = fevalCount;
	%
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  %10.3e,  %3d,  %3d;  %10.3e,  %10.3e;  %10.3e.", ...
	  norm(vecF0), sizeV, sizePhi, ...
	  rcond(matWTW), rcond(matH), ...
	  norm(vecFModel_newton) ) );
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
	v0 = norm(vecV);
	if (0.0==v0)
		return;
	endif
	orthoTol = mygetfield( prm, "orthoTol", 1.0e-10 );
	for n=1:numPasses
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( v <= orthoTol*v0 )
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
	debugMode = mygetfield( prm, "debugMode", false );
	fevalCount = 0;
	sizeF = size(vecF0,1);
	sizeV = size(matV,2);
	matWTW = matW'*matW;
	[ matPsi, matLambda ] = eig( matWTW );
	if (debugMode)
		assert( isrealarray(matPsi,[sizeV,sizeV]) );
	endif
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
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	debugMode = mygetfield( prm, "debugMode", false );
	fevalCount = 0;
	%msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Performing __findGoodStep()." );
	if (~isempty(stepConstraintDat))
		msgif( verbLev >= VERBLEV__WARN, __FILE__, __LINE__, "WARNING: Incoming step constraints are not yet supported here." );
	endif
	%
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
	endif
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		localModelDat
		localModelDat.matV
		error( "__findGoodStep() was reached with localModelDat.matV being empty. How is this possible?" );
	endif
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
	dtr_c0 = mygetfield( prm, "dtr_c0", 0.9 );
	if ( norm(vecF_trial) < dtr_c0*norm(vecF0) + (1.0-dtr_c0)*norm(vecFModel_trial) )
		trialResult = "accept";
		return;
	endif
	trialResult = "reject";
	return;
return;
endfunction



function stepConstraintDat = __updateSCD_addWall( vecX_trial, vecF_trial, vecFModel_trial, localModelDat, stepConstraintDat, prm )
	wallBTFactor = mygetfield( prm, "wallBTFactor", 0.5 );
	assert( wallBTFactor < 1.0 );
	assert( wallBTFactor > 0.0 );
	vecX0 = localModelDat.vecX0;
	numWalls = mygetfield( stepConstraintDat, "numWalls", 0 );
	numWalls++;
	stepConstraintDat.numWalls = numWalls;
	stepConstraintDat.vecDeltaOfWall(:,numWalls) = wallBTFactor*(vecX_trial-vecX0);
return;
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM ACCESS FUNCTIONS
%


function vecFModel_trial = __calcFModel( vecX_trial, localModelDat, prm )
	debugMode = mygetfield( prm, "debugMode", false );
	if (isempty(vecX_trial))
		vecFModel_trial = [];
		return;
	endif
	vecX0 = localModelDat.vecX0;
	vecF0 = localModelDat.vecF0;
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
	endif
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		% We don't have a valid model (yet).
		vecFModel_trial = [];
		return;
	endif
	matW = localModelDat.matW;
	sizeV = size(matV,2);
	if (debugMode)
		assert( 0 < sizeV );
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeF,sizeV]) );
		assert( reldiff(matV'*matV,eye(size(matV,2)),eps) <= sqrt(eps) );
	endif
	%
	assert( isrealarray(vecX_trial,[sizeX,1]) );
	vecD = vecX_trial - vecX0;
	vecFModel_trial = vecF0 + matW*(matV'*vecD);
	%
	matPhi = localModelDat.matPhi;
	matGamma = localModelDat.matGamma;
	sizePhi = size(matPhi,2);
	if ( 0 < sizePhi )
		if (debugMode)
			assert( isrealarray(matPhi,[sizeX,sizePhi]) );
			assert( isrealarray(matGamma,[sizeF,sizePhi]) );
			assert( reldiff(matPhi'*matPhi,eye(sizePhi),eps) <= sqrt(eps) );
			assert( reldiff(matPhi'*matPhi,eye(size(matPhi,2)),eps) <= sqrt(eps) );
		endif
		%
		for n=1:sizePhi
			vecFModel_trial += 0.5 * matGamma(:,n) * (matPhi(:,n)'*vecD)^2;
		endfor
	endif
return;
endfunction


function vecFModelOfPts = __calcFModelOfDelta( vecDeltaOfPts, localModelDat, prm )
	if (isempty(vecDeltaOfPts))
		vecFModel_trial = [];
		return;
	endif
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		vecFModel_trial = [];
		return;
	endif
	vecFModelOfPts = localModelDat.vecF0 + localModelDat.matW*(matV'*vecDeltaOfPts); % Automatic broadcasting.
	matPhi = mygetfield( localModelDat, "matPhi", [] );
	if (~isempty(matPhi))
		sizePhi = size(matPhi,2);
		matGammaO2 = localModelDat.matGamma / 2.0;
		for n=1:sizePhi
			vecFModelOfPts += matGammaO2 * (matPhi(:,n)'*vecDeltaOfPts).^2;
		endfor
	endif
return;
endfunction


function vecFModel_trial = __calcFModelOfY( vecY, localModelDat, prm )
	if (isempty(vecX_trial))
		vecFModel_trial = [];
		return;
	endif
	%
	matV = mygetfield( localModelDat, "matV", [] );
	if (isempty(matV))
		% We don't have a valid model (yet).
		vecFModel_trial = [];
		return;
	endif
	%
	vecFModel_trial = localModelDat.vecF0 + localModelDat.matW*vecY;
	%
	matVPhi = mygetfield( localModelDat, "matVPhi", [] );
	if (~isempty(matVPhi))
		sizePhi = size(localModelDat.matVPhi,2);
		for n=1:sizePhi
			vecFModel_trial += 0.5 * localModelDat.matGamma(:,n) * (localModelDat.matVPhi(:,n)'*vecY)^2;
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
