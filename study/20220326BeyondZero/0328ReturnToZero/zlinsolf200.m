function [ vecX_best, vecF_best, retCode, fevalCount, stepsCount, datOut ] = zlinsolf200( funchF, vecX_initial, vecF_initial=[], prmIn=[] )
	% INIT
	mydefs;
	vecX_best = [];
	vecF_best = [];
	retCode = RETCODE__NOT_SET;
	fevalCount = 0;
	stepsCount = 0;
	datOut = [];
	%
	[ retCode, fevalIncr, vecF_initial, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn );
	fevalCount += fevalIncr; clear fevalIncr;
	if ( 0~= retCode )
		msgretcodeif( true, __FILE__, __LINE__, retCode );
		return;
	endif
	vecF_best = vecF_initial;
	%
	%
	% MAIN LOOP
	vecX = vecX_initial;
	vecF = vecF_initial;
	vecX_next = [];
	vecF_next = [];
	fModelDat = [];
	iterCount = 0;
	while (1)
		% Tier 0a actions: simple stoping criteria.
		if ( norm(vecF_best) <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: norm(vecF_best) <= prm.fTol." );
			retCode = RETCODE__SUCCESS;
			break;
		elseif ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterMax." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		elseif ( fevalCount >= prm.fevalMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: fevalCount >= prm.fevalMax." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		elseif ( stepsCount >= prm.stepsMax )
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
		% Tier 0b actions: mandatory actions.
		if ( isempty(fModelDat) )
			[ retCode, fevalIncr, fModelDat ] = __initModel( funchF, vecX, vecF, prm );
			fevalCount += fevalIncr; clear fevalIncr;
			if ( 0~=retCode )
				msg( __FILE__, __LINE__, "Passing an error." );
				break;
			endif
			continue;
		endif
		%
		[ retCode, fevalIncr, fModelDat ] = __analyzeModel( fModelDat, prm );
		fevalCount += fevalIncr; clear fevalIncr;
		if ( 0~= retCode )
			msgretcodeif( true, __FILE__, __LINE__, retCode );
			break
		endif
		%
		%
		% Tier 1 actions: rather clear-cut.
		
		% If we have a "next" candidate and it looks sufficiently good, move to it.
		
		% If high variation omega of zero variation (bounded) step is below cnvg tol, try it.
		
		% If zero variation omega of (zero variation) unbounded step is large, expand subspace.
		% We'll be particularly inclined to do this if we have yet to take a step.
		% But we cannot do this is we've already explored the full X space.
		
		% If zero variation omega of zero variation (bounded) step offers only a too small decrease,
		%  declare that the trust region size has gotten too small.
		
		% If hivar omega of hivar (b) step offer suffic reduction, try it.
		
	endwhile
return;
endfunction


function [ retCode, fevalIncr, vecF_initial, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn )
	% If an error is thrown during init, let it be thrown for useful reporting to the command window.
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
	%
	prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__ZERO; % Production.
	prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration.
	prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Performance testing.
	prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	prm.iterMax = ceil( 100 + 10*sqrt(sizeX+sizeF) + sizeX );
	prm.fevalMax = prm.iterMax;
	prm.stepsMax = 100;
	prm.fTol = sizeF*100.0*eps;
	prm.epsFD = 1.0e-3;
	prm.orthoTol = 1.0e-10;
	prm.curveType = "lev";
	prm.curveScaling = "b";
	prm.precon_funchPrecon = [];
	prm.precon_matL = [];
	prm.precon_matU = [];
	prm.precon_matJA = [];
	prm.fModelDat_initial = [];
	prm.matC = [];
	prm.cholTol = sqrt(eps);
	prm.epsRegu = sqrt(eps);
	prm.candidateStepSizeRelTol = 0.2;
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
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( isrealscalar(prm.verbLev) );
		assert( isrealscalar(prm.valdLev) );
		%
		assert( isrealarray(vecX_initial,[sizeX,1]) );
		assert( isrealarray(vecF_initial,[sizeF,1]) );
		%
		assert( isrealscalar(prm.iterMax) );
		assert( isrealscalar(prm.fevalMax) );
		assert( isrealscalar(prm.stepsMax) );
		assert( abs(prm.iterMax-round(prm.iterMax)) < sqrt(eps) );
		assert( abs(prm.fevalMax-round(prm.fevalMax)) < sqrt(eps) );
		assert( abs(prm.stepsMax-round(prm.stepsMax)) < sqrt(eps) );
		assert( 0 <= prm.iterMax );
		assert( 0 <= prm.fevalMax );
		assert( 0 <= prm.stepsMax );
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
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __initModel( funchF, vecX, vecF, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	fModelDat = [];
	%
	vecRhoF = vecF;
	vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
	vecV = __calcOrthonorm( vecU, [], prm );
	if ( norm(vecV) <= sqrt(eps) )
		msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "Krylov failed to generate an linearly independent vector." );
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
	fevalIncr ++;
	%
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matALo = [ 0.0 ];
	fModelDat.matAHi = [ 0.0 ];
	fModelDat.matB = [ 0.0 ];
	%
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matVLocal = [ vecV ];
	fModelDat.matWLocal = [ vecW ]; % Only for dev/debug?
	%
	fModelDat.matALo_frozen = [ 0.0 ]; % Only for dev/debug?
	fModelDat.matAHi_frozen = [ 0.0 ]; % Only for dev/debug?
	fModelDat.matB_frozen = [ 0.0 ]; % Only for dev/debug?
	%
	retCode = RETCODE__SUCCESS;
	return;
endfunction


function vecU = __applyPrecon( vecRhoF, prm, vecX, vecF )
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
	u0 = norm(vecU);
	if (0.0==u0)
		vecV = zeros(size(vecU));
		return;
	elseif (isempty(matV))
		vecV = vecU/u0;
		return;
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
	mydefs;
	v = norm(vecV);
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < v );
	endif
	vecFP = funchF( vecX + prm.epsFD*vecV );
	vecW = ( vecFP - vecF ) / (prm.epsFD * v);
	return;
endfunction


function [ retCode, fevalIncr, fModelDat ] = __analyzeModel( fModelDat, prm )
	mydefs;
	retCode = RETCODE__NOT_SET;
	fevalIncr = 0;
	%
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matALo = fModelDat.matALo; % Low estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matAHi = fModelDat.matAHi; % High estimate for Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; (bound) candidate steps must satify ||B*y|| <= 1.
	vecX = fModelDat.vecX; % Current guess.
	vecF = fModelDat.vecF; % Function at current guess.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Only for dev/debug?
	matALo_frozen = fModelDat.matALo_frozen; % Only for dev/debug?
	matAHi_frozen = fModelDat.matAHi_frozen; % Only for dev/debug?
	matB_frozen = fModelDat.matB_frozen; % Only for dev/debug?
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,1);
	sizeVLocal = size(matVLocal,2);
	if ( prm.valdLev >= VALDLEV__HIGH )
		assert( sizeB == sizeV ); % Not strictly required, but it's what we're doing.
		assert( isrealarray(matV,[sizeX,sizeV]) );
		assert( isrealarray(matW,[sizeX,sizeV]) );
		assert( isrealarray(matALo,[sizeV,sizeV]) );
		assert( isrealarray(matAHi,[sizeV,sizeV]) );
		assert( isrealarray(matB,[sizeB,sizeV]) );
		assert( isrealarray(matVLocal,[sizeX,sizeVLocal]) );
		assert( isrealarray(matWLocal,[sizeX,sizeVLocal]) );
		assert( isrealarray(matALo_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matAHi_frozen,[sizeV,sizeV]) );
		assert( isrealarray(matB_frozen,[sizeB,sizeV]) );
	endif
	%
	matH = matW'*matW;
	vecG = matW'*vecF;
	matBTB = matB'*matB;
	matIV = eye(sizeV,sizeV);
	% "eta" is estimate for cost function;
	% "omega" is observed values.
	funchEta_zeroV = @(y)( sumsq(vecF)/2.0 + vecG'*y + (y'*matH*y)/2.0 );
	funchEta_loVar = @(y)( funchEta_zeroV(y) + (y'*matALo*y)/2.0 );
	funchEta_hiVar = @(y)( funchEta_zeroV(y) + (y'*matAHi*y)/2.0 );
	%
	if ( ~isempty(prm.matC) )
		matC = prm.matC;
	else
		switch ( tolower(prm.curveScaling) )
		case { "1", "eye", "i" }
			matC = matIV;
		case { "b", "btb", "boundary", "optimal" }
			matC = matBTB;
		case { "ddbtb" }
			matC = diag(diag(matBTB));
		case { "wtw", "newton" }
			matC = matH;
		case { "m", "marq", "marquardt", "ddwtw" }
			matC = diag(diag(matH));
		otherwise
			error( "Invalid value of curveScaling." );
		endswitch
		cScale = norm(diag(matC));
		if ( 0.0 == cScale )
			msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "Curve scaling matrix was zero; setting to I." );
			matC = matIV;
		else
			[ matRC, cholFlag ] = chol(matC);
			if ( 0 ~= cholFlag || min(diag(matRC)) < prm.cholTol * max(abs(diag(matRC))) )
				msgif( prm.verbLev >= VERBLEV__DETAILS, __FILE__, __LINE__, "Curve scaling matrix was non positive-definite; applying regularization." );
				matC += cScale * prm.epsRegu * matIV;
			endif
			cScale = norm(diag(matC));
			matRC = chol(matC);
		endif
	endif
	%
	vecY_ideal = __findCandidateStep( matH, vecG, matC, [], [], prm );
	vecY_zeroV = __findCandidateStep( matH, vecG, matC, matB, 1.0, prm );
	vecY_loVar = __findCandidateStep( matH + matALo, vecG, matC, matB, 1.0, prm );
	verY_hiVar = __findCandidateStep( matH + matAHi, vecG, matC, matB, 1.0, prm );
	if ( prm.valdLev >= VALDLEV__MEDIUM )
		assert( norm(matB*vecY_zeroV) <= 1.0 + prm.candidateStepSizeRelTol );
		assert( norm(matB*vecY_loVar) <= 1.0 + prm.candidateStepSizeRelTol );
		assert( norm(matB*vecY_hiVar) <= 1.0 + prm.candidateStepSizeRelTol );
		assert( funchEta_zeroV(vecY_ideal) <= (1.0+sqrt(eps))*abs(funchEta_zeroV(vecY_zeroV)) + eps*sumsq(vecF) );
		assert( funchEta_zeroV(vecY_zeroV) <= (1.0+sqrt(eps))*abs(funchEta_zeroV(vecY_loVar)) + eps*sumsq(vecF) );
		assert( funchEta_zeroV(vecY_loVar) <= (1.0+sqrt(eps))*abs(funchEta_zeroV(vecY_hiVar)) + eps*sumsq(vecF) );
		assert( funchEta_zeroV(vecY_hiVar) <= (1.0+sqrt(eps))*sumsq(vecF)/2.0 );
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
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
			for n=1:100
				vecY_temp = yNorm + 0.1*(2.0*rand(sizeX,1)-1.0);
				b_temp = norm(matB*vecY_temp);
				eta_temp = funchEta(vecY_temp);
				assert( eta_temp >= eta || b_temp >= b );
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
	endif
	%
	fModelDat.vecY_ideal = vecY_ideal;
	fModelDat.vecY_zeroV = vecY_zeroV;
	fModelDat.vecY_loVar = vecY_loVar;
	fModelDat.vecY_hiVar = vecY_hiVar;
	%
	fModelDat.etaZeroV_ideal = funchEta_zeroV( vecY_ideal );
	fModelDat.etaZeroV_zeroV = funchEta_zeroV( vecY_zeroV );
	fModelDat.etaZeroV_loVar = funchEta_zeroV( vecY_loVar );
	fModelDat.etaZeroV_hiVar = funchEta_zeroV( vecY_hiVar );
	%
	fModelDat.etaLoVar_ideal = funchEta_loVar( vecY_ideal );
	fModelDat.etaLoVar_zeroV = funchEta_loVar( vecY_zeroV );
	fModelDat.etaLoVar_loVar = funchEta_loVar( vecY_loVar );
	fModelDat.etaLoVar_hiVar = funchEta_loVar( vecY_hiVar );
	%
	fModelDat.etaHiVar_ideal = funchEta_hiVar( vecY_ideal );
	fModelDat.etaHiVar_zeroV = funchEta_hiVar( vecY_zeroV );
	fModelDat.etaHiVar_loVar = funchEta_hiVar( vecY_loVar );
	fModelDat.etaHiVar_hiVar = funchEta_hiVar( vecY_hiVar );
	%
	return;
endfunction
%prm.curveType = "lev";
%prm.curveScaling = "b";



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vecY = __calcStep( matH, vecMG, prm )
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
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
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
	if ( isempty(matB) || 0.0==max(max(abs(matB))) )
		vecY = __calcStep( matH, vecMG, prm );
		return;
	endif
	%
	if ( mygetfield(prm,"useDogLeg",false) )
		if (~prm.useBBall)
			error( "useDogLeg requires useBBall!" );
		endif
		%
		% DRaburn 2022-05-12-2233
		%  I now believe only S_curve = B^T * B is meaningful;
		%   this should produce the same point as generating points along a curve and pulling to the surface,
		%   if not better.
		%  However, in order to allow a fair comparison to the existing Levenberg curve,
		%   I'll allow B and S_curve to be independent.
		%
		sizeV = size(matH,1);
		rc = rcond( matH );
		if ( rc > sqrt(eps) )
			matHRegu = matH;
		else
			epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
			hScale = max(max(abs(matH)));
			sScale = max(max(abs(matSCurve)));
			assert( 0.0 < hScale );
			assert( 0.0 < sScale );
			epsRegu = epsRelRegu*hScale/sScale;
			matHRegu = matH + epsRegu*eye(sizeV,sizeV);
		endif
		matR = chol(matHRegu);
		vecDeltaNewton = matR \ ( matR' \ vecMG );
		%
		% If Newton step is in bounds, take it.
		if ( sumsq(matB*vecDeltaNewton) <= 1.0 )
			vecY = vecDeltaNewton;
			return;
		endif
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "Restricting step to bound (%g).",  sumsq(matB*vecDeltaNewton) ) );
		%
		assumeDiagonal = false;
		if (assumeDiagonal)
		vecS = diag(matSCurve);
		matS = diag(vecS);
		assert( reldiff(matSCurve,matS) < sqrt(eps) );
		assert( min(vecS) > eps*max(abs(vecS)) );
		vecYCauchyDir = vecMG./vecS;
		else
		matRCurve = chol(matSCurve);
		vecYCauchyDir = matRCurve \ (matRCurve'\vecMG);
		endif
		%
		ythy = vecYCauchyDir'*matHRegu*vecYCauchyDir;
		assert( ythy >= 0.0 );
		ythmg = vecYCauchyDir'*vecMG;
		assert( ythmg > 0.0 );
		pCauchy = ythmg / ythy;
		vecDeltaCauchy = pCauchy*vecYCauchyDir;
		%
		% If Cauchy step goes OOB, take its intersection with boundary.
		vecBC = matB*vecDeltaCauchy;
		bcsq = sumsq( vecBC );
		if ( bcsq >= 1.0 );
			vecY = vecDeltaCauchy / sqrt(bcsq);
			assert( reldiff( sumsq(matB*vecY), 1.0 ) < sqrt(eps) );
			return;
		endif
		%
		% Find where the second leg intersects the boundary.
		vecDelta2 = vecDeltaNewton - vecDeltaCauchy;
		vecB2 = matB*vecDelta2;
		% We can't use calcLinishRootOfQuad() here!
		% We want positive root of quad!
		%   Using y = vecDeltaCauchy + t * vecDelta2,
		%   ||B*y|| = (t^2)*(b2'*b2) + t*(2.0*b2'*bc) + (bc'*bc) - 1.0,
		%   where b2 = B*vecDelta2 and bc = B*vecDeltaCauchy.
		a = sumsq(vecB2);
		b = 2.0*(vecBC'*vecB2);
		c = bcsq-1.0;
		discrim = (b^2) - (4.0*a*c);
		assert( discrim >= 0.0 );
		assert( 0.0 < a );
		t = (-b+sqrt(discrim))/(2.0*a); % Because a must be positive.
		vecY = vecDeltaCauchy + (t*vecDelta2);
		assert( reldiff( sumsq(matB*vecY), 1.0 ) < sqrt(eps) );
		assert( t >= 0.0 );
		assert( t <= 1.0 );
		return;
	endif
	%
	if (1)
		% DRaburn 2022-05-22...
		%  New to zlinsolf150.
		if (~prm.useBBall)
			error( "useBBall is required." );
		endif
		%%%assert( reldiff(matB'*matB,matSCurve) < sqrt(eps) );
		% In princple, we should be able to short-circuit form here.
		% DRaburn 2022-05-22:
		%  Looks like findLevPt_0522 makes things worse,
		%  perhaps because it semi-requires the constrant matrix to be positive-definite?
		vecY = findLevPt_0522( -vecMG, matH, 1.0, matB );
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
	if ( rc > 10.0*sqrt(eps) )
		s1 = 1.0;
	else
		epsRelRegu = mygetfield( prm, "epsRelRegu", 10.0*sqrt(eps) );
		hScale = max(max(abs(matH)));
		sScale = max(max(abs(matSCurve)));
		assert( 0.0 < hScale );
		assert( 0.0 < sScale );
		s1 = 1.0 - (epsRelRegu*hScale/(epsRelRegu*hScale+sScale));
	endif
	assert( s1 >= 0.0 );
	%
	funchYOfS = @(s)( ( s*matH + (1.0-s)*matSCurve ) \ (s*vecMG) );
	assert( rcond(matSCurve) > eps^1.5 );
	assert( rcond( s1*matH + (1.0-s1)*matSCurve ) > eps^1.5 );
	%
	if (prm.useBBall)
	funchBOfY = @(y)( sumsq(matB*y) );
	else
	% NOTE: IF NOT USING DOG-LEG,
	% WE SHOULD GENERATE A BUNCH OF POINTS AND PULL THEM IN TO THE SURFACES!
	msgif( prm.msgCopious, __FILE__, __LINE__, "NextVer: Generate many points and pull to satisfy B as in slinsolf200." );
	funchBOfY = @(y)( max(abs(y'*matB)) );
	endif
	%
	vecY1 = funchYOfS(s1);
	b1 = funchBOfY(vecY1);
	if ( b1 <= 1.0 )
		vecY = vecY1;
		return;
	endif
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, "Restricting step to bound!" );
	funchBM1OfS = @(s)( funchBOfY(funchYOfS(s)) - 1.0 );
	%
	%
	if (1)
		% DRaburn 2022-05-22...
		%  New to zlinsolf150.
		if (~prm.useBBall)
			error( "useBBall is required." );
		endif
		%%%assert( reldiff(matB'*matB,matSCurve) < sqrt(eps) );
		vecY = findLevPt_0522( -vecMG, matH, 1.0, matB );
	else
		s = fzero( funchBM1OfS, [0.0, s1] );
		vecY = funchYOfS(s);
	endif
	%
	if (prm.useBBall)
	%msg( __FILE__, __LINE__, sprintf( "BBall: %f, %f.", sumsq(matB*vecY1), sumsq(matB*vecY) ) );
	assert( reldiff(norm(matB*vecY),1.0) < sqrt(eps) );
	assert( reldiff(sumsq(matB*vecY),1.0) < sqrt(eps) );
	else
	assert( reldiff(max(abs(vecY'*matB)),1.0) < sqrt(eps) );
	endif
return;
endfunction


function [ fModelDat, datOut ] = __expandModel( vecRhoF, funchF, fModelDat, prm )
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
	fevalCount = 0;
	% Unpack.
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matA0 = fModelDat.matA0;
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,2);
	%
	%
	assert( 0.0 < norm(vecRhoF) );
	vecU = __precon( vecRhoF, prm, vecX, vecF );
	assert( 0.0 < norm(vecU) );
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
	if (isempty(matVLocal))
		fModelDat.matVLocal = [ vecV ];
		fModelDat.matWLocal = [ vecW ];
	else
		fModelDat.matVLocal = [ matVLocal, vecV ];
		fModelDat.matWLocal = [ matWLocal, vecW ];
	endif
	%
	fModelDat.matV = [ matV, vecV ];
	fModelDat.matW = [ matW, vecW ];
	fModelDat.matA = zeros(sizeV+1,sizeV+1);
	fModelDat.matA(1:sizeV,1:sizeV) = (matA'+matA)/2.0;
	if (prm.useBBall)
	fModelDat.matB = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = (matB'+matB)/2.0;
	else
	fModelDat.matB = [ matB; zeros(1,sizeB) ];
	endif
	fModelDat.matA0 = zeros(sizeV+1,sizeV+1);
	fModelDat.matA0(1:sizeV,1:sizeV) = (matA0'+matA0)/2.0;
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction


function fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm )
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matA0 = fModelDat.matA0;
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
	%%%msgif( prm.msgCopious, __FILE__, __LINE__, "TODO: consider partial quadratic update (OSQU)." );
	%%%msgif( prm.msgCopious, __FILE__, __LINE__, "  Here, only linear (Broyden) update is applied." );
	vecFModel_trial = vecF + matW*vecY;
	vecRho = vecF_trial - vecFModel_trial;
	useQuadUpdate = mygetfield( prm, "useQuadUpdate", true );
	if (~useQuadUpdate)
	matW_plus = matW + vecRho * (vecY')/(yNorm^2);
	else
	rhoSumsq = sumsq(vecRho);
	ytay = vecY'*matA*vecY;
	if ( rhoSumsq <= ytay )
		s = 1.0;
	else
		s = ytay/rhoSumsq;
	endif
	matW_plus = matW + (2.0-s) * vecRho * (vecY')/(yNorm^2);
	endif
	%
	%
	%
	%%%
	%%% STUDY ME!
	vecYHat = vecY/yNorm;
	if (0)
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
	coeffD = mygetfield( prm, "coeffD", 10.0 );
	s*=coeffD;
	matA0 = matE*( matA + s*matD )*matE;
	endif
	%%%matA0 = matA + s*matD;
	%
	matWTW = matW'*matW;
	matA0 = 100.0*( diag(diag(matWTW)) + eye(sizeV,sizeV)*max(max(matWTW))*0.1 );
	%%%matA0 += matA + 0.00001*diag(diag(matWTW));
	matA0 = (matA0'+matA0)/2.0;
	%%%
	%%%
	%
	if (0)
		msg( __FILE__, __LINE__, "Data dump..." );
		vecD = diag(matWTW);
		matD = diag(vecD);
		dAvg = sum(vecD)/sizeV;
		dSqAvg = sum(vecD.^2)/sizeV;
		dVar = sqrt( dSqAvg - dAvg^2 );
		dMin = min(vecD);
		dMax = max(vecD);
		[ sumsq(vecRho), dAvg, dVar, dMin, dMax ]
		[ sumsq(vecRho)/sumsq(vecY), sumsq(vecRho)*dAvg/sumsq(matW*vecY), sumsq(vecRho)*dAvg/(vecY'*matD*vecY) ]
		error( "HALT!" );
	endif
	%
	%
	%
	fModelDat.matVLocal = [];
	fModelDat.matWLocal = [];
	fModelDat.vecX = vecX_trial;
	fModelDat.vecF = vecF_trial;
	fModelDat.matW = matW_plus;
	fModelDat.matA = matA0;
	fModelDat.matA0 = matA0;
return;
endfunction


function [ fModelDat, datOut ] = __refresh( vecY, funchF, fModelDat, prm )
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
	fevalCount = 0;
	% Unpack.
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal;
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA0 = fModelDat.matA0; % Hessian variation matrix, < (delta W)' * (delta W) >.
	%matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeVLocal = size(matVLocal,2);
	%sizeB = size(matB,2);
	%
	assert( sizeVLocal < sizeV );
	%
	yNorm = norm(vecY);
	assert( 0.0 < yNorm );
	vecYHat = vecY/yNorm;
	vecU = matV*vecYHat;
	vecV = __calcOrthonorm( vecU, matVLocal, prm );
	%%
	doSaveMeHack20220511 = true;
	if (doSaveMeHack20220511)
		clear vecYHat;
		if ( 0.0 == norm(vecV) )
			vecU = matV*fModelDat.vecYPB;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		if ( 0.0 == norm(vecV) )
			vecU = matV*fModelDat.vecYIB;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		if ( 0.0 == norm(vecV) )
			vecU = matV*fModelDat.vecYIU;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		% Expand model?
		if ( 0.0 == norm(vecV) )
			vecU = matW'*fModelDat.vecF;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		% Try a random vector???
		if ( 0.0 == norm(vecV) )
			vecU = (1:sizeV)';
			vecV = __calcOrthonorm( vecU, matVLocal, prm )
		endif
		if ( 0.0 == norm(vecV) )
			echo__matV = matV
			echo__matVLocal = matVLocal
			error( "Failed to generate local subspace basis vector." );
		endif
	endif
	%%
	if ( 0.0 == norm(vecV) )
		error( "Failed to generate local subspace basis vector." );
	endif
	if (prm.debugMode)
		prmTemp.orthoTol = 10.0*sqrt(eps);
		foo = __calcOrthonorm( vecV, matV, prmTemp );
		if ( norm(foo) != 0.0 )
			msg( __FILE__, __LINE__, "Data dump..." );
			[ norm(vecU), norm(vecV), norm(foo) ]
			[ norm(vecU), norm(vecV), norm(foo) ] - 1.0
			[ vecV'*vecU, foo'*vecU, foo'*vecV ]
			[ norm( matV'*vecU ), norm( matV'*vecV ), norm( matV'*foo ) ]
			[ norm( vecV - matV*(matV'*vecV) ), norm( vecU - matV*(matV'*vecU) ), norm( foo - matV*(matV'*foo) ) ]
			[ norm( vecV - matV*(matV'*vecV) ), norm( vecU - matV*(matV'*vecU) ), norm( foo - matV*(matV'*foo) ) ] - 1.0
			error( "Subspace basis vector was thrown out of own space?!?!" );
		endif
	endif
	%
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecW) )
		error( "Jacobian along local subspace basis vector is zero." );
	endif
	%
	matVLocal = [ matVLocal, vecV ];
	matWLocal = [ matWLocal, vecW ];
	sizeVLocal++;
	if ( prm.debugMode )
		assert( reldiff( matVLocal, matV*(matV'*matVLocal), eps ) < sqrt(eps) );
		assert( reldiff( eye(sizeVLocal,sizeVLocal), matVLocal'*matVLocal, eps ) < sqrt(eps) );
	endif
	matVLocal = matV*(matV'*matVLocal);
	for n=1:sizeVLocal
		matVLocal(:,n) /= norm(matVLocal(:,n));
	endfor
	%
	matE = eye(sizeV,sizeV) - matV'*matVLocal*(matVLocal')*matV;
	%
	fModelDat.matVLocal = matVLocal;
	fModelDat.matWLocal = matWLocal;
	%
	%%% The yHats are NOT perpendicular!
	%%%fModelDat.matW = matW + (vecW - matW*vecYHat)*(vecYHat');
	%%% This z thing also works..
	%%%vecZ = matV'*vecV;
	%%%fModelDat.matW = matW + (vecW - matW*vecZ)*(vecZ');
	% But, to be safe, let's do this.
	fModelDat.matW = matW + ( matWLocal - matW*(matV')*matVLocal)*(matVLocal'*matV); % Yeah?
	%
	matA = matE' * matA0 * matE;
	%matA += sqrt(eps)*max(max(abs(matA)))*eye(sizeV,sizeV);
	matA = (matA'+matA)/2.0;
	fModelDat.matA = matA;
	if (0)
		vecLambda0 = eig(matA0);
		vecLambda = eig(matA);
		if ( min(vecLambda0) < 0.0 || min(vecLambda) < 0.0 )
			vecLambda0
			vecLambda
			assert( min(vecLambda0) >= 0.0 );
			assert( min(vecLambda) >= 0.0 );
		endif
	endif
	%
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function fModelDat = __addB( vecY, fModelDat, prm )
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
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
	if (prm.useBBall)
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "Increasing B (%g).", sumsq(matB*vecY) ) );
		assert( sumsq(matB*vecY) < 1.0+sqrt(eps) );
		yNorm = norm(vecY);
		assert( 0.0 < yNorm );
		vecYHat = vecY/yNorm;
		sizeV = size(matB,2);
		matEY = eye(sizeV,sizeV) - vecYHat*(vecYHat');
		matB = matEY*matB*matEY + vecYHat*(vecYHat'/yNorm);
		assert( reldiff( sumsq(matB*vecY), 1.0 ) < sqrt(eps) );
		fModelDat.matB = (matB'+matB)/2.0;
		return;
	endif
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
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
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
	if (prm.useBBall)
		if ( sumsq(matB*vecY) < 1.0 )
			% vecY is already in the trust region; nothing to do.
			return;
		endif
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "Decreasing B (%g).", sumsq(matB*vecY) ) );
		assert( sumsq(matB*vecY) > 1.0 - sqrt(eps) );
		yNorm = norm(vecY);
		assert( 0.0 < yNorm );
		vecYHat = vecY/yNorm;
		sizeV = size(matB,2);
		matEY = eye(sizeV,sizeV) - vecYHat*(vecYHat');
		matB = matEY*matB*matEY + vecYHat*(vecYHat'/yNorm);
		if ( reldiff(sumsq(matB*vecY),1.0) >= sqrt(eps) )
			msg( __FILE__, __LINE__, "Data dump..." );
			matB
			vecY
			sumsq(matB*vecY)
			reldiff(sumsq(matB*vecY),1.0)
		endif
		assert( reldiff( sumsq(matB*vecY), 1.0 ) < sqrt(eps) );
		fModelDat.matB = (matB'+matB)/2.0;
		return;
	endif
	%
	if ( isempty(matB) )
		return;
	endif
	msk = logical( abs(vecY'*matB) >= 1.0 );
	fModelDat.matB = matB( :, msk );
return;
endfunction

function __dumpModel( fModelDat, prm )
	error( "DRABURN: 2022-05-26: HIT OLDVER CODE!" );
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
