% DRaburn 2022.05.07
%  zlinsolf100
%  First concrete attempt at MVP.

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
	%
	while (1)
		iterCount++;
		fModelDat = __analyzeModel( fModelDat, prm );
		msgif( prm.msgCopious, __FILE__, __LINE__, "vvvvv Data dump..." );
		omega = fModelDat.omega
		omegaModelAvgIU = fModelDat.omegaModelAvgIU
		omegaModelAvgPlusVarIU = fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU
		matVLocal = fModelDat.matVLocal
		matA = fModelDat.matA
		msgif( prm.msgCopious, __FILE__, __LINE__, "^^^^^ End data dump." );
		%
		% Simple stoping criteria.
		if ( norm(vecF_best) <= prm.fTol )
			msgif( prm.msgMain, __FILE__, __LINE__, "SUCCESS: sumsq(vecF_best) <= prm.fTol^2." );
			break;
		endif
		%
		if ( iterCount > prm.iterLimit )
			msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterLimit." );
			break;
		endif
		%
		if ( stopsignalpresent() )
			msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		%
		if ( fModelDat.omegaModelAvgIU > prm.omegaTol )
			% Conceptually, we could also consider "refreshing" something already in our space.
			% This is an area for future analysis.
			msgif( prm.msgCopious, __FILE__, __LINE__, "Expanding subspace." );
			[ fModelDat, datOut_expandModel ] = __expandModel( fModelDat.vecFModelIU, funchF, fModelDat, prm );
			fevalCount += datOut_expandModel.fevalCount;
			continue;
		endif
		%
		if ( fModelDat.bIU <= 1.0  && fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU <= prm.omegaTol )
		%%%%%if ( fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU <= prm.omegaTol )
			msgif( prm.msgCopious, __FILE__, __LINE__, "Trying ideal unbound step." );
			vecX_trial = fModelDat.vecXIU;
			%%%%%vecX_trial = fModelDat.vecXIB;
			vecF_trial = funchF( vecX_trial );
			msg( __FILE__, __LINE__, "TODO: Decrease B if model is very good!" );
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
				msgif( prm.msgNotice, __FILE__, __LINE__, "Current trial is worse than earlier candidate; moving to earlier candidate." );
				fModelDat = __moveTo( vecX_cand, vecF_cand, fModelDat, prm );
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
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  actual fall = %g.", fModelDat.omega - omega_trial ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  expected fall = %g.", fModelDat.omega - (fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU) ) );
			if ( omega_trial <= fModelDat.omega - avefaThresh*( fModelDat.omega - (fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU) ) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Accepting step." );
				fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm );
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
			bUpFactor = mygetfield( prm, "bUpFactor", 0.5 );
			fModelDat = __increaseB( bUpFactor*fModelDat.vecYIU, fModelDat, prm );
			%
			clear vecX_trial;
			clear vecF_trial;
			clear omega_trial;
			continue;
		endif
		%
		minRelFallThresh = mygetfield( prm, "minRelFallThresh", 1.0E-4 );
		if ( fModelDat.omegaModelAvgIB > fModelDat.omega*(1.0-minRelFallThresh) )
			msgif( prm.msgMain, "We seem to have no way to reduce omega much." );
			msgif( prm.msgMain, "  This is expected to happen near a bad local minimum." );
			msgif( prm.msgMain, "  Todo: add handling for this case." );
			break;
		endif
		%
		%practicalRelFallThresh = mygetfield( prm, "practicalRelFallThresh", 0.5 );
		%if ( fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB <= fModelDat.omega - practicalRelFallThresh*(fModelDat.omega-fModelDat.omegaModelAvgIB) )
		%	error( "Not implemented." );
		%	continue;
		%endif
		%
		[ fModelDat, datOut_refresh ] = __refresh( fModelDat.vecYPB, funchF, fModelDat, prm );
		fevalCount += datOut_refresh.fevalCount;
		continue;
		%
		error( "Reached end of main loop without having selected an action." );
	endwhile
	%
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction


function prm = __initPrm( vecX, vecF, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	prm.msgCopious = mygetfield( prm, "msgCopious", verbLev >= VERBLEV__COPIOUS );
	prm.msgProgress = mygetfield( prm, "msgProgress", verbLev >= VERBLEV__PROGRESS );
	prm.msgMain = mygetfield( prm, "msgMain", verbLev >= VERBLEV__MAIN );
	prm.msgNotice = mygetfield( prm, "msgNotice", verbLev >= VERBLEV__NOTICE );
	prm.debugMode = mygetfield( prm, "debugMode", true );
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	fTol = max([ sqrt(eps)*norm(vecF), eps ]);
	prm.fTol = mygetfield( prm, "fTol", fTol );
	prm.iterLimit = mygetfield( prm, "iterLimit", ceil(20 + 10*sqrt(sizeX) + 0.01*sizeX) );
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
	%
	bScale0 = mygetfield( prm, "bScale0", eps );
	assert( 0.0 < bScale0 );
	fModelDat.matB = [ bScale0 ];
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
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify y'*B*y <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	%sizeV = size(matV,2);
	%
	% Basics.
	matWTW = matW'*matW;
	vecMG = -(matW'*vecF);
	%
	vecYIU = __calcStep( matWTW, vecMG, prm );
	vecYIB = __calcBoundStep( matWTW, matB, vecMG, prm );
	vecYPB = __calcBoundStep( matWTW + matA, matB, vecMG, prm );
	%
	if (0)
		msg( __FILE__, __LINE__, "Doing tests on __calcStep() and __calcBoundStep()." );
		matW = [ 1.0, 2.0; 3.0, 6.0 ]
		vecF = [ 1.0; 1.0 ]
		matH = matW'*matW
		vecG = matW'*vecF
		matB = [ 1.0, 0.0; 0.0, 1.0  ]*100.0
		vecYU = __calcStep( matH, -vecG, prm )
		vecGU = vecG + matH*vecYU
		bu = vecYU'*matB*vecYU
		vecYB = __calcBoundStep( matH, matB, -vecG, prm )
		vecGB = vecG + matH*vecYB
		bb = vecYB'*matB*vecYB
		error( "PLEASE SEE DEBUG RESULTS ABOVE." );
	endif
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
	fModelDat.omegaModeAvglPB = sumsq(vecFModelPB)/2.0;
	%
	fModelDat.omegaModelVarIU = vecYIU'*matA*vecYIU;
	fModelDat.omegaModelVarIB = vecYIB'*matA*vecYIB;
	fModelDat.omegaModelVarPB = vecYPB'*matA*vecYPB;
	%
	fModelDat.bIU = vecYIU'*matB*vecYIU;
	fModelDat.bIB = vecYIB'*matB*vecYIB;
	fModelDat.bPB = vecYPB'*matB*vecYPB;
	%
	fModelDat.omega = sumsq(vecF)/2.0;
return;
endfunction


function vecY = __calcStep( matH, vecMG, prm )
	[ matR, cholFlag ] = chol( matH );
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	if ( 0 == cholFlag && min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		vecY = matR \ (matR'\vecMG);
	else
		msgif( prm.msgNotice, __FILE__, __LINE__, "Extrapolating step to singular point." );
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


function vecY = __calcBoundStep( matH, matB, vecMG, prm )
	[ matR, cholFlag ] = chol( matH );
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	if ( 0 == cholFlag && min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		s1 = 1.0;
	else
		epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
		hScale = max(max(abs(matH)));
		bScale = max(max(abs(matB)));
		assert( 0.0 < hScale );
		assert( 0.0 < bScale );
		s1 = 1.0 - (epsRelRegu*hScale/bScale);
	endif
	%
	funchYOfS = @(s)( ( s*matH + (1.0-s)*matB ) \ (s*vecMG) );
	funchBOfY = @(y)( y'*matB*y );
	vecY1 = funchYOfS(s1);
	b1 = funchBOfY(vecY1);
	if ( b1 <= 1.0 )
		vecY = vecY1;
		return;
	endif
	%
	funchBM1OfS = @(s)( funchBOfY(funchYOfS(s)) - 1.0 );
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
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify y'*B*y <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
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
	%
	bScale0 = mygetfield( prm, "bScale0", eps );
	bScale1 = mygetfield( prm, "bScale1", 1.0e-4 );
	assert( 0.0 < bScale0 );
	assert( 0.0 < bScale1 );
	fModelDat.matB = zeros(sizeV+1,sizeV+1);
	fModelDat.matB(1:sizeV,1:sizeV) = matB;
	fModelDat.matB(sizeV+1,sizeV+1) = bScale0 + bScale1*sum(sum(matB.^2))/(sizeV^2);
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
	msg( __FILE__, __LINE__, "TODO: Update matA in a REASONABLE fashion!" );
	matW = fModelDat.matW;
	fModelDat.matA += diag(abs(diag(matW'*matW)));
	msg( __FILE__, __LINE__, "TODO: Update matW when move!" );
	fModelDat.matVLocal = [];
	fModelDat.vecX = vecX_trial;
	fModelDat.vecF = vecF_trial;
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
	%matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify y'*B*y <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
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
	%%%matE = eye(sizeV,sizeV) - vecYHat*(vecYHat');
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



function fModelDat = __increaseB( vecY, fModelDat, prm )
	% Unpack.
	%vecX = fModelDat.vecX;
	%vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%matV = fModelDat.matV; % Subspace basis matrix.
	%matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	%matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify y'*B*y <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	%sizeV = size(matV,2);
	%
	%
	%error( "Frick. Had wrong definition of matB..." );
	msg( __FILE__, __LINE__, "This is a crude hack for __increaseB()." );
	b = vecY'*matB*vecY;
	assert( 0.0 < b );
	matB /= b;
	%
	fModelDat.matB = matB;
	if (0)
		fModelDat.hackCounter = mygetfield( fModelDat, "hackCounter", 0 );
		fModelDat.hackCounter++;
		if ( 10 <= fModelDat.hackCounter )
			error( "BREAK!" );
		endif
	endif
return;
endfunction

