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
	haveCand = false; % Have candidate for next X?
	vecXCand = [];
	vecFCand = [];
	vecX_best = vecX_initial;
	vecF_best = vecF_initial;
	%
	%
	while (1)
		iterCount++;
		fModelDat = __analyzeModel( fModelDat, prm );
		%
		% Simple stoping criteria.
		if ( sumsq(vecF_best) <= prm.fTol^2 )
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
			[ fModelDat, datOut_expandModel ] = __expandModel( fModelDat.vecFModelIU, funchF, fModelDat, prm );
			fevalCount += datOut_expandModel.fevalCount;
			continue;
		endif
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
	fModelDat.matVLoc = [ vecV ];
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
	%vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify y'*B*y <= 1.
	%matVLoc = fModelDat.matVLoc;
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
	matVLoc = fModelDat.matVLoc; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify y'*B*y <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%
	uNorm = norm(vecU);
	assert( 0.0 < uNorm );
	vecV = __calcOrthonorm( vecU, fModelDat.matV, prm );
	if ( 0.0 == norm(vecV) )
		error( "Failed to generate new subspace basis vector." );
	endif
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecV) )
		error( "Jacobian along new subspace basis vector is zero." );
	endif
	%
	fModelDat.matVLoc = [ matVLoc, vecV ];
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
