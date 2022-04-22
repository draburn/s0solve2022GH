% Dev
%  9xx = approaching new JFNK (with phi-patch, etc),
%   for prototyping structure and refereshing memory.
%  900 = JFNK (w strict start) + AP (w OSQU) + Lev*minscan
%  904 = 900 but replace Lev*minscan with BT/TR per _444 and coast.


% 20220422 TODO
%   Confirm phi patch is analytically correct.
%   Put phi-patch in loop with linsolf?
%   Make BT smarter...?

function [ vecXF, vecFF, datOut ] = findZero_904( vecX0, funchF, prm=[] )
	time0 = time();
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = mygetfield( prm, "vecF0", [] );
	if (isempty(vecF0))
		vecF0 = funchF( vecX0 );
		fevalCount++;
	endif
	sizeF = size(vecF0,1);
	assert( 0 < sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( 0~=norm(vecF0) );
	%
	%
	vecX_best = vecX0;
	vecF_best = vecF0;
	%
	%
	%
	%fModelType = mygetfield( prm, "fModelType", F_MODEL_TYPE__CONVENTIONAL );
	%fModelDat.modelType = fModelType;
	fModelDat.vecX = vecX0;
	fModelDat.vecF = vecF0;
	fModelDat.matA = eye(sizeF,sizeX); % Model for full-space Jacobian.
	%fModelDat.matV = zeros(sizeX,0); % Most recent local analysis V.
	%fModelDat.matW = zeros(sizeF,0); % Most recent local analysis JV, plus updates.
	fModelDat = mygetfield( prm, "fModelDat0", fModelDat );
	%assert( fModelDat.modelType == fModelType );
	assert( reldiff(vecX0,fModelDat.vecX,eps) < eps );
	assert( reldiff(vecF0,fModelDat.vecF,eps) < eps );
	assert( isrealarray(fModelDat.matA,[sizeF,sizeX]) );
	%
	%
	%
	iterCount = 0;
	vecX = vecX0;
	vecF = vecF0;
	datOut.iterCountVals(iterCount+1) = iterCount;
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.fNormVals(iterCount+1) = norm(vecF_best);
	datOut.vecXVals(:,iterCount+1) = vecX;
	datOut.vecFVals(:,iterCount+1) = vecF;
	%
	step_tol = sqrt(eps); % Use a tight solve on first iteration to get a large subspace.
	step_prm = mygetfield( prm, "step_prm", [] );
	stepSearchDat = [];
	[ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, stepSearchDat, step_tol, step_prm );
	fevalCount += step_datOut.fevalCount;
	%
	%
	%
	while (1);
		iterCount++;
		if ( norm(vecF_next) < norm(vecF_best) )
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
		datOut.iterCountVals(iterCount+1) = iterCount;
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.fNormVals(iterCount+1) = norm(vecF_best);
		datOut.vecXVals(:,iterCount+1) = vecX_next;
		datOut.vecFVals(:,iterCount+1) = vecF_next;
		%
		msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
		  time()-time0, iterCount, fevalCount, ...
		  norm(vecX_next-vecX0), norm(vecX_next-vecX0)-norm(vecX-vecX0), norm(vecX_next-vecX), ...
		  norm(vecF_next), norm(vecF)-norm(vecF_next), norm(vecF-vecF_next) ) );
		%
		%
		%
		fTol = mygetfield( prm, "fTol", eps );
		iterMax = mygetfield( prm, "iterMax", 50 );
		if ( norm(vecF_best) <= fTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: norm(vecF_best) <= fTol." );
			break;
		elseif ( norm(vecF_next) + fTol >= norm(vecF) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecF_next) + fTol >= norm(vecF)." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		vecX = vecX_next;
		vecF = vecF_next;
		fModelDat = fModelDat_next;
		%
		%
		%
		% TODO: Add coasting.
		% If successful, "continue" back to start of loop.
		% If not, do following...
		%
		%
		%
		step_tol = 0.1*norm(vecF)/norm(vecF0);
		step_prm = mygetfield( prm, "step_prm", [] );
		[ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, stepSearchDat, step_tol, step_prm );
		fevalCount += step_datOut.fevalCount;
	endwhile
	%
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
return;
endfunction



function [ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, stepSearchDat, step_tol, step_prm )
	fevalCount = 0;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	matA = fModelDat.matA;
	%
	%
	%
	epsFD = mygetfield( step_prm, "epsFD", eps^0.4 );
	funchMatJProd = @(v)( ( funchF(vecX+epsFD*v) - vecF ) / epsFD );
	vecXGuess0 = zeros(sizeX,1);
	linsolf_prm = mygetfield( step_prm, "linsolf_prm", [] );
	linsolf_prm.matP = inv(matA);
	linsolf_prm.tol = step_tol;
	%
	[ vecSSDeltaN, linsolf_datOut ] = linsolf( funchMatJProd, -vecF, vecXGuess0, linsolf_prm );
	fevalCount += linsolf_datOut.fevalCount;
	%
	matV = linsolf_datOut.matV;
	sizeV = size(matV,2);
	assert( isrealarray(matV,[sizeX,sizeV]) );
	matW = linsolf_datOut.matW;
	assert( isrealarray(matW,[sizeF,sizeV]) );
	matA += ( matW - (matA*matV) ) * (matV'); % Update per subspace.
	%[ matLambda ] = eig( matW'*matW ); msg( __FILE__, __LINE__, sprintf( "phiMin = %e.", min(matLambda)/max(abs(matLambda)) ) );
	%
	%
	%
	fModelDat_curr = fModelDat;
	fModelDat_curr.matV = matV;
	fModelDat_curr.matW = matW;
	fModelDat_curr.matA = matA;
	%
	% Phi-patch concept was to do *inside* linsolf, but, whatever.
	% Also, concept would allow for multiple patches, but we use only one here.
	usePostLinsolfPhiPatch = mygetfield( step_prm, "usePostLinsolfPhiPatch", true );
	if ( usePostLinsolfPhiPatch )
		fModelDat_curr.vecPhi = [];
		fModelDat_curr.vecGamma = [];
		phiPatchTol = mygetfield( step_prm, "phiPatchTol", 1.0 );
		matWTW = matW'*matW;
		[ matPhi, matLambda ] = eig( matWTW );
		vecLambda = diag(matLambda);
		[ lambdaMin, indexOfLambdaMin ] = min( vecLambda );
		if ( lambdaMin / max(abs(vecLambda)) < phiPatchTol )
			vecPhi = matPhi(:,indexOfLambdaMin);
			phiPatchOrder = mygetfield( step_prm, "phiPatchOrder", 2 );
			switch (phiPatchOrder)
			case 1
				error( "phiPatchOrder 1 not implemented." );
			case 2
				vecXP = vecX + epsFD*matV*vecPhi;
				vecXM = vecX - epsFD*matV*vecPhi;
				vecFP = funchF( vecXP ); fevalCount++;
				vecFM = funchF( vecXM ); fevalCount++;
				vecGamma = ( vecFP + vecFM - 2.0*vecF ) / (epsFD^2);
			otherwise
				error( "Invalid value of phiPatchOrder." );
			endswitch
			fModelDat_curr.vecPhi = vecPhi;
			fModelDat_curr.vecGamma = vecGamma;
		endif
	endif
	%
	stepSearchDat_curr = stepSearchDat;
	stepSearch_prm = [];
	stepSearch_prm.usePostLinsolfPhiPatch = usePostLinsolfPhiPatch;
	%
	[ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, stepSearch_datOut ] = __searchForStep_tr444( ...
	  funchF, vecX, vecF, fModelDat_curr, stepSearchDat_curr, stepSearch_prm );
	fevalCount += stepSearch_datOut.fevalCount;
	%
	%
	%
	step_datOut.fevalCount = fevalCount;
	step_datOut.linsolf_datOut = linsolf_datOut;
	step_datOut.sizeV = sizeV;
	step_datOut.stepSearch_datOut = stepSearch_datOut;
endfunction



function [ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, stepSearch_datOut ] = __searchForStep_minscan( ...
	  funchF, vecX, vecF, fModelDat, stepSearchDat, stepSearch_prm )
	%
	fevalCount = 0;
	matV = fModelDat.matV;
	matW = fModelDat.matW;
	matA = fModelDat.matA;
	sizeV = size(matV,2);
	%
	vecG = matW'*vecF;
	matHRaw = matW'*matW;
	usePostLinsolfPhiPatch = mygetfield( stepSearch_prm, "usePostLinsolfPhiPatch", true );
	if ( usePostLinsolfPhiPatch )
		matHRaw += (vecF'*fModelDat.vecGamma) * fModelDat.vecPhi * ( fModelDat.vecPhi' );
	endif
	matHRegu = calcHRegu(matHRaw);
	%
	matIV = eye(sizeV,sizeV);
	funchYOfP = @(p)( ( p*matHRegu + (1.0-p)*matIV ) \ (-p*vecG) );
	funchDeltaOfP = @(p)( matV * funchYOfP(p) );
	%
	funchFNormOfP = @(p)( norm(funchF(vecX+funchDeltaOfP(p))) );
	fminbnd_options = optimset( "TolX", 1.0E-4, "TolFun", norm(vecF)*1.0E-2 );
	fminbnd_options = mygetfield( stepSearch_prm, "fminbnd_options", fminbnd_options );
	%
	[ fminbnd_x, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchFNormOfP, 0.0, 1.0, fminbnd_options );
	fevalCount += fminbnd_output.funcCount;
	p = fminbnd_x;
	%
	vecX_next = vecX + funchDeltaOfP(p);
	vecF_next = funchF(vecX_next);
	fevalCount++;
	%
	fooX = vecX_next - vecX;
	assert( norm(fooX) > 0.0 );
	fooF = vecF_next - ( vecF + matA*fooX );
	matA += 2.0 * fooF * (fooX') / (fooX'*fooX); % Quadratic upate.
	fModelDat_next = fModelDat;
	fModelDat_next.matA = matA;
	%
	stepSearchDat_next = [];
	stepSearch_datOut.fevalCount = fevalCount;
endfunction



function [ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, stepSearch_datOut ] = __searchForStep_tr444( ...
	  funchF, vecX, vecF, fModelDat, stepSearchDat, stepSearch_prm )
	%
	setVerbLevs;
	verbLev = mygetfield( stepSearch_prm, "verbLev", VERBLEV__MAIN );
	fevalCount = 0;
	matV = fModelDat.matV;
	matW = fModelDat.matW;
	matA = fModelDat.matA;
	sizeV = size(matV,2);
	%
	vecG = matW'*vecF;
	matHRaw = matW'*matW;
	usePostLinsolfPhiPatch = mygetfield( stepSearch_prm, "usePostLinsolfPhiPatch", true );
	if ( usePostLinsolfPhiPatch )
	if (~isempty(fModelDat.vecGamma))
		matHRaw += (vecF'*fModelDat.vecGamma) * fModelDat.vecPhi * ( fModelDat.vecPhi' );
	endif
	endif
	matHRegu = calcHRegu(matHRaw);
	%
	matIV = eye(sizeV,sizeV);
	funchYOfP = @(p)( ( p*matHRegu + (1.0-p)*matIV ) \ (-p*vecG) );
	funchDeltaOfP = @(p)( matV * funchYOfP(p) );
	%
	dTreg = mygetfield( stepSearchDat, "dTreg", [] );
	if (isempty(dTreg))
		pMax = 1.0;
	else
		pMax = __findPOfDeltaNorm( dTreg, funchDeltaOfP  );
	endif
	vecY_pMax = funchYOfP( pMax );
	vecFModel_pMax = vecF + matW*vecY_pMax;
	%
	deltaNormTol = mygetfield( stepSearch_prm, "deltaNormTol", 100.0*eps );
	btMax = mygetfield( stepSearch_prm, "btMax", 30 );
	btCount = 0;
	vecDelta_rejected = [];
	p = pMax;
	while (1)
		vecY = funchYOfP( p );
		vecDelta = matV*vecY;
		vecX_next = vecX + vecDelta;
		vecF_next = funchF( vecX_next );
		fevalCount++;
		%
		if ( norm(vecF_next) < 0.5*norm(vecF) + 0.5*norm(vecFModel_pMax) )
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
	vecFModel_next = vecF + matW*vecY;
	rhoThresh0 = mygetfield( stepSearch_prm, "rhoThresh0", 0.05 );
	rhoThresh1 = mygetfield( stepSearch_prm, "rhoThresh1", 0.30 );
	rho = norm(vecF_next-vecFModel_next)/norm(vecF);
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
	%
	%
	fooX = vecX_next - vecX;
	assert( norm(fooX) > 0.0 );
	fooF = vecF_next - ( vecF + matA*fooX );
	matA += 2.0 * fooF * (fooX') / (fooX'*fooX); % Quadratic upate.
	fModelDat_next = fModelDat;
	fModelDat_next.matA = matA;
	%
	stepSearchDat_next = [];
	stepSearch_datOut.fevalCount = fevalCount;
endfunction



function p = __findPOfDeltaNorm( deltaNormMax, funchDeltaOfP )
	fzero_fun = @(p_dummy)( norm(funchDeltaOfP(p_dummy)) - deltaNormMax );
	if ( fzero_fun(1.0) <= 0.0 )
		p = 1.0;
		return;
	endif
	p = fzero( fzero_fun, [0.0, 1.0] );
return;
endfunction
