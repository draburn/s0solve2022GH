% Dev
%  9xx = approaching new JFNK (with phi-patch, etc),
%   for prototyping structure and refereshing memory.
%  900 = JFNK (w strict start) + AP (w OSQU) + Lev*minscan
%  904 = 900 but replace Lev*minscan with BT/TR per _444 and coast.

function [ vecXF, vecFF, datOut ] = findZero_900( vecX0, funchF, prm=[] )
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
	[ vecX_next, vecF_next, fModelDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, step_tol, step_prm );
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
		[ vecX_next, vecF_next, fModelDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, step_tol, step_prm );
		fevalCount += step_datOut.fevalCount;
	endwhile
	%
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
return;
endfunction



function [ vecX_next, vecF_next, fModelDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, step_tol, step_prm )
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
	%
	%
	%
	vecG = matW'*vecF;
	matH = matW'*matW;
	[ foo, cholFlag ] = chol(matH);
	assert( 0 == cholFlag ); % Singular matH not supported here.
	matIV = eye(sizeV,sizeV);
	funchYOfP = @(p)( ( p*matH + (1.0-p)*matIV ) \ (-p*vecG) );
	funchDeltaOfP = @(p)( matV * funchYOfP(p) );
	%
	funchFNormOfP = @(p)( norm(funchF(vecX+funchDeltaOfP(p))) );
	fminbnd_options = optimset( "TolX", 1.0E-4, "TolFun", norm(vecF)*1.0E-2 );
	fminbnd_options = mygetfield( step_prm, "fminbnd_options", fminbnd_options );
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
	fModelDat_next.matA = matA;
	%
	%
	%
	step_datOut.fevalCount = fevalCount;
	step_datOut.linsolf_datOut = linsolf_datOut;
	step_datOut.sizeV = sizeV;
endfunction
