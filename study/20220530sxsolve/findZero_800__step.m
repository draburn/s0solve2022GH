	%
	matH = matW'*matW;
	vecG = matW'*vecF;
	hNorm = sqrt(sum(sumsq(matH))/sizeX);
	assert( 0 ~= hNorm );
	cholTol = mygetfield( prm, "cholTol", 1e-6 );
	[ matR, cholFlag ] = chol(matH);
	%if ( 0==cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) ) % But "\" warning triggers off rcond!
	if ( rcond(matH) > eps^0.8 && 0==cholFlag )
		matHRegu = matH;
	else
		sizeV = size(matW,2);
		matHRegu = matH + hNorm*sqrt(eps)*eye(sizeV);
		[ matR, cholFlag ] = chol(matHRegu);
		assert( 0 == cholFlag );
		assert( min(diag(matR)) > max(abs(diag(matR)))*cholTol );
	endif
	%
	funchDeltaOfP = @(p) (matV*( ( p*matHRegu + (1.0-p)*eye(size(matHRegu)) ) \ (-p*vecG) )); % <<< THIS WORKS.
	%%%funchDeltaOfP = @(p) ( matV * (__funcSSDeltaOfP( p, matHRegu, vecG )) ); %%% <<< THIS DOES NOT WORK!
	% DRaburn 2022.05.10: Shouldn't that be hScale*eye?
	%
	%
	pMax = __findPOfDeltaNorm( dTreg, funchDeltaOfP  );
	vecY_pMax = __funcSSDeltaOfP( pMax, matHRegu, vecG );
	vecDelta_pMax = matV*vecY_pMax;
	vecFModel_pMax = vecF + matW*vecY_pMax;
	%
	%
	if ( isempty(initialFallRatio) )
		% This step is using freshly calculated JV.
	else
		% This step is using possibly stale JV.
		% If model says too little merit, then bail.
		fallRatioModelExp = mygetfield( prm, "fallRatioModelExp", 0.5 );
		if ( norm(vecFModel_pMax)/norm(vecF) > initialFallRatio^fallRatioModelExp )
			vecX_next = vecX;
			vecF_next = vecF;
			return;
		endif
		%msgif( verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "Coasting!" );
	endif
	assert( norm(vecFModel_pMax) <= norm(vecF) );
	%
	%
	deltaNormTol = mygetfield( prm, "deltaNormTol", sizeX*100.0*eps );
	btMax = mygetfield( prm, "btMax", 3 );
	btCount = 0;
	vecDelta_rejected = [];
	p = pMax;
	while (1)
		vecY = __funcSSDeltaOfP( p, matHRegu, vecG );
		vecDelta = matV*vecY;
		vecX_next = vecX + vecDelta;
		vecF_next = funchF(vecX_next); fevalCount++;
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
	vecFModel_next = vecF + matW*vecY;
	%
	rhoThresh0 = mygetfield( prm, "rhoThresh0", 0.05 );
	rhoThresh1 = mygetfield( prm, "rhoThresh1", 0.30 );
	rho = norm(vecF_next-vecFModel_next)/norm(vecF);
	%
	if ( ~isempty(vecDelta_rejected) )
		dTreg = min([ dTreg, norm(vecDelta_rejected) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  step: Have a rejected step. Set dTreg = %0.3e.", dTreg ) )
	endif
	if ( rho < rhoThresh0 )
		% Model is very accurate at the point.
		dTreg = max([ dTreg, 2.0*norm(vecDelta) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  step: Model was very accurate. Set dTreg = %0.3e.", dTreg ) )
	elseif ( rho > rhoThresh1 )
		% Model was inaccurate at the point.
		dTreg = min([ dTreg, norm(vecDelta) ]); % "min([ dTreg," should be superfluous.
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  step: Model was very inaccurate. Set dTreg = %0.3e.", dTreg ) )
	endif
