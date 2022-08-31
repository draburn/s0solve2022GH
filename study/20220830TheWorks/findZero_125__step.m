	matH = matJ'*matJ;
	vecG = matJ'*vecF;
	hNorm = sqrt(sum(sumsq(matH))/sizeX);
	assert( 0 ~= hNorm );
	[ matR, cholFlag ] = chol(matH);
	if ( 0==cholFlag )
		matHRegu = matH;
	else
		matHRegu = matH + hNorm*sqrt(eps)*matIX;
		matR = chol(matHRegu); % Ensure it's pos-def.
	endif
	%
	funchDeltaOfP = @(p) ( p*matHRegu + (1.0-p)*matIX ) \ (-p*vecG);
	%
	pMax = __findPOfDeltaNorm( dTreg, funchDeltaOfP );
	vecDelta_pMax = funchDeltaOfP(pMax);
	vecFModel_pMax = vecF + matJ*vecDelta_pMax;
	assert( norm(vecFModel_pMax) < norm(vecF) );
	%
	funchFNormOfP = @(p) norm(funchF(vecX+funchDeltaOfP(p)));
	fminbnd_options = optimset( "TolX", 1.0E-3, "TolFun", norm(vecF)*1.0E-4 );
	[ fminbnd_x, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchFNormOfP, 0.0, pMax, fminbnd_options );
	fevalCount += fminbnd_output.funcCount;
	p = fminbnd_x;
	%
	vecDelta = funchDeltaOfP(p);
	vecFModel_next = vecF + matJ*vecDelta;
	vecX_next = vecX + vecDelta;
	vecF_next = funchF(vecX_next); fevalCount++;
	%
	rhoThresh0 = mygetfield( prm, "rhoThresh0", 0.05 );
	rhoThresh1 = mygetfield( prm, "rhoThresh1", 0.30 );
	rho = norm(vecF_next-vecFModel_next)/norm(vecF);
	%
	if ( p + sqrt(eps) < pMax );
		% We know the full step was rejected.
		dTreg = min([ dTreg, norm(vecDelta_pMax) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "pMax was rejected. Set dTreg = %0.3e.", dTreg ) )
	endif
	if ( rho < rhoThresh0 )
		% Model is very accurate at the point.
		dTreg = max([ dTreg, 2.0*norm(vecDelta) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Model was very accurate. Set dTreg = %0.3e.", dTreg ) )
	elseif ( rho > rhoThresh1 )
		% Model was inaccurate at the point.
		dTreg = min([ dTreg, norm(vecDelta) ]);
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Model was very inaccurate. Set dTreg = %0.3e.", dTreg ) )
	endif
