function vecX = mycholdiv( matA, vecB, requireWellBehavedSolution=true, prm=[] )
	sz = size(vecB,1);
	assert( isrealarray(vecB,[sz,1]) );
	b = norm(vecB);
	if ( 0.0 == b )
		msg( __FILE__, __LINE__, "WARNING: vecB is zero." );
	endif
	epsCasual = eps^0.2;
	epsStrict = eps^0.8;
	epsIntermed = sqrt( epsCasual * epsStrict );
	cholTol = eps^0.4;
	%
	assert( isrealarray(matA,[sz,sz]) );
	assert( issymmetric(matA,epsCasual) );
	aScl = max(max(abs(matA)));
	if ( 0.0 == aScl )
		msg( __FILE__, __LINE__, "WARNING: matA is zero." );
		vecX = vecB;
		matR = zeros(sz,sz);
		return;
	endif
	aDiagMin = min(diag(matA));
	%
	if ( aDiagMin > 0.0 )
		msg( __FILE__, __LINE__, "Attempting 'clearly posititive-definite' solution..." );
		[ matR, cholFlag ] = chol( matA );
		if ( 0 == cholFlag )
			if ( min(diag(matR)) > cholTol*sqrt(aScl) )
				vecX = matR \ ( matR' \ vecB );
				msg( __FILE__, __LINE__, "Accepting 'clearly posititive-definite' solution." );
				return;
			endif
			msg( __FILE__, __LINE__, "  Cholesky factorization was unreliable." );
		else
			msg( __FILE__, __LINE__, "  Cholesky factorization failed." );
		endif
		clear matR;
		clear cholFlag;
	endif
	%
	validateExtrapolation = mygetfield( prm, "validateExtrapolation", true );
	epsExtrap = epsIntermed;
	if ( aDiagMin >= 0.0 )
		msg( __FILE__, __LINE__, "Attempting 'extrapolated posititive-semi-definite' solution..." );
		[ matR1, cholFlag ] = chol( matA + epsExtrap*aScl*eye(sz,sz) );
		if ( 0 == cholFlag )
			% Since above chol() did not fail, following chol() should never fail.
			matR2 = chol( matA + 2.0*epsExtrap*aScl*eye(sz,sz) );
			vecX1 = matR1 \ ( matR1' \ vecB );
			vecX2 = matR2 \ ( matR2' \ vecB );
			vecX = ( (2.0*vecX1) - vecX2 );
			if ( validateExtrapolation )
				% We could pretty much skip this if 0.0 == b, but, meh.
				matR3 = chol( matA + 3.0*epsExtrap*aScl*eye(sz,sz) );
				vecX3 = matR3 \ ( matR3' \ vecB );
				vecXAlt = ( (1.5*vecX1) - (0.5*vecX3) );
				if ( reldiff( vecX, vecXAlt ) <= epsCasual )
					msg( __FILE__, __LINE__, "Accepting validated 'extrapolated posititive-semi-definite' solution." );
					return;
				endif
				msg( __FILE__, __LINE__, "  Extrapolation was inconsistent." );
			else
				msg( __FILE__, __LINE__, "Accepting unvalidated 'extrapolated posititive-semi-definite' solution." );
				return;
			endif
		endif
		clear matR1;
		clear cholFlag;
		msg( __FILE__, __LINE__, "  Cholesky factorization failed." );
	endif
	if (requireWellBehavedSolution)
		error( "Failed to find a well-behaved solution." );
	endif
	%
	% We're on to the "has (at least one) negative (eigenvalue)" case.
	msg( __FILE__, __LINE__, "Generating perturbed positive-definite solution." );
	vecPosDefDiagMin = sum(abs(matA),2) + epsStrict*aScl; % Scalar autobroadcast.
	vecADiag = diag(matA);
	vecmaskModifyMe = (vecPosDefDiagMin>vecADiag);
	vecADiagMod = vecADiag;
	vecADiagMod(vecmaskModifyMe) = vecPosDefDiagMin(vecmaskModifyMe);
	matAMod = matA + diag(vecADiagMod - vecADiag);
	%
	% vecAMod should now be strictly positive-definite, beyond an finite precision issues.
	matR = chol( matAMod );
	vecX = matR \ ( matR' \ vecB );
	return;
return;
endfunction
