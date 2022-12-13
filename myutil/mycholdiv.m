function vecX = mycholdiv( matA, vecB, forceSolution=false, prm=[] )
	sz = size(vecB,1);
	assert( isrealarray(vecB,[sz,1]) );
	b = norm(vecB);
	if ( 0.0 == b )
		msg( __FILE__, __LINE__, "WARNING: vecB is zero." );
	endif
	epsCasual = eps^0.2;
	epsStrict = eps^0.8;
	epsIntermed = sqrt( epsCasual * epsStrict );
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
		msg( __FILE__, __LINE__, "Attempting 'clearly posititive-definite' solution." );
		[ matR, cholFlag ] = chol( matA );
		if ( 0 == cholFlag )
		if ( min(diag(matR)) > epsCasual*sqrt(aScl) ) % Agressively distrust this solution.
			vecX = matR \ ( matR' \ vecB );
			%min(diag(matR))
			%epsCasual*sqrt(aScl)
			%assert( reldiff( matA * vecX, vecB ) < epsCasual );
			msg( __FILE__, __LINE__, "Accepting 'clearly posititive-definite' solution." );
			return;
		endif
		endif
		clear matR;
		clear cholFlag;
	endif
	msg( __FILE__, __LINE__, "Rejected 'clearly posititive-definite' solution." );
	%
	validateExtrapolation = mygetfield( prm, "validateExtrapolation", true );
	epsExtrap = epsIntermed;
	if ( aDiagMin >= 0.0 )
		msg( __FILE__, __LINE__, "Attempting 'extrapolated posititive-semi-definite' solution." );
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
				if ( reldiff( vecX, vecXAlt ) > epsCasual )
					if (forceSolution)
						msg( __FILE__, __LINE__, "WARNING: positive-semi-definite extrapolation is inconsistent." );
						msg( __FILE__, __LINE__, "  ( This may be because gradient is not in span of Hessian. )" );
					else
						msg( __FILE__, __LINE__, "ERROR: positive-semi-definite extrapolation is inconsistent." );
						msg( __FILE__, __LINE__, "  ( This may be because gradient is not in span of Hessian. )" );
						error("Positive-semi-definite extrapolation is inconsistent.");
					endif
				endif
			endif
			msg( __FILE__, __LINE__, "Accepting 'extrapolated posititive-semi-definite' solution." );
			return;
		endif
		clear matR1;
		clear cholFlag;
	endif
	msg( __FILE__, __LINE__, "Rejected 'extrapolated posititive-semi-definite' solution." );
	%
	% We're on to the "has (at least one) negative (eigenvalue)" case.
	if (forceSolution)
		msg( __FILE__, __LINE__, "WARNING: matA appears to have a negative eigenvalue." );
	else
		msg( __FILE__, __LINE__, "ERROR: matA appears to have a negative eigenvalue." );
		error( "matA appears to have a negative eigenvalue." );
	endif
	%
	vecPosDefDiagMin = sum(abs(matA),2) + epsStrict*aScl
	; % Scalar autobroadcast.
	vecADiag = diag(matA);
	vecmaskModifyMe = (vecPosDefDiagMin>vecADiag);
	vecADiagMod = vecADiag;
	vecADiagMod(vecmaskModifyMe) = vecPosDefDiagMin(vecmaskModifyMe);
	matAMod = matA + diag(vecADiagMod - vecADiag);
	%
	% vecAMod should now be strictly positive-definite, beyond an finite precision issues.
	matR = chol( matAMod );
	vecX = matR \ ( matR' \ vecB );
	%vecX = matA \ vecB;
	return;
return;
endfunction
