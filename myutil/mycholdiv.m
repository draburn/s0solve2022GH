function vecX = mycholdiv( matA, vecB, requireWellBehavedSolution=true, prm=[] )
	warning( "This code is deprecated; we should simply use eig() if the first chol() is no good." );
	debugMode = mygetfield( prm, "debugMode", false );
	sz = size(vecB,1);
	if ( debugMode )
		assert( isrealarray(vecB,[sz,1]) );
		assert( isrealarray(matA,[sz,sz]) );
		assert( issymmetric(matA,1.0e-4) );
	endif
	b = norm(vecB);
	if ( 0.0 == b )
		msg( __FILE__, __LINE__, "WARNING: vecB is zero." );
	endif
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
		msgif( debugMode, __FILE__, __LINE__, "Attempting 'clearly posititive-definite' solution..." );
		[ matR, cholFlag ] = chol( matA );
		if ( 0 == cholFlag )
			cholTol = mygetfield( prm, "cholTol", 1.0e-4 );
			if ( min(diag(matR)) > cholTol*sqrt(aScl) ) % We'll be rather strict here.
				vecX = matR \ ( matR' \ vecB );
				msgif( debugMode,__FILE__, __LINE__, sprintf( "  Accepting solution ( %0.3e vs %0.3e ).", min(diag(matR)), cholTol*sqrt(aScl) ) );
				return;
			endif
			msgif( debugMode,__FILE__, __LINE__, sprintf( "  Cholesky factorization was unreliable ( %0.3e vs %0.3e ).", min(diag(matR)), cholTol*sqrt(aScl) ) );
		else
			msgif( debugMode,__FILE__, __LINE__, sprintf( "  Cholesky factorization failed ( %d ).", cholFlag ) );
		endif
		clear matR;
		clear cholFlag;
	endif
	%
	if ( aDiagMin >= 0.0 )
		epsExtrap = mygetfield( prm, "epsExtrap", 1.0e-12 );
		msgif( debugMode,__FILE__, __LINE__, "Attempting 'extrapolated posititive-semi-definite' solution..." );
		[ matR1, cholFlag ] = chol( matA + epsExtrap*aScl*eye(sz,sz) );
		if ( 0 == cholFlag )
			% Since above chol() did not fail, following chol() should never fail.
			matR2 = chol( matA + 2.0*epsExtrap*aScl*eye(sz,sz) );
			vecX1 = matR1 \ ( matR1' \ vecB );
			vecX2 = matR2 \ ( matR2' \ vecB );
			vecX = ( (2.0*vecX1) - vecX2 );
			validateExtrapolation = mygetfield( prm, "validateExtrapolation", true );
			if ( validateExtrapolation )
				% We could pretty much skip this if 0.0 == b, but, meh.
				extrapTol = mygetfield( prm, "extrapTol", 0.1 );
				matR3 = chol( matA + 3.0*epsExtrap*aScl*eye(sz,sz) );
				vecX3 = matR3 \ ( matR3' \ vecB );
				vecXAlt = ( (1.5*vecX1) - (0.5*vecX3) );
				if ( reldiff( vecX, vecXAlt ) <= extrapTol ) % We'll be very loose here.
					msgif( debugMode,__FILE__, __LINE__, sprintf( "  Accepting validated solution ( %0.3e ).", reldiff( vecX, vecXAlt ) ) );
					return;
				else
					msgif( debugMode,__FILE__, __LINE__, sprintf( "  Extrapolation was inconsistent ( %0.3e ).", reldiff( vecX, vecXAlt ) ) );
				endif
				clear matR3;
				clear vecX3;
				clear vecXAlt;
			else
				msgif( debugMode,__FILE__, __LINE__, "  Accepting unvalidated solution." );
				return;
			endif
			clear matR2;
			clear vecX1;
			clear vecX2;
			clear vecX;
		else
			msgif( debugMode,__FILE__, __LINE__, sprintf( "  Cholesky factorization failed ( %d ).", cholFlag ) );
		endif
		clear matR1;
		clear cholFlag;
	endif
	if (requireWellBehavedSolution)
		%dump__matA = matA
		%dump__vecB = vecB
		%dump__aDiagMin = aDiagMin
		%dump__aScl = aScl
		%dump__matR1 = matR1
		%dump__matR2 = matR2
		%dump__matR3 = matR3
		%dump__vecX_vecXAlt = [ vecX, vecXAlt ]
		%dump__rd = reldiff( vecX, vecXAlt )
		error( "Failed to find a well-behaved solution." );
	endif
	%
	% We're on to the "has (at least one) negative (eigenvalue)" case.
	msgif( debugMode,__FILE__, __LINE__, "Generating perturbed positive-definite solution." );
	epsPerturb = mygetfield( prm, "epsPerturb", 1.0e-8 );
	vecPosDefDiagMin = sum(abs(matA),2) - diag(abs(matA)) + epsPerturb*sz*aScl; % Scalar autobroadcast.
	vecADiag = diag(matA);
	vecADiagMod = vecADiag;
	vecADiagMod( vecADiag < vecPosDefDiagMin ) = vecPosDefDiagMin( vecADiag < vecPosDefDiagMin );
	matAMod = matA + diag( vecADiagMod - vecADiag );
	%
	% vecAMod should now be strictly positive-definite, beyond an finite precision issues.
	matR = chol( matAMod );
	vecX = matR \ ( matR' \ vecB );
	return;
return;
endfunction

%!test
%!	setprngstates();
%!	sizeX = 2 + ceil(20*rand());
%!	size1 = sizeX-1;
%!	size2 = 1;
%!	foo = randn(size1,sizeX); matA1 = foo'*foo; clear foo;
%!	foo = randn(size2,sizeX); matA2 = foo'*foo; clear foo;
%!	vecX = randn(sizeX,1);
%!	vecBMod = randn(sizeX,1);
%!	epsA = 1.0e-6;
%!	%
%!	msg( __FILE__, __LINE__, "Strongly positive-definite test..." );
%!	matA = matA1 + 1.0*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB, true );
%!	vecBCalc = matA * vecX;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	assert( reldiff(vecX,vecXCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Weakly positive-definite test..." );
%!	matA = matA1 + epsA*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB, true );
%!	vecBCalc = matA * vecX;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	assert( reldiff(vecX,vecXCalc) < eps^0.1 );
%!	%
%!	msg( __FILE__, __LINE__, "Convergent positive-semi-definite test..." );
%!	matA = matA1 + 0.0*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB, true );
%!	vecBCalc = matA * vecX;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Non-convergent positive-semi-definite test with and without validation..." );
%!	matA = matA1 + 0.0*matA2;
%!	vecB = matA * vecX + 0.1*vecBMod;
%!	vecXCalc = mycholdiv( matA, vecB, false );
%!	prm = [];
%!	prm.validateExtrapolation = false;
%!	vecXCalc = mycholdiv( matA, vecB, false, prm );
%!	%
%!	msg( __FILE__, __LINE__, "Marginal 'has negative' test..." );
%!	matA = matA1 - epsA*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB, false );
%!	%
%!	msg( __FILE__, __LINE__, "Strongly 'has negative' test..." );
%!	matA = -( matA1 + matA2 );
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB, false );
%!	%
%!	msg( __FILE__, __LINE__, "Note: Testing of requireWellBehavedSolution has not been implemented." );
