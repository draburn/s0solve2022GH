function [ vecX, matR ] = mycholdiv( matA, vecB, prm=[] )
	assert( 2 <= nargin );
	assert( nargin <= 3 );
	sz = size(vecB,1);
	assert( isrealarray(vecB,[sz,1]) );
	assert( isrealarray(matA,[sz,sz]) );
	debugMode = mygetfield( prm, "debugMode", false );
	symTol = mygetfield( prm, "symTol", sqrt(eps) );
	issymmetric( matA, symTol );
	%
	epsChol = mygetfield( prm, "epsChol", sqrt(eps)*10.0 );
	[ matR, cholFlag ] = chol(matA);
	if ( 0 == cholFlag && min(diag(matR)) > epsChol*max(abs(diag(matR))) )
		if ( debugMode )
			msg( __FILE__, __LINE__, "Original Cholesky factorization was successful." );
			cholFlag
			min(diag(matR))
			max(abs(diag(matR)))
		endif
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	if ( debugMode )
		msg( __FILE__, __LINE__, "Original Cholesky factorization failed." );
		cholFlag
		min(diag(matR))
		max(abs(diag(matR)))
	endif
	clear matR;
	clear cholFlag;
	%
	% matA may be pos-semi-def. It could also be "has-neg".
	epsExtrap = mygetfield( prm, "epsExtrap", sqrt(eps) );
	matE = mygetfield( prm, "matE", eye(sz,sz) );
	aScl = max(max(abs(matA)));
	if ( 0.0 == aScl )
		if ( debugMode )
			msg( __FILE__, __LINE__, "0.0 == aScl. There is clearly no well-behaved solution." );
		endif
		vecX = [];
		matR = [];
		return;
	endif
	%
	[ matR1, cholFlag ] = chol( matA + epsExtrap*aScl*matE );
	if ( 0~= cholFlag )
		if ( debugMode )
			msg( __FILE__, __LINE__, "Extrapolation failed; matrix appears to indefinite/negative." );
			size_matR1 = size(matR1)
			cholFlag
		endif
		vecX = [];
		matR = [];
		return;
	endif
	matR2 = chol( matA + 2.0*epsExtrap*aScl*matE ); % This should never fail; so, if it somehow does, just throw the error.
	vecX1 = matR1 \ ( matR1' \ vecB );
	vecX2 = matR2 \ ( matR2' \ vecB );
	vecX = (2.0*vecX1) - vecX2;
	%
	extrapTol = mygetfield( prm, "extrapTol", sqrt(sqrt(eps)) );
	rd = reldiff( vecX, vecX2 );
	if ( rd > extrapTol )
		if ( debugMode )
			msg( __FILE__, __LINE__, "Extrapolated solution falied consistency check." );
			reldiff_vecX_vecX2 = rd
			extrapTol
		endif
		vecX = [];
		matR = [];
		return;
	endif
	if ( debugMode )
		msg( __FILE__, __LINE__, "Extrapolated solution passed consistency check." );
		reldiff_vecX_vecX2 = rd
		extrapTol
	endif
	%
	matR = (2.0*matR1) - matR2;
	return;
return;
endfunction

%!test
%!	setprngstates(85023360);
%!	sizeX = 2 + ceil(20*rand())
%!	%sizeX = 5
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
%!	vecXCalc = mycholdiv( matA, vecB );
%!	vecBCalc = matA * vecXCalc;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	assert( reldiff(vecX,vecXCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Weakly positive-definite test..." );
%!	matA = matA1 + epsA*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB );
%!	vecBCalc = matA * vecXCalc;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	assert( reldiff(vecX,vecXCalc) < eps^0.1 );
%!	%
%!	msg( __FILE__, __LINE__, "Convergent positive-semi-definite test..." );
%!	matA = matA1 + 0.0*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB );
%!	vecBCalc = matA * vecXCalc;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Non-convergent positive-semi-definite test..." );
%!	matA = matA1 + 0.0*matA2;
%!	vecB = matA * vecX + 0.1*vecBMod;
%!	vecXCalc = mycholdiv( matA, vecB );
%!	if ( ~isempty(vecXCalc) )
%!		msg( __FILE__, __LINE__, "Expected vecXCalc to be empty, but it's not! This should not be typical." );
%!	endif
%!	%
%!	msg( __FILE__, __LINE__, "Strongly 'has negative' test..." );
%!	matA = -( matA1 + matA2 );
%!	vecB = matA * vecX;
%!	vecXCalc = mycholdiv( matA, vecB, false );
%!	assert( isempty(vecXCalc) );
