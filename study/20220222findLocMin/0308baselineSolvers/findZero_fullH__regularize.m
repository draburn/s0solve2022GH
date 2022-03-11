	matHRegu = matH;
	[ matR, cholFlag ] = chol( matHRegu );
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	assert( isrealscalar(cholSafeTol) );
	assert( 0.0 < cholSafeTol  );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
	assert( isrealscalar(epsRelRegu) );
	assert( 0.0 < epsRelRegu  );
	matHRegu = matH + ( epsRelRegu * hNorm * matIX );
	[ matR, cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	msg( __FILE__, __LINE__, "Hessian appears to have a negative eigenvalue." );
	msg( __FILE__, __LINE__, " Performing full eigendecomposition." );
	msg( __FILE__, __LINE__, " This may be slow;" );
	msg( __FILE__, __LINE__, "  future optimization may be possible." );
	%
	[ matPsi_eig, matLambda_eig ] = eig( matH );
	muCrit = -min(diag(matLambda_eig));
	mu = muCrit + epsRelRegu * ( muCrit + hNorm );
	matHRegu = matH + ( mu * matIX );
	[ matR, cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	error( "Regularization of Hessian failed." );
