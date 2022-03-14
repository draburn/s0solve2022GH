	matH_input = matHScaled;
	%
	matHRegu = matH_input;
	[ matR_hRegu, cholFlag ] = chol( matHRegu );
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	assert( isrealscalar(cholSafeTol) );
	assert( 0.0 < cholSafeTol  );
	if ( 0 == cholFlag )
	if ( min(diag(matR_hRegu)) > cholSafeTol*max(abs(diag(matR_hRegu))) )
		return;
	endif
	endif
	%
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Hessian is not positive definite; attempting epsilon regularization..." );
	%
	epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
	assert( isrealscalar(epsRelRegu) );
	assert( 0.0 < epsRelRegu  );
	hNorm_input = sqrt( sum(sumsq(matH_input))/sizeX );
	matHRegu = matH_input + ( epsRelRegu * hNorm_input * matIX );
	[ matR_hRegu, cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR_hRegu)) > cholSafeTol*max(abs(diag(matR_hRegu))) )
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Epsilon regularization successful." );
		return;
	endif
	endif
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Epsilon regularization failed." );
	%
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Hessian appears to have a negative eigenvalue." );
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, " Performing full eigendecomposition." );
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, " This may be slow; run-time faster approaches may be possible." );
	%
	[ matPsi_eig, matLambda_eig ] = eig( matH_input );
	muCrit = -min(diag(matLambda_eig));
	mu = muCrit + epsRelRegu * ( muCrit + hNorm_input );
	matHRegu = matH_input + ( mu * matIX );
	[ matR_hRegu, cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR_hRegu)) > cholSafeTol*max(abs(diag(matR_hRegu))) )
		return;
	endif
	endif
	%
	error( "Regularization of Hessian failed." );
