	% Consider letting an inital guess for mu be passed in.
	mu = 0.0;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with mu zero failed; trying with a small mu." );
	mu = muReguCoeff*hMaxAbs;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with small mu failed; finding muCrit." );
	msgif( debugMode, __FILE__, __LINE__, "Calling eig(). This may be slow. Faster approaches may be possible, such as:" );
	msgif( debugMode, __FILE__, __LINE__, "  start with mu = upper bound for eigenvalue of H, and target omega = 0.0; or," );
	msgif( debugMode, __FILE__, __LINE__, "  increase mu exponentially until chol() works." );
	[ matPsi_eig, matLambda_eig ] = eig( matH );
	muCrit = -min(diag(matLambda_eig));
	mu = muCrit + muReguCoeff * ( muCrit + hNorm );
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( debugMode )
		clear matPsi_eig;
		clear matLambda_eig;
		clear muCrit;
	endif
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
		return;
	endif
	endif
	error( "Cholesky factorization failed even for a very large mu; this should be impossible!" );
