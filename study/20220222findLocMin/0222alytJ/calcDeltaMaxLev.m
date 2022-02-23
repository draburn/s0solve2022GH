% Function...
%  [ vecDelta, datOut ] = calcDeltaMaxLev( vecG, matH, prm=[] )

function [ vecDelta, datOut ] = calcDeltaMaxLev( vecG, matH, prm=[] )
	%
	%
	sizeX = size(vecG,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if ( debugMode )
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG,[sizeX,1]) );
		assert( isrealarray(matH,[sizeX,sizeX]) );
		assert( issymmetric(matH) );
	endif
	hNorm = sqrt(sum(sumsq(matH)));
	matI = eye(sizeX,sizeX);
	mu0 = 0.0;
	muReguCoeff = 1.0e-5;
	if ( ~isempty(prm) )
		mu0 = mygetfield( prm, "mu0", mu0 );
		muReguCoeff = mygetfield( prm, "muReguCoeff", muReguCoeff );
	endif
	if ( debugMode )
		assert( 0~=hNorm );
		assert( isrealscalar(mu0) );
		assert( isrealscalar(muReguCoeff) );
		assert( 0.0 < muReguCoeff );
	endif
	if ( nargout >= 2 )
		datOut = [];
	endif
	% Not bothering to check if vecG = 0; also, user could possibly want mu and matR.
	%
	%
	mu = mu0;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
		vecDelta = -( matR \ (matR'\vecG) );
		if ( nargout >= 2 )
			datOut.mu = mu;
			datOut.matR = matR;
		endif
		return;
	endif
	endif
	msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with mu0 failed." );
	%
	%
	if ( mu0 < 0.0 )
		mu = 0.0;
		matM = matH + mu*matI;
		[ matR, cholFlag ] = chol( matM );
		if ( 0 == cholFlag )
		if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
			vecDelta = -( matR \ (matR'\vecG) );
			if ( nargout >= 2 )
				datOut.mu = mu;
				datOut.matR = matR;
			endif
			return;
		endif
		endif
		msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with mu = 0.0 failed." );
	endif
	%
	%
	mu = max([ 0.0, mu0 ]) + muReguCoeff*hMaxAbs;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
		vecDelta = -( matR \ (matR'\vecG) );
		if ( nargout >= 2 )
			datOut.mu = mu;
			datOut.matR = matR;
		endif
		return;
	endif
	endif
	msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with 'slightly larger' my failed." );
	%
	%
	msgif( debugMode, __FILE__, __LINE__, "Finding muCrit using eig()." );
	msgif( debugMode, __FILE__, __LINE__, "This may be slow. Faster approaches may be possible, such as:" );
	msgif( debugMode, __FILE__, __LINE__, "  start with mu = upper bound for eigenvalue of H, and target omega = 0.0; or," );
	msgif( debugMode, __FILE__, __LINE__, "  increase mu exponentially until chol() works." );
	[ matPsi_eig, matLambda_eig ] = eig( matH );
	muCrit = -min(diag(matLambda_eig));
	mu = muCrit + muReguCoeff * ( muCrit + hNorm );
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > sqrt(eps)*max(abs(diag(matR))) )
		vecDelta = -( matR \ (matR'\vecG) );
		if ( nargout >= 2 )
			datOut.mu = mu;
			datOut.matR = matR;
		endif
		return;
	endif
	endif
	%
	%
	error( "Cholesky factorization failed even for a mu just beyond muCrit; this should be impossible!" );
	%
	%
endfunction
