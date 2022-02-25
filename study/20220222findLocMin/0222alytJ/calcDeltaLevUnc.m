% Function...
%  [ vecDelta, datOut ] = calcDeltaLevUnc( vecG, matH, prm=[] )
% Calculates vecDelta corresponding to the unconstrained minimization of the omega model.
% If Hessian is not positive definite, the result may be somewhat arbitrarily large.

function [ vecDelta, datOut ] = calcDeltaLevUnc( vecG, matH, prm=[] )
	%
	%
	%
	sizeX = size(vecG,1);
	debugMode = mygetfield( prm, "debugMode", false );
	if ( debugMode )
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isrealarray(vecG,[sizeX,1]) );
		assert( isrealarray(matH,[sizeX,sizeX]) );
		assert( issymmetric(matH) );
	endif
	hNorm = sqrt(sum(sumsq(matH)));
	matI = eye(sizeX,sizeX);
	muReguCoeff = 1.0e-5;
	if ( ~isempty(prm) )
		muReguCoeff = mygetfield( prm, "muReguCoeff", muReguCoeff );
	endif
	if ( debugMode )
		assert( isrealscalar(muReguCoeff) );
		assert( 0.0 < muReguCoeff );
	endif
	%
	% Default return values.
	vecDelta = zeros(sizeX,1);
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	% Not bothering to check if vecG = 0; also, user could possibly want mu and matR.
	if ( 0.0 == hNorm )
		error( "Input Hessian is zero." );
	end
	%
	%
	%
	%%%safeRelTol = sqrt(eps); % Pre "case 99041968".
	safeRelTol = eps^0.35;
	mu = 0.0;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > safeRelTol*max(abs(diag(matR))) )
		vecDelta = -( matR \ (matR'\vecG) );
		if ( nargout >= 2 )
			datOut.mu = mu;
			datOut.matR = matR;
		endif
		return;
	endif
	endif
	msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with mu = 0.0 failed." );
	%
	%
	mu = muReguCoeff*hNorm;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 == cholFlag )
	if ( min(abs(diag(matR))) > safeRelTol*max(abs(diag(matR))) )
		vecDelta = -( matR \ (matR'\vecG) );
		if ( nargout >= 2 )
			datOut.mu = mu;
			datOut.matR = matR;
		endif
		return;
	endif
	endif
	msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with small positive mu failed." );
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
	if ( min(abs(diag(matR))) > safeRelTol*max(abs(diag(matR))) )
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


%!test
%!	msg( __FILE__, __LINE__, "~~~ Positive Definite Test ~~~ " );
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = []
%!	[ vecDelta, datOut ] = calcDeltaLevUnc( vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Positive Semi-Definite Test ~~~ " );
%!	vecG = [ 1.0; 0.0 ]
%!	matH = [ 1.0, 0.0; 0.0, 0.0 ]
%!	prm = []
%!	[ vecDelta, datOut ] = calcDeltaLevUnc( vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Negative Definite / Indefinite Test ~~~ " );
%!	vecG = [ 1.0; 0.0 ]
%!	matH = -eye(2,2)
%!	prm = []
%!	[ vecDelta, datOut ] = calcDeltaLevUnc( vecG, matH, prm )
