% Function...

function [ matHRegu, datOut ] = calcHRegu( matH, prm=[] )
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	assert( isrealscalar(cholSafeTol) );
	assert( 0.0 < cholSafeTol  );
	datOut = [];
	%
	matHRegu = matH;
	[ matR, cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
	assert( isrealscalar(epsRelRegu) );
	assert( 0.0 < epsRelRegu  );
	sizeX = size(matH,1);
	%
	hNorm = sqrt( sum(sumsq(matH))/sizeX );
	assert( 0.0 ~= hNorm );
	matHRegu = matH + ( epsRelRegu * hNorm * eye(sizeX,sizeX) );
	[ matR, cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	% There may be better ways to find muCrit.
	[ matPsi, matLambda ] = eig( matH );
	muCrit = -min(diag(matLambda));
	mu = muCrit + epsRelRegu * ( muCrit + hNorm );
	matHRegu = matH + ( mu * eye(sizeX,sizeX) );
	[ matR cholFlag ] = chol( matHRegu );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		return;
	endif
	endif
	%
	error( "Regularization of Hessian failed." );
	%
endfunction


%!test
%!	matHRegu = calcHRegu( [1.0,0.0;0.0,0.0] )
%!	inv(matHRegu)
