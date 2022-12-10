function [ vecX, matR, vecLambda ] = mycholdiv( matH, vecB, cholTol=sqrt(eps), prm=[] )
	vecLambda = [];
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	endif
	%
	assert( isnumeric(matH) );
	assert( isreal(matH) );
	assert( issymmetric(matH,sqrt(eps)) );
	assert( sum(sum(abs(matH))) > 0.0 );
	%
	epsHCoeff = mygetfield( prm, "epsHCoeff", sqrt(eps) );
	matH += (epsHCoeff * max(abs(diag(matH)))) * eye(size(matH));
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	endif
	%
	if ( mygetfield( prm, "useEig", false ) )
		vecLambda = eig( (matH'+matH)/2.0 );
		minLambda = min(vecLambda);
		maxAbsLambda = max(abs(vecLambda));
		assert( minLambda < sqrt(eps)*maxAbsLambda ); % Really, minLambda should be negative.
		matH += ( abs(minLambda) + epsHCoeff*maxAbsLambda) * eye(size(matH));
	else
		matH += sum(sum(abs(matH))) * eye(size(matH));
	endif
	%
	matR = chol(matH);
	vecX = matR \ ( matR' \ vecB );
return;
endfunction
