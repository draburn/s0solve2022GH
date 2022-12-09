function vecX = mycholdiv( matH, vecB, cholTol=sqrt(eps), prm=[] )
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	%
	assert( isnumeric(matH) );
	assert( isreal(matH) );
	assert( issymmetric(matH,sqrt(eps)) );
	assert( sum(sum(abs(matH))) > 0.0 );
	%
	epsHCoeff = mygetfield( prm, "epsHCoeff", sqrt(eps) );
	matH += (epsHCoeff * (max(abs(diag(matR)))^2)) * eye(size(matH));
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	%
	if ( mygetfield( prm, "useEig", false ) )
		lambda = eig( matH );
		minLambda = min(lambda);
		maxAbsLambda = max(abs(lambda));
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
