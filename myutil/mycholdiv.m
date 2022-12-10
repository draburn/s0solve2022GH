function [ vecX, matR, vecLambda ] = mycholdiv( matH, vecB, forceSolution=false, prm=[] )
	vecLambda = [];
	cholTol = mygetfield( prm, "cholTol", sqrt(eps) );
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	endif
	%
	warningLevel = mygetfield( prm, "warningLevel", 10 );
	assert( isnumeric(matH) );
	assert( isreal(matH) );
	assert( issymmetric(matH,sqrt(eps)) );
	assert( sum(sum(abs(matH))) > 0.0 );
	if ( min(diag(matH)) <= 0.0 )
		if ( 1 <= warningLevel )
			msg( __FILE__, __LINE__, "WARNING: matrix has a non-positive diagonal element." );
		endif
	endif
	if ( 20 <= warningLevel )
		msg( __FILE__, __LINE__, "WARNING: chol() on original matrix was unsuccessful." );
	endif
	%
	hScale = max(sum(matH.^2,1));
	epsHCoeff = mygetfield( prm, "epsHCoeff", sqrt(eps) );
	matH += (epsHCoeff * hScale) * eye(size(matH));
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
		vecX = matR \ ( matR' \ vecB );
		return;
	endif
	if ( 10 <= warningLevel )
		msg( __FILE__, __LINE__, "WARNING: chol() on weakly perturbed matrix was unsuccessful." );
		if ( 0 )
			echo__minsizeH = min(size(matH))
			echo__mindiagH = min(diag(matH))
			echo__cholFlag = cholFlag
			echo__minsizeR = min(size(matR))
			echo__mindiagR = min(diag(matR))
			echo__maxabsdiagR = max(abs(diag(matR)))
		endif
	endif
	%
	if (~forceSolution)
		% Matrix pretty clearly has a negative eigenvalue.
		% Require caller to deal with it.
		vecX = [];
		return;
	endif
	%
	switch ( tolower(mygetfield( prm, "strongPerturbMethod", "basic" )) )
	case { "eig" }
		vecLambda = eig( (matH'+matH)/2.0 );
		minLambda = min(vecLambda);
		maxAbsLambda = max(abs(vecLambda));
		assert( minLambda < sqrt(eps)*maxAbsLambda ); % Really, minLambda should be negative.
		hStrongPerturb = abs(minLambda) + epsHCoeff*maxAbsLambda
	case { "basic" }
		hStrongPerturb = hScale * (1.0+epsHCoeff);
	otherwise
		error( "Unsupported value of 'strongPerturbMethod'." );
	endswitch
	matH += hStrongPerturb * eye(size(matH));
	matR = chol(matH);
	vecX = matR \ ( matR' \ vecB );
return;
endfunction
