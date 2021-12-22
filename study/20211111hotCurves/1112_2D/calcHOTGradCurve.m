function matY = calcHOTGradCurve( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTGradCurve";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	minStepSize = mygetfield( prm, "minStepSize", 0.002 );
	maxStepSize = mygetfield( prm, "maxStepSize", 0.02 );
	assert( isrealscalar(minStepSize) );
	assert( isrealscalar(maxStepSize) );
	assert( 0.0 < minStepSize );
	assert( minStepSize <= maxStepSize );
	%
	maxNumIter = mygetfield( prm, "maxNumIter", 10000 );
	assert( isposintscalar(maxNumIter) );
	%
	vecX = vecX0;
	[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
	assert( 0.0 < norm(vecG) );	
	vecG0 = vecG;
	omegaPrev = 1.0+omega;
	%
	n = 1;
	while (1)
		if ( norm(vecG) <= eps*norm(vecG0) )
			break;
		elseif ( omega >= omegaPrev )
			break;
		elseif ( n > maxNumIter )
			break;
		end
		matY(:,n) = vecX;
		%
		n++;
		omegaPrev = omega;
		vecDelta = -minStepSize*vecG/norm(vecG);
		[ matR, cholFlag ] = chol( matH );
		if ( 0==cholFlag )
			vecDeltaN = -( matR \ (matR'\vecG) );
			if ( norm(vecDeltaN)<=maxStepSize )
				vecDelta = vecDeltaN;
			end
		else
			%msg( thisFile, __LINE__, "Hessian is non-positive definite! HACK: FORCING STOP." );
			%return;
		end
		vecX += vecDelta;
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
	end
	return;
end
