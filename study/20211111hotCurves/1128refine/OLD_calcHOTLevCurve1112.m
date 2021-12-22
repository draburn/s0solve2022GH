function matY = OLD_calcHOTLevCurve1112( modelFuncPrm, vecX0, prm=[] )
	thisFile = "OLD_calcHOTLevCurve1112";
	msg( thisFile, __LINE__, "This function deserves revision." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	%
	numSteps = mygetfield( prm, "numSteps", 100 );
	assert( isposintscalar(numSteps) );
	finalStepEps = mygetfield( prm, "finalStepEps", 0.5 );
	assert( isrealscalar(finalStepEps) );
	assert( 0.0 <= finalStepEps );
	%
	matI = eye(sizeX,sizeX);
	n = 1;
	vecX = vecX0;
	while (1)
		matY(:,n) = vecX;
		n++;
		if ( n > numSteps )
			break;
		end
		s = n/(numSteps+finalStepEps);
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		vecG = s*vecG + (1.0-s)*(vecX-vecX0);
		matH = s*matH + (1.0-s)*matI;
		[ matR, cholFlag ] = chol( matH );
		if ( 0 ~= cholFlag )
			msg( thisFile, __LINE__, "Hessian is non-positive-definite." );
			%error( "Hessian is non-positive-definite." );
			return;
		else
			vecDelta = -( matR \ (matR'\vecG) );
		end
		vecX += vecDelta;
	end
	%
return;
end
