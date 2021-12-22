function matY = calcHOTLevCurveRK4( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTLevCurveRK4";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	%
	maxNumIter = mygetfield( prm, "maxNumIter", 1000 );
	assert( isposintscalar(maxNumIter) );
	stepSize = mygetfield( prm, "stepSize", 0.1 );
	assert( isrealscalar(stepSize) );
	gNormTol = mygetfield( prm, "gNormTol", eps );
	assert( isrealscalar(gNormTol) );
	deltaNormRelTol = mygetfield( prm, "deltaNormRelTol", 0.01 );
	assert( isrealscalar(deltaNormRelTol) );
	%
	n = 1;
	vecX = vecX0;
	while (1)
		matY(:,n) = vecX;
		n++;
		if ( n > maxNumIter )
			break;
		end
		s = 2.0*n/(0.1+maxNumIter);
		s = min([ s, 1.0 ]);
		%
		vecXTemp = vecX;
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecXTemp, modelFuncPrm );
		vecG = s*vecG + (1.0-s)*(vecXTemp-vecX0);
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD1 = -vecG/gNorm;
		%
		vecXTemp = vecX+(vecD1*stepSize/2.0);
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecXTemp, modelFuncPrm );
		vecG = s*vecG + (1.0-s)*(vecXTemp-vecX0);
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD2 = -vecG/gNorm;
		%
		vecXTemp = vecX+(vecD2*stepSize/2.0);
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecXTemp, modelFuncPrm );
		vecG = s*vecG + (1.0-s)*(vecXTemp-vecX0);
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD3 = -vecG/gNorm;
		%
		vecXTemp = vecX+(vecD3*stepSize);
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecXTemp, modelFuncPrm );
		vecG = s*vecG + (1.0-s)*(vecXTemp-vecX0);
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD4 = -vecG/gNorm;
		%
		vecXTemp = vecX;
		vecDelta = stepSize*( (vecD1) + (2.0*vecD2) + (2.0*vecD3) + (vecD4) )/6.0;
		vecX += vecDelta;
	end
	%
	[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
		vecDelta = -( matR \ (matR'\vecG) );
		vecX += vecDelta;
		matY(:,n) = vecX;
	end
	%
return;
end
