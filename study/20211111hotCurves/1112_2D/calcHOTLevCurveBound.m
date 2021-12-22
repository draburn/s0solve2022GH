function matY = calcHOTLevCurveBound( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTLevCurveBound";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	%
	maxNumIter = mygetfield( prm, "maxNumIter", 10000 );
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
	ptIndex = 1;
	totalStepSize = stepSize;
	matY(:,ptIndex) = vecX;
	while (1)
		n++;
		if ( n > maxNumIter )
			break;
		end
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		gNorm = norm(vecG);
		if ( 0.0 == gNorm )
			return;
		end
		vecX -= stepSize*vecG/norm(vecG);
		vecD = vecX-vecX0;
		dNorm = norm(vecD);
		if ( dNorm > totalStepSize + stepSize )
			vecD = (totalStepSize-0.5*stepSize)*vecD/dNorm;
			vecX = vecX0 + vecD;
			totalStepSize += 0.01*stepSize;
		end
		%
		ptIndex++;
		matY(:,ptIndex) = vecX;
	end
	%
return;
end
