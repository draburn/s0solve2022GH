function matY = calcHOTGradCurveRK4( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTGradCurveRK4";
	msg( thisFile, __LINE__, "This function deserves revision." );
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
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD1 = -vecG/gNorm;
		%
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX+(vecD1*stepSize/2.0), modelFuncPrm );
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD2 = -vecG/gNorm;
		%
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX+(vecD2*stepSize/2.0), modelFuncPrm );
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		%
		vecD3 = -vecG/gNorm;
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX+(vecD3*stepSize), modelFuncPrm );
		gNorm = norm(vecG);
		if ( gNorm < gNormTol )
			break;
		end
		vecD4 = -vecG/gNorm;
		%
		vecDelta = stepSize*( (vecD1) + (2.0*vecD2) + (2.0*vecD3) + (vecD4) )/6.0;
		if ( norm(vecDelta) < deltaNormRelTol*stepSize*(norm(vecD1)+2.0*norm(vecD2)+2.0*norm(vecD3)+norm(vecD4))/6.0 )
			break;
		end
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
