function matY = OLD_calcHOTGradCurve2Heun( modelFuncPrm, vecX0, prm=[] )
	thisFile = "OLD_calcHOTGradCurve2Heun";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
%%%	minStepSize = mygetfield( prm, "minStepSize", 0.01 );
%%%	assert( isrealscalar(minStepSize) );
%%%	assert( 0.0 < minStepSize );
%%%	assert( minStepSize <= maxStepSize );
	maxStepSize = mygetfield( prm, "maxStepSize", 0.2 );
	assert( isrealscalar(maxStepSize) );
	assert( 0.0 < maxStepSize );
	stepTol = mygetfield( prm, "stepTol", 0.2 );
	assert( isrealscalar(stepTol) );
	assert( 0.0 < stepTol );
	%
	maxNumIter = mygetfield( prm, "maxNumIter", 1000 );
	assert( isposintscalar(maxNumIter) );
	maxBTIter = mygetfield( prm, "maxBTIter", 10 );
	assert( isposintscalar(maxBTIter) );
	useRK4 = mygetfield( prm, "useRK4", false );
	assert( isscalar(useRK4) );
	assert( isbool(useRK4) );
	%
	matI = eye(sizeX,sizeX);
	vecX = vecX0;
	[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
	assert( 0.0 < norm(vecG) );	
	vecG0 = vecG;
	g0Norm = norm(vecG0);
	omegaPrev = 1.0+omega;
	%
	n = 1;
	while (1)
		gNorm = norm(vecG);
		if ( gNorm <= (eps^0.75)*g0Norm )
			return;
		elseif ( n > maxNumIter )
			return;
		end
		matY(:,n) = vecX;
		%
		n++;
		omegaPrev = omega;
		%
		vecGHat = vecG/gNorm;
		vecHGHat = matH*vecGHat;
		normHGHat = norm(vecHGHat);
		deltaT = maxStepSize;
		if ( normHGHat*maxStepSize > stepTol * gNorm )
			deltaT = stepTol * gNorm / normHGHat;
		end
		if (useRK4)
			assert(0);
			vecK1 = -vecGHat;
			[ foo1, foo2, foo3, vecG_temp, foo4 ] = testFunc_evalDeriv( vecX+0.5*deltaT*vecK1, modelFuncPrm );
			vecK2 = -vecG_temp/norm(vecG_temp);
			[ foo1, foo2, foo3, vecG_temp, foo4 ] = testFunc_evalDeriv( vecX+0.5*deltaT*vecK2, modelFuncPrm );
			vecK3 = -vecG_temp/norm(vecG_temp);
			[ foo1, foo2, foo3, vecG_temp, foo4 ] = testFunc_evalDeriv( vecX+deltaT*vecK3, modelFuncPrm );
			vecK4 = -vecG_temp/norm(vecG_temp);
			vecDelta = deltaT*(vecK1 + 2.0*vecK2 + 2.0*vecK3 + vecK4)/6.0;
		else
		%%%vecDelta = ( matI + 0.5*deltaT*(matH-vecGHat*(vecGHat'*matH))/gNorm ) \ (-deltaT*vecGHat);
		%deltaT *= 2.0;
		vecDeltaFWD = ( matI + 0.5*deltaT*(matH-vecGHat*(vecGHat'*matH))/gNorm ) \ (-deltaT*vecGHat);
		%
		vecX_temp = vecX+vecDeltaFWD;
		[ omega_temp, vecF_temp, matJ_temp, vecG_temp, matH_temp ] = testFunc_evalDeriv( vecX_temp, modelFuncPrm );
		gNorm_temp = norm(vecG_temp);
		assert( 0.0 < gNorm_temp );
		vecGHat_temp = vecG_temp/gNorm_temp;
		%%%
		%%%vecDeltaBWD =  ( matI + 0.5*deltaT*(matH_temp-vecGHat_temp*(vecGHat_temp'*matH_temp))/gNorm_temp ) \ (-deltaT*vecGHat_temp);
		vecDeltaBWD =  ( matI - 0.5*deltaT*(matH_temp-vecGHat_temp*(vecGHat_temp'*matH_temp))/gNorm_temp ) \ (-deltaT*vecGHat_temp);
		%%%
		vecDelta = (vecDeltaFWD + vecDeltaBWD)/2.0;
		%deltaT /= 2.0;
		end
		%
		% But, if Newton step is good and close, use that instead.
		[ matR, cholFlag ] = chol( matH );
		if ( 0==cholFlag )
			vecDeltaN = -( matR \ (matR'\vecG) );
			if ( norm(vecDeltaN)<=maxStepSize )
				vecDelta = vecDeltaN;
			end
		end
		%
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX+vecDelta, modelFuncPrm );
		btIter = 0;
		while ( omega >= omegaPrev )
			btIter++;
			if ( btIter > maxBTIter )
				return;
			end
			vecDelta /= 2.0;
			[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX+vecDelta, modelFuncPrm );
		end
		%
		vecX += vecDelta;
	end
	return;
end
