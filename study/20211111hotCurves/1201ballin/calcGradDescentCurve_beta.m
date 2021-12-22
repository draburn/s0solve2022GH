function matX = calcGradDescentCurve_beta( funchOmega, funchGrad, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve_beta";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	stepSize = mygetfield( prm, "stepSize", 0.1 );
	assert( isrealscalar(stepSize) );
	assert( 0.0 < stepSize );
	%
	iterLimit = mygetfield( prm, "iterLimit", 1000 );
	assert( isrealscalar(iterLimit) );
	%
	doLev = mygetfield( prm, "doLev", false );
	assert( isscalar(doLev) );
	%
	minStepSize = mygetfield( prm, "minStepSize", stepSize/10.0 );
	assert( isrealscalar(minStepSize) );
	assert( 0.0 < minStepSize );
	%
	%
	matX(:,1) = vecX0;
	omega0 = funchOmega(vecX0);
	%
	vecG = funchGrad(vecX0);
	assert( 0.0 ~= norm(vecG) );
	vecX = vecX0 - 0.01*stepSize*vecG/norm(vecG);
	omega = omega0;
	iterCount = 0;
	while (1)
		% Update iter count.
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iter limit." );
			return;
		end
		%
		obfmbPrm.forceOnSurf = true;
		if (doLev)
			vecX_trial = omegaBall_findMin_beta( funchOmega, funchGrad, vecX0, iterCount*stepSize, vecX, obfmbPrm );
		else
			vecX_trial = omegaBall_findMin_beta( funchOmega, funchGrad, vecX, stepSize, vecX, obfmbPrm );
		end
		omega_trial = funchOmega(vecX_trial);
		if ( omega_trial >= omega )
			msg( thisFile, __LINE__, "Failed to decrease omega." );
			return;
		end
		matX(:,iterCount+1) = vecX_trial;
		if ( norm(vecX_trial-vecX) < minStepSize )
			msg( thisFile, __LINE__, "Reached minStepSize." );
			return;
		end
		vecX = vecX_trial;
		omega = omega_trial;
		%
		% Make selection of next vecX smarter:
		% instead of just prevX, extrapolate based on curve.
	end
return;
end
