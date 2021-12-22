function matX = calcGradDescentCurve_eulerMin( funchOmega, funchGrad, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve_eulerMin";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	omega0 = funchOmega( vecX0 );
	assert( isrealscalar(omega0) );
	assert( 0.0 <= omega0 );
	vecG0 = funchGrad( vecX0 );
	assert( isrealarray(vecG0,[sizeX,1]) );
	%
	iterLimit = 1000;
	stepSize0 = 0.1;
	omegaTol = 1e-8;
	gradNormTol = 1e-8;
	btLimit = 3;
	minStepSize = 1e-4;
	accelCoeff = 1.1;
	%
	%
	matX(:,1) = vecX0;
	iterCount = 0;
	stepSize = stepSize0;
	vecX = vecX0;
	omega = omega0;
	vecG = vecG0;
	while (1)
		if ( omega <= omegaTol )
			msg( thisFile, __LINE__, "Reached omegaTol." );
			return;
		end
		if ( norm(vecG) <= gradNormTol )
			msg( thisFile, __LINE__, "Reached gradNormTol." );
			return;
		end
		if ( iterCount >= iterLimit )
			msg( thisFile, __LINE__, "Reached iterLimit." );
			return;
		end
		iterCount++;
		%
		gtg = vecG'*vecG;
		deltaT = stepSize0/sqrt(gtg);
		%
		vecDelta = -deltaT*vecG;
		vecX_trial = vecX + vecDelta;
		omega_trial = funchOmega(vecX_trial);
		vecG_trial = funchGrad(vecX_trial);
		%
		% omegaModel = omega - norm(vecG)^2*deltaT + c2*(deltaT^2);
		c2 = ( omega_trial - omega + gtg*deltaT ) / deltaT^2;
		if ( c2 > 0.0 )
			deltaTOfMin = gtg/(2.0*c2);
			if ( deltaTOfMin < deltaT )
				vecDelta = -deltaTOfMin*vecG;
				vecX_trial = vecX + vecDelta;
				omega_trial = funchOmega(vecX_trial);
				vecG_trial = funchGrad(vecX_trial);
			end
		end
		%
		assert( omega_trial < omega );
		vecX = vecX_trial;
		omega = omega_trial;
		vecG = vecG_trial;
		matX(:,iterCount+1) = vecX;
	end
	%
return;
end
