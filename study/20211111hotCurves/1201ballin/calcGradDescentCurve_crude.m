function matX = calcGradDescentCurve_crude( funchOmega, funchGrad, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve_crude";
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
	stepOrder = mygetfield( prm, "stepOrder", 1 );
	%
	iterLimit = 1000;
	stepSize0 = 0.1;
	omegaTol = 1e-8;
	gradNormTol = 1e-8;
	decelCoeff = 0.5;
	accelCoeff = 1.1;
	btLimit = 10;
	minStepSize = 1e-4;
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
		if (stepSize > stepSize0)
			stepSize = stepSize0;
		end
		deltaT = stepSize/norm(vecG);
		switch (stepOrder)
		case 1
			vecDelta = -deltaT*vecG;
		case 2
			vecDelta = -deltaT*0.5*( vecG + funchGrad(vecX-deltaT*vecG) );
		case 4
			vecG0 = vecG;
			vecG1 = funchGrad(vecX-deltaT*vecG0/2.0);
			vecG2 = funchGrad(vecX-deltaT*vecG1/2.0);
			vecG3 = funchGrad(vecX-deltaT*vecG2);
			vecDelta = -deltaT*( vecG0 + 2.0*vecG1 + 2.0*vecG2 + vecG3 ) / 6.0;
		otherwise
			error( "Invalid value of stepOrder." );
		end
		btCount = 0;
		while (1)
			vecX_trial = vecX + vecDelta;
			omega_trial = funchOmega(vecX_trial);
			if ( omega_trial < omega )
				stepSize *= accelCoeff;
				break;
			end
			btCount++;
			if ( btCount >= btLimit )
				msg( thisFile, __LINE__, "Reached btLimit." );
				return;
			end
			stepSize *= decelCoeff;
			deltaT = stepSize/norm(vecG);
			switch (stepOrder)
			case 1
				vecDelta = -deltaT*vecG;
			case 2
				vecDelta = -deltaT*0.5*( vecG + funchGrad(vecX-deltaT*vecG) );
			case 4
				vecG0 = vecG;
				vecG1 = funchGrad(vecX-deltaT*vecG0/2.0);
				vecG2 = funchGrad(vecX-deltaT*vecG1/2.0);
				vecG3 = funchGrad(vecX-deltaT*vecG2);
				vecDelta = -deltaT*( vecG0 + 2.0*vecG1 + 2.0*vecG2 + vecG3 ) / 6.0;
			otherwise
				error( "Invalid value of stepOrder." );
			end
			if ( norm(vecDelta) < minStepSize )
				msg( thisFile, __LINE__, "Reached minStepSize." );
				return;
			end
		end
		assert( omega_trial < omega );
		vecX = vecX_trial;
		omega = omega_trial;
		vecG = funchGrad(vecX);
		matX(:,iterCount+1) = vecX;
	end
	%
return;
end
