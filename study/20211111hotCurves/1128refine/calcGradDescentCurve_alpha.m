function matX = calcGradDescentCurve_alpha( funchGOmega, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve_alpha";
	msg( thisFile, __LINE__, "WORK-IN-PROGRESS." );
	msg( thisFile, __LINE__, " - Use stepAccuracyTarget to update deltaT?" );
	msg( thisFile, __LINE__, " - Consider larger step along delta!" );
	msg( thisFile, __LINE__, " - Measure step accuracy using RK5/4(?) instead of RK4/Heunesque?" );
	msg( thisFile, __LINE__, " - Use mygetfield, validate prm." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	%
	% RK4 + BT...
	fevalCountMax = 100000;
	numStepsMax = 10000;
	gradNormTol = 0.0;
	stepSizeTarget = 0.1;
	stepSizeTol = 1e-4;
	omegaDecreaseTol0 = 1e-6;
	omegaDecreaseTol1 = 1e-3;
	decelCoeff = 0.5;
	accelCoeff = 1.1;
	stepAccuracyTolLo = 0.05;
	stepAccuracyTolHi = 0.2;
	%
	matX(:,1) = vecX0;
	[ vecG0, omega0 ] = funchGOmega( vecX0 );
	fevalCount = 1;
	if ( gradNormTol >= norm(vecG0) )
		msg( thisFile, __LINE__, sprintf( ...
		  "Initial gradient is too small ( %0.3e >= %0.3e ).", gradNormTol, norm(vecG0) ) );
		return;
	end
	deltaT0 = stepSizeTarget/norm(vecG0);
	%
	vecX = vecX0;
	vecG = vecG0;
	omega = omega0;
	deltaT = deltaT0;
	%
	numSteps = 0;
	while (1) % main loop
		if ( numSteps >= numStepsMax )
			break;
		end
		if ( gradNormTol >= norm(vecG) )
			break;
		end
		if ( deltaT > stepSizeTarget/norm(vecG)  )
			deltaT = stepSizeTarget/norm(vecG);
		end
		%
		while (1) % BT search.
			if ( fevalCount >= fevalCountMax )
				msg( thisFile, __LINE__, "Reached fevalCountMax." );
				return;
			end
			%
			vecGD1 = -vecG;
			vecGD2 = -funchGOmega( vecX + 0.5*deltaT*vecGD1 );
			vecGD3 = -funchGOmega( vecX + 0.5*deltaT*vecGD2 );
			vecGD4 = -funchGOmega( vecX + deltaT*vecGD3 );
			fevalCount += 3;
			vecDelta = deltaT*( vecGD1 + (2.0*vecGD2) + (2.0*vecGD3) + vecGD4 ) / 6.0;
			vecDeltaHeunEsque = deltaT*( vecGD1 + vecGD4 ) / 2.0;
			if ( norm(vecDelta-vecDeltaHeunEsque) > stepAccuracyTolHi*(norm(vecDelta)+norm(vecDeltaHeunEsque)) )
				deltaT *= decelCoeff;
				continue;
			end
			%
			vecX_trial = vecX + vecDelta;
			[ vecG_trial, omega_trial ] = funchGOmega( vecX_trial );
			fevalCount++;
			if ( omega_trial > omega )
				deltaT *= decelCoeff;
				continue;
			end
			%
			if ( norm(vecDelta-vecDeltaHeunEsque) < stepAccuracyTolLo*(norm(vecDelta)+norm(vecDeltaHeunEsque)) )
				deltaT *= accelCoeff;
			end
			break;
		end
		if ( omega_trial > omega )
			error( "Omega increased; this should be impossible." );
		end
		omegaDecrease = omega-omega_trial;
		%
		vecX = vecX_trial;
		vecG = vecG_trial;
		omega = omega_trial;
		numSteps++;
		matX(:,numSteps+1) = vecX;
		if ( stepSizeTol >= norm(vecDelta) )
			msg( thisFile, __LINE__, "Reached step size tolerance." );
			break;
		end
		if ( omegaDecrease < omegaDecreaseTol1*omega && omegaDecrease < omegaDecreaseTol0*omega0 )
			msg( thisFile, __LINE__, "Reached omega decrease tolerances." );
			break;
		end
	end
	return;
	%
return;
end
