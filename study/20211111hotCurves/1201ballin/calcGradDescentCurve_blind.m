function matX = calcGradDescentCurve_blind( funchOmega, funchGrad, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve_blind";
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
	stepSize = mygetfield( prm, "stepSize", 0.1 );
	iterLimit = mygetfield( prm, "iterLimit", 1000 );
	omegaTol = 1e-8;
	gradNormTol = 1e-8;
	%
	vecX = vecX0;
	omega = omega0;
	vecG = vecG0;
	numPts = 1;
	matX(:,numPts) = vecX;
	iterCount = 0;
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
		deltaT = stepSize/norm(vecG);
		switch (stepOrder)
		case 1
			vecDelta = -deltaT*vecG;
		case -1
			if (1==iterCount)
				vecDelta = -deltaT*vecG;
			else
				% Adams-Bashforth (multistep method)?
				% But, doesn't work very well (blindly).
				vecDelta = -0.5*deltaT*(3.0*vecG-vecGPrev);
			end
		otherwise
			error( "Invalid value of stepOrder." );
		end
		omegaPrev = omega;
		vecGPrev = vecG;
		%
		vecX = vecX + vecDelta;
		omega = funchOmega(vecX);
		if ( omega > 10.0*omega0 )
			return;
		end
		%if ( omega > omegaPrev )
		%	return;
		%end
		numPts++;
		matX(:,numPts) = vecX;
		vecG = funchGrad( vecX );
	end
end

