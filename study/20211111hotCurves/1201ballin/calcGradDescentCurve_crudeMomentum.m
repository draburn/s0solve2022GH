function matX = calcGradDescentCurve_crudeMomentum( funchOmega, funchGrad, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve_crudeMomentum";
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
	vecV0 = zeros(sizeX,1);
	iterLimit = 1000;
	omegaTol = 1e-8;
	gradNormTol = 1e-8;
	%btLimit = 10;
	%btCoeff = 0.5;
	massCoeff = 10.0;
	dampingCoeff = 0.0;
	deltaT = 0.1;
	%
	%
	matX(:,1) = vecX0;
	iterCount = 0;
	vecX = vecX0;
	omega = omega0;
	vecG = vecG0;
	vecV = vecV0;
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
		vecA = -(vecG/massCoeff);
		vecX_trial = vecX + vecV*deltaT + 0.5*vecA*(deltaT^2);
		vecV = (1.0-dampingCoeff)*vecV + vecA*deltaT;
		omega_trial = funchOmega(vecX_trial);
		if ( omega_trial > omega )
			vecX_trial = vecX;
			vecV *= 0.0;
		end
		%
		%btCount = 0;
		%while (1)
		%	vecX_trial = vecX + deltaT*vecV;
		%	omega_trial = funchOmega(vecX_trial);
		%	if ( omega_trial < omega )
		%		break;
		%	end
		%	%
		%	btCount++;
		%	if ( btCount >= btLimit )
		%		msg( thisFile, __LINE__, "Reached btLimit." );
		%		return;
		%	end
		%	%
		%	vecV *= btCoeff;
		%end
		%assert( omega_trial < omega );
		vecX = vecX_trial;
		omega = omega_trial;
		vecG = funchGrad(vecX);
		matX(:,iterCount+1) = vecX;
	end
	%
return;
end
