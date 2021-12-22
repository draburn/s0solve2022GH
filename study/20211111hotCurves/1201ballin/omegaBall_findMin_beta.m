function vecX = omegaBall_findMin_beta( funchOmega, funchGrad, vecXC, bigR, vecX0, prm=[] )
	thisFile = "omegaBall_findMin_beta";
	%
	sizeX = size(vecXC,1);
	assert( isrealarray(vecXC,[sizeX,1]) );
	assert( isrealscalar(bigR) );
	assert( 0.0 < bigR );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	assert( isrealscalar(iterLimit) );
	%
	stepSizeTol = mygetfield( prm, "stepSizeTol", 1e-4 );
	fallTolRel = mygetfield( prm, "fallTolRel", 1e-4 );
	fallTolAbs = mygetfield( prm, "fallTolAbs", 1e-8 );
	assert( isrealscalar(stepSizeTol) );
	assert( isrealscalar(fallTolRel) );
	assert( isrealscalar(fallTolAbs) );
	assert( 0.0 < stepSizeTol );
	assert( 0.0 < fallTolRel );
	assert( 0.0 < fallTolAbs );
	%
	btLimit = mygetfield( prm, "btLimit", 10 );
	assert( isrealscalar(btLimit) );
	%
	maxStepSize = mygetfield( prm, "maxStepSize", 1.0 );
	btFactor = mygetfield( prm, "btFactor", 0.3 );
	assert( isrealscalar(maxStepSize) );
	assert( isrealscalar(btFactor) );
	assert( stepSizeTol < maxStepSize );
	assert( 0.0 < btFactor );
	assert( btFactor < 1.0 );
	%
	forceOnSurf = mygetfield( prm, "forceOnSurf", false );
	assert( isscalar(forceOnSurf) );
	%
	vecD0 = vecX0 - vecXC;
	% Auto-detect if on surf.
	if ( norm(vecD0) > bigR + 0.1*stepSizeTol )
		onSurf = true;
		vecX = vecXC + vecD0*bigR/norm(vecD0);
	else
		onSurf = false;
		vecX = vecX0;
	end
	%
	% HACK to force onSurf.
	if (forceOnSurf)
	if ( 0.0 == norm(vecD0) )
		% Do nothing in this case.
	else
		onSurf = true;
		vecD = vecX - vecXC;
		vecX = vecXC + bigR*vecD/norm(vecD);
	end
	end
	%
	%
	%
	omega = funchOmega(vecX);
	assert( isrealscalar(omega) );
	assert( omega >= 0.0 );
	if ( 0.0 == omega )
		msg( thisFile, __LINE__, "Starting omega is zero." );
		return;
	end
	%
	iterCount = 0;
	while (1)
		% Update iter count.
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iter limit during BT." );
			return;
		end
		%
		% Calculate gradient.
		vecGrad = funchGrad(vecX);
		assert( isrealvector(vecGrad,[sizeX,1]) );
		%
		% Validate current point and "onSurf" status.
		vecD = vecX - vecXC;
		if ( onSurf && ~forceOnSurf )
			assert( abs(norm(vecD)-bigR) < 10.0*stepSizeTol );
			if ( vecGrad'*vecD > 0.0 )
				% Gradient points outward / descent points inward: fall off surface.
				onSurf = false;
			end
		end
		%	
		% Pick our step direction. (Magnitude doesn't matter.)
		if ( onSurf )
			vecS = vecD*(vecD'*vecGrad)/(norm(vecD)^2) - vecGrad;
		else
			vecS = -vecGrad;
		end
		if ( norm(vecS) < sqrt(eps)*norm(vecGrad) )
			%msg( thisFile, __LINE__, "Hit zero (effective) gradient." );
			return;
		end
		%
		% Look at how omega varies along this direction.
		% vecDelta = t*vecS;
		% omegaModel = omega + t*vecS'*vecG + t^2*h;
		% Note: a more accurate model would be to use something like...
		%  vecDelta = R*sin(theta)*vecSHat + R*(cos(theta)-1)*vecRHat.
		% That may later be useful for optimization, but, not now.
		t = stepSizeTol/norm(vecS);
		if ( onSurf )
			vecXTemp = vecX + t*vecS;
			vecDTemp = vecXTemp - vecXC;
			assert( norm(vecDTemp) ~= 0.0 );
			vecXTemp = vecXC + vecDTemp*bigR/norm(vecDTemp);
			omegaTemp = funchOmega( vecXTemp );
			clear vecXTemp;
			clear vecDTemp;
		else
			omegaTemp = funchOmega( vecX + t*vecS );
		end
		stg = vecS'*vecGrad;
		h = ( omegaTemp - omega - t*stg ) / (t^2);
		t = getQuadGoodPt( h, stg, omega ); % Returns min of omegaModel or appropriate zero.
		if ( t <= 0.0 )
			msg( thisFile, __LINE__, "" );
			msg( thisFile, __LINE__, "" );
			msg( thisFile, __LINE__, "**********     DEBUG ME!     **********" );
			msg( thisFile, __LINE__, "" );
			msg( thisFile, __LINE__, "" );
			%echo__h = h
			%echo__stg = stg
			%echo__omega = omega
			%echo__vecX0 = vecX0
			%echo__vecXC = vecXC
			return;
		end
		%assert( t > 0.0 );
		%
		% Apply maxStepSize, if necessary.
		if ( t*norm(vecS) > maxStepSize )
			t = maxStepSize/norm(vecS);
		end
		%
		% If not on surf, do not go past surf.
		if (~onSurf)
			dts = vecD'*vecS;
			discrim = (dts)^2 + ( bigR^2 - norm(vecD)^2 ) * norm(vecS)^2;
			assert( discrim >= 0.0 );
			tMax = ( -dts + sqrt(discrim) ) / (norm(vecS)^2);
			if ( t > tMax )
				t = tMax;
				tIsMaxValue = true;
			else
				tIsMaxValue = false;
			end
			clear discrim;
			clear tMax;
		end
		%
		% Take the step.
		% Use geometric backtracking if omega fails to decrease.
		btCount = 0;
		while (1)
			vecDelta = t*vecS;
			vecX_trial = vecX + vecDelta;
			vecD_trial = vecX_trial - vecXC;
			if ( onSurf )
				% Stay on surface.
				assert( 0.0 ~= norm(vecD_trial) );
				vecX_trial = vecXC + (vecD_trial*bigR/norm(vecD_trial));
				omega_trial = funchOmega(vecX_trial);
			else
				assert( norm(vecD_trial) < bigR + stepSizeTol );
				omega_trial = funchOmega(vecX_trial);
			end
			%
			% If omega decreases, accept the point.
			if ( omega_trial < omega )
				if ( ~onSurf )
				if ( tIsMaxValue )
					onSurf = true;
					clear tIsMaxValue;
				end
				end
				break;
			end
			%
			if ( norm(vecDelta) < stepSizeTol )
				%msg( thisFile, __LINE__, "Hit BT stepSize limit." );
				return;
			end
			%
			t *= btFactor;
			if ( ~onSurf )
				tIsMaxValue = false;
			end
			%
			btCount++;
			if ( btCount > btLimit )
				msg( thisFile, __LINE__, "Hit BT trial." );
				return;
			end
			%msg( thisFile, __LINE__, "Backtracking..." );
		end
		assert( omega_trial < omega );
		omegaFall = omega - omega_trial;
		%
		%msg( thisFile, __LINE__, "Found next guess." );
		vecX = vecX_trial;
		omega = omega_trial;
		%
		if ( norm(vecDelta) <= stepSizeTol )
			%msg( thisFile, __LINE__, "Reached step size tol." );
			return;
		end
		if (  omegaFall <= omega*fallTolRel  &&  omegaFall <= fallTolAbs )
			%msg( thisFile, __LINE__, "Reached fall tol." );
			return;
		end
	end
	%
return
