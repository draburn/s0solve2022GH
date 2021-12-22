function vecX = omegaBall_findMin_alpha( funchOmega, funchG, vecXC, bigR, bigC, vecX0, prm=[] )
	thisFile = "omegaBall_findMin_alpha";
	%
	sizeX = size(vecXC,1);
	assert( isrealarray(vecXC,[sizeX,1]) );
	assert( isrealscalar(bigR) );
	assert( 0.0 < bigR );
	assert( isrealscalar(bigC) );
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
	assert( 0.0 < maxStepSize );
	assert( 0.0 < btFactor );
	assert( btFactor < 1.0 );
	%
	funchOmegaBall = @(x) omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, x );
	funchGBall = @(x) omegaBall_evalGrad( funchG, vecXC, bigR, bigC, x );
	%
	vecX = vecX0;
	iterCount = 0;
	while (1)
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iter limit." );
			return;
		end
		%
		omega = funchOmegaBall(vecX);
		assert( isrealscalar(omega) );
		assert( omega >= 0.0 );
		if ( 0.0 == omega )
			msg( thisFile, __LINE__, "Hit zero cost." );
			return;
		end
		%
		vecG = funchGBall(vecX);
		assert( isrealvector(vecG,[sizeX,1]) );
		gtg = vecG'*vecG;
		if ( 0.0 == gtg )
			msg( thisFile, __LINE__, "Hit zero gradient." );
			return;
		end
		%
		s = stepSizeTol/sqrt(gtg);
		omegaP = funchOmegaBall( vecX - s*vecG );
		%
		% omegaModel = omega - s * (vecG'*vecG) + c2 * s^2...
		c2 = ( omegaP - omega + s * gtg ) / (s^2);
		s = getQuadGoodPt( c2, -gtg, omega ); % Returns min of omegaModel or appropriate zero.
		assert( s > 0.0 );
		%
		vecDelta = -s*vecG;
		if ( norm(vecDelta) > maxStepSize )
			s = maxStepSize/sqrt(gtg);
			vecDelta = -s*vecG;
			assert( norm(vecDelta) < maxStepSize*1.0001 );
		end
		%
		btCount = 0;
		while (1)
			%
			vecDelta = -s*vecG;
			vecX_trial = vecX + vecDelta;
			omega_trial = funchOmegaBall(vecX_trial);
			if ( omega_trial < omega )
				break;
			end
			%
			if ( norm(vecDelta) < stepSizeTol )
				msg( thisFile, __LINE__, "Hit BT stepSize limit." );
				return;
			end
			%
			s *= btFactor;
			btCount++;
			if ( btCount > btLimit )
				msg( thisFile, __LINE__, "Hit BT trial." );
				return;
			end
			msg( thisFile, __LINE__, "Backtracking..." );
		end
		assert( omega_trial < omega );
		omegaFall = omega - omega_trial;
		%
		msg( thisFile, __LINE__, "Found next guess." );
		vecX = vecX_trial;
		omega = omega_trial;
		%
		if ( norm(vecDelta) <= stepSizeTol )
			msg( thisFile, __LINE__, "Reached step size tol." );
			return;
		end
		if (  omegaFall <= omega*fallTolRel  &&  omegaFall <= fallTolAbs )
			msg( thisFile, __LINE__, "Reached fall tol." );
			return;
		end
	end
	%
return
