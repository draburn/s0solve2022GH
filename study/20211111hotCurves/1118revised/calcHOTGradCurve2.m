function matX = calcHOTGradCurve( funchOmegaG, vecX0, prm=[] )
	thisFile = "calcHOTGradCurve";
	msg( thisFile, __LINE__, "WORK-IN-PROGRESS." );
	fevalCount = 0;
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	[ omega0, vecG0 ] = funchOmegaG( vecX0 ); fevalCount++;
	assert( isrealscalar(omega0) );
	assert( 0.0 <= omega0 );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( 0.0 < norm(vecG0) );
	%
	deltaXNormMin = mygetfield( prm, "deltaXNormMin", 1e-4 );
	deltaXNormMax = mygetfield( prm, "deltaXNormMax", 1.0 );
	assert( isrealscalar(deltaXNormMin) );
	assert( isrealscalar(deltaXNormMax) );
	assert( 0.0 < deltaXNormMin );
	assert( deltaXNormMin < deltaXNormMax );
	%
	resTolLo = mygetfield( prm, "resTolLo", 0.01 );
	resTolHi = mygetfield( prm, "resTolHi", 1000.0 );
	assert( isrealscalar(resTolLo) );
	assert( isrealscalar(resTolHi) );
	assert( 0.0 < resTolLo );
	assert( resTolLo < resTolHi );
	%
	brakeFactor = mygetfield( prm, "brakeFactor", 0.8 );
	accelFactor = mygetfield( prm, "accelFactor", 1.01 );
	assert( isrealscalar(brakeFactor) );
	assert( isrealscalar(accelFactor) );
	assert( 0.0 < brakeFactor );
	assert( brakeFactor < 1.0 );
	assert( 1.0 < accelFactor );
	assert( brakeFactor*accelFactor < 1.0 ); % Since we'll immedieately re-increase.
	%
	fevalLimit = mygetfield( prm, "fevalLimit", 10000 );
	trialLimit = mygetfield( prm, "trialLimit", 10 );
	assert( isposintscalar(fevalLimit) );
	assert( isposintscalar(trialLimit) );
	%
	ptIndex = 1;
	vecX = vecX0;
	omega = omega0;
	vecG = vecG0;
	matX(:,1) = vecX0;
	%msg( thisFile, __LINE__, sprintf( "Accepting point:  %3d,  %10.3e,  %10.3e,  %10.3e.", ptIndex, 0.0, norm(vecX-vecX0), omega ) );
	%
	%deltaT = 1.9 * sqrt( deltaXNormMax * deltaXNormMin )
	deltaT = 0.1*deltaXNormMax;
	searchingForNextPt = true;
	while (searchingForNextPt)
		ptIndex++;
		%
		%
		%
		searchingForGoodStep = true;
		trialCount = 0;
		while (searchingForGoodStep)
			% Here, we'll decrease deltaT as needed.
			trialCount++;
			if ( trialCount > trialLimit )
				msg( thisFile, __LINE__, sprintf( "Reached trial limit ( %d > %d ).", trialCount, trialLimit ) );
				return;
			end
			%
			% Calculate Heun step with Euler residual.
			if ( 0.0 == norm(vecG) )
				msg( thisFile, __LINE__, sprintf( "Gradient is too small ( %0.3e < %0.3e ) for pt %d.", norm(vecG), norm(vecG0), ptIndex ) );
				return;
			end
			vecDeltaX1 = -deltaT*vecG/norm(vecG);
			vecX1 = vecX + vecDeltaX1;
			[ omega1, vecG1 ] = funchOmegaG( vecX1 ); fevalCount++;
			assert( isrealscalar(omega1) );
			assert( 0.0 <= omega1 );
			assert( isrealarray(vecG1,[sizeX,1]) );
			if ( 0.0 == norm(vecG1) )
				msg( thisFile, __LINE__, sprintf( "Gradient is too small ( %0.3e < %0.3e ) for pt %d.", norm(vecG1), norm(vecG0), ptIndex ) );
				return;
			end
			vecDeltaX2 = -deltaT*vecG1/norm(vecG1);
			vecDeltaX = 0.5 * ( vecDeltaX1 + vecDeltaX2 );
			vecRes = 0.5 * abs( vecDeltaX1 - vecDeltaX2 );
			clear vecDeltaX2;
			clear vecG1;
			clear omega1;
			clear vecX1;
			clear vecDeltaX1;
			%
			%msg( thisFile, __LINE__, sprintf( " %3d, %3d;  %10.3e;  %10.3e, %10.3e.", ...
			%  ptIndex, trialCount, deltaT, norm(vecDeltaX), norm(vecRes) ) );
			%
			%
			%
			if ( norm(vecRes) > resTolHi * norm(vecDeltaX) )
				msg( thisFile, __LINE__, sprintf( "Step was too inaccurate ( %0.3e > %0.3e * %0.3e ) for pt %d.", norm(vecRes), resTolHi, norm(vecDeltaX), ptIndex ) );
				deltaT *= brakeFactor;
				continue;
			end
			%
			if ( norm(vecDeltaX) > deltaXNormMax )
				msg( thisFile, __LINE__, sprintf( "Step size was too large ( %0.3e > %0.3e ) for pt %d.", norm(vecDeltaX), deltaXNormMax, ptIndex ) );
				deltaT *= brakeFactor;
				continue;
			end
			%
			vecX_trial = vecX + vecDeltaX;
			[ omega_trial, vecG_trial ] = funchOmegaG( vecX_trial ); fevalCount++;
			assert( isrealscalar(omega) );
			assert( 0.0 <= omega );
			assert( isrealarray(vecG,[sizeX,1]) );
			%
			if ( omega_trial >= omega )
				msg( thisFile, __LINE__, sprintf( "Step did not go downhill ( %0.3e >= %0.3e ) for pt %d.",  omega_trial, omega, ptIndex ) );
				deltaT *= brakefactor;
				continue;
			end
			%
			searchingForGoodStep = false;
			%
			% Instead of THIS, use PREDICTION based on vecRes and vecDeltaX.
			if ( norm(vecRes) < resTolLo * norm(vecDeltaX) )
				msg( thisFile, __LINE__, sprintf( "Step is excessively accurate ( %0.3e < %0.3e * %0.3e ) for pt %d.", norm(vecRes), resTolLo, norm(vecDeltaX), ptIndex ) );
				deltaT *= accelFactor;
			end
		end
		%
		vecX = vecX_trial;
		omega = omega_trial;
		vecG = vecG_trial;
		matX(:,ptIndex) = vecX;
		%msg( thisFile, __LINE__, sprintf( "Accepting point:  %3d,  %10.3e,  %10.3e,  %10.3e.", ptIndex, norm(vecDeltaX), norm(vecX-vecX0), omega ) );
		%
		%
		if ( norm(vecDeltaX) < deltaXNormMin )
			msg( thisFile, __LINE__, sprintf( "Step size was too small ( %0.3e < %0.3e ) for pt %d.", norm(vecDeltaX), deltaXNormMin, ptIndex ) );
			return;
		end
		if ( fevalCount >= fevalLimit )
			msg( thisFile, __LINE__, sprintf( "Reached feval limit ( %d >= %d ) for pt %d.",  fevalCount, fevalLimit, ptIndex ) );
			return;
		end
	end
	%
return;
end
