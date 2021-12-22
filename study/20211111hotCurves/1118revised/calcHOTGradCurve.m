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
	deltaXNormMin = mygetfield( prm, "deltaXNormMin", 0.01 );
	deltaXNormMax = mygetfield( prm, "deltaXNormMax", 0.1 );
	assert( 0.0 < deltaXNormMin );
	assert( deltaXNormMin < deltaXNormMax );
	%
	deltaT = mygetfield( prm, "deltaT", 1.9 * sqrt( deltaXNormMax * deltaXNormMin )  );
	assert( 0.0 < deltaT );
	%
	fevalLimit = mygetfield( prm, "fevalLimit", 10000 );
	assert( 0.0 < fevalLimit );
	%
	ptIndex = 1;
	vecX = vecX0;
	omega = omega0;
	vecG = vecG0;
	matX(:,1) = vecX0;
	%msg( thisFile, __LINE__, sprintf( "Accepting point:  %3d,  %10.3e,  %10.3e,  %10.3e.", ptIndex, 0.0, norm(vecX-vecX0), omega ) );
	%
	searchingForNextPt = true;
	while (searchingForNextPt)
		ptIndex++;
		%
		%
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
		%
		%
		if ( norm(vecDeltaX) < deltaXNormMin )
			msg( thisFile, __LINE__, sprintf( "Step size was too small ( %0.3e < %0.3e ) for pt %d.", norm(vecDeltaX), deltaXNormMin, ptIndex ) );
			return;
		end
		if ( norm(vecDeltaX) > deltaXNormMax )
			msg( thisFile, __LINE__, sprintf( "Step size was too large ( %0.3e > %0.3e ) for pt %d.", norm(vecDeltaX), deltaXNormMax, ptIndex ) );
			return;
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
			return;
		end
		%
		vecX = vecX_trial;
		omega = omega_trial;
		vecG = vecG_trial;
		matX(:,ptIndex) = vecX;
		%msg( thisFile, __LINE__, sprintf( "Accepting point:  %3d,  %10.3e,  %10.3e,  %10.3e.", ptIndex, norm(vecDeltaX), norm(vecX-vecX0), omega ) );
		%
		if ( fevalCount >= fevalLimit )
			msg( thisFile, __LINE__, sprintf( "Reached feval limit ( %d >= %d ) for pt %d.",  fevalCount, fevalLimit, ptIndex ) );
			return;
		end
	end
	%
return;
end
