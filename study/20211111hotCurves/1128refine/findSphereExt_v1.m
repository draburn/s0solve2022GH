%function vecX = findSphereExt_v1( funchG, vecXC, bigR, vecX0, prm=[] )
	thisFile = "findSphereExt_v1";
	%
	sizeX = size(vecXC,1);
	assert( isrealarray(vecXC,[sizeX,1]) );
	assert( isrealscalar(bigR) );
	assert( 0.0 < bigR );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	iterLimit = mygetfield( prm, "iterLimit", 10000 );
	trialLimit = mygetfield( prm, "trialLimit", 40 );
	mu0 = mygetfield( prm, "mu0", 1e-4 );
	muFactor = mygetfield( prm, "muFactor", 2.0 );
	assert( isposintscalar(iterLimit) );
	assert( isposintscalar(trialLimit) );
	assert( isrealscalar(mu0) );
	assert( 0.0 < mu0 );
	assert( isrealscalar(muFactor) );
	assert( 1.0 < muFactor );
	normDeltaTol = 1e-4;
	normFallTol = 1e-8;
	%
	vecR0 = vecX0 - vecXC;
	if ( 0.0 == norm(vecR0) )
		vecG0 = funchG( vecX0 );
		assert( isrealarray(vecG0,[sizeX,1]) );
		if ( 0.0 == norm(vecG0) )
			error( "Initial guess is sphere center and initial gradient is zero." );
			return;
		end
		vecX0 = vecX0 + vecG0;
		vecR0 = vecX0 - vecXC;
		clear vecG0;
	end
	assert( 0.0 ~= norm(vecR0) );
	vecX0 = vecXC + (bigR * vecR0 / norm(vecR0));
	%
	%
	%
	funchGamma = @(vecXDummy)( findSphereExt__calcGamma( funchG, vecXC, bigR, vecXDummy, prm ) );
	vecX = vecX0;
	%
	vecGamma = funchGamma(vecX);
	iterCount = 0;
	while (1)
		figure(2)
		hold on;
		plot( vecX(1), vecX(2), 'k+', 'markersize', 10*(iterCount+1), 'linewidth', 3 );
		plot( vecX(1), vecX(2), 'ko', 'markersize', 10*(iterCount+1), 'linewidth', 3 );
		hold off;
		%
		iterCount++;
		if ( iterCount > iterLimit )
			return;
		end
		%
		matV = eye(sizeX,sizeX);
		for k=1:sizeX
			vecV = matV(:,k);
			assert( 0.0 ~= norm(vecV) );
			epsFD = 0.01;
			vecHV = ( funchGamma( vecX + epsFD*vecV ) - funchGamma( vecX - epsFD*vecV ) ) / ( 2.0*epsFD*norm(vecV) );
			matHV(:,k) = vecHV;
		end
		sizeK = sizeX;
		matVTHV = matV'*matHV
		%[ matR, cholFlag ] = chol(matVTV);
		%if ( 0 ~= cholFlag )
		%	matVTHV = eye(sizeK,sizeK)*mu0;
		%end
		vecVTGamma = matV'*vecGamma % Should always be norm(vecGamma)*vecE1.
		matIK = eye(sizeK,sizeK);
		%
		trialCount = 0;
		while (1)
			trialCount++;
			if ( trialCount >= trialLimit )
				msg( thisFile, __LINE__, "Failed to find an acceptable step." );
				return;
			elseif ( 0 == trialCount )
				mu = 0.0;
			elseif ( 1 == trialCount )
				mu = mu0;
			else
				mu *= muFactor;
			end
			echo__mu = mu
			%
			matM = matVTHV + (mu*matIK);
			[ matR, cholFlag ] = chol(matM);
			if ( 0 == cholFlag )
				vecDeltaPre_trial = -(matV*(matR\(matR'\vecVTGamma)))
				vecXPre_trial = vecX + vecDeltaPre_trial;
				vecRPre_trial = vecXPre_trial - vecXC;
				if ( 0.0 == norm(vecRPre_trial) )
					msg( thisFile, __LINE__, "Weird BT." );
					continue;
				end
				vecX_trial = vecXC + (bigR*vecRPre_trial/norm(vecRPre_trial));
				vecGamma_trial = funchGamma( vecX_trial );
				if ( norm(vecGamma_trial) < norm(vecGamma) )
					% Accept the step!
					msg( thisFile, __LINE__, "Ist good!" );
					break;
				end
				%
				msg( thisFile, __LINE__, "||g|| increased." );
			else
				msg( thisFile, __LINE__, "NPD." );
			end
		end
		%
		makeThisLastStep = false;
		if ( norm(vecX_trial-vecX) < normDeltaTol )
			makeThisLastStep = true;
		elseif ( norm(vecGamma_trial-vecGamma) < normFallTol )
			makeThisLastStep = true;
		end
		vecX = vecX_trial;
		vecGamma = vecGamma_trial;
		figure(2)
		hold on;
		plot( vecX(1), vecX(2), 'k+', 'markersize', 10*(iterCount+1), 'linewidth', 3 );
		plot( vecX(1), vecX(2), 'ko', 'markersize', 10*(iterCount+1), 'linewidth', 3 );
		hold off;
		if (makeThisLastStep)
			return;
		end
	end
return;
%end
