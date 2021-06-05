function [ mu, retCode, datOut ] = extFit_findMuOfOmega( omegaTarget, omega0, vecG, matH, matR, prm=[] )
	%
	thisFile = "extFit_findMuOfOmega";
	commondefs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	%
	size1 = size(matH,1);
	size2 = size(matH,2);
	assert( size1 == size2 );
	assert( isrealscalar(omegaTarget) );
	assert( isrealscalar(omega0) );
	assert( omegaTarget < omega0 );
	assert( isrealarray(vecG,[size1,1]) );
	assert( isrealarray(matH,[size1,size2]) );
	assert( isrealarray(matR,[size1,size2]) );
	%
	%
	omegaTol = mygetfield( prm, "omegaTol", sqrt(eps)*(abs(omega0)+abs(omegaTarget)) );
	assert( isrealscalar(omegaTol) );
	assert( omegaTol > 0.0 );
	assert( omegaTol < abs(omega0-omegaTarget) );
	%
	symTol = mygetfield( prm, "symTol", sqrt(eps) );
	assert( isrealscalar(symTol) );
	assert( symTol > 0.0 );
	%
	hScale = sqrt( sum(sum(matH.^2)) / (size1*size2) );
	rScale = sqrt( sum(sum(matR.^2)) / (size1*size2) );
	assert( rScale > 0.0 );
	muScale = hScale / rScale;
	epsMu = mygetfield( prm, "epsMu", sqrt(eps) * muScale );
	assert( isrealscalar(epsMu) );
	assert( epsMu > 0.0 );
	%
	% Require (and enforce) R and H sym.
	assert( sum(sum((matH'-matH).^2)) < (symTol*hScale*size1)^2 );
	assert( sum(sum((matR'-matR).^2)) < (symTol*rScale*size1)^2 );
	matR = 0.5*(matR'+matR);
	matH = 0.5*(matH'+matH);
	%
	% Require R pos-def.
	% There may be a better way to do this.
	rvecEigValsR = eig(matR);
	assert( isrealarray(rvecEigValsR) );
	assert( min(rvecEigValsR) > 0.0 );
	%
	% Confirm that omega is decreasing for a sufficiently small delta.
	assert( vecG'*(matR\vecG) > 0.0 );
	%
	%
	% Default return values.
	mu = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	% Find min eigenvalue.
	rvecLamda = eig(matH,matR); % Returns eig( matR^-1 * matH ), right?
	lambdaMin = min(rvecLamda);
	%
	% Confirm that we can hit below omegaTarget.
	mu = max([ -lambdaMin + epsMu, 0.0 ]);
	vecDelta = -(matH + mu*matR)\vecG;
	omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
	msg_main( verbLev, thisFile, __LINE__, sprintf( "Calculation %d, %g, (%g), %g.", 0, mu, norm(vecDelta), omega ) );
	%
	if ( abs(omega-omegaTarget) <= omegaTol )
		msg_main( verbLev, thisFile, __LINE__, "Success." );
		retCode = RETCODE__SUCCESS;
		return;
	elseif ( omega > omegaTarget )
		msg_main( verbLev, thisFile, __LINE__, "omegaTarget appears to be unreachable." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	%
	%
	% Prep for iterative solver.
	bracketBumper = mygetfield( prm, "bracketBumper", 0.01 );
	assert( isrealscalar(bracketBumper) );
	assert( bracketBumper > 0.0 );
	assert( bracketBumper <= 0.5 );
	%
	trialCountLimit = mygetfield( prm, "trialCountLimit", 100 );
	assert( isrealscalar(trialCountLimit) );
	assert( trialCountLimit > 0 );
	%
	muLo = mu;
	muHi = [];
	haveHi = false;
	%
	% Initial guess.
	f0 = vecG' * ( matR \ vecG );
	mu = f0 / ( omega0 - omegaTarget );
	if ( mu < muLo )
		mu = muLo + muScale;
	end
	%
	trialCount = 0;
	while (1)
		assert( mu > muLo );
		if ( haveHi )
			assert( mu < muHi );
		end
		trialCount++;
		vecDelta = -(matH + mu*matR)\vecG;
		omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
		msg_main( verbLev, thisFile, __LINE__, sprintf( "Calculation %d, %g, (%g), %g.", trialCount, mu, norm(vecDelta), omega ) );
		%
		if ( abs(omega-omegaTarget) <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, "Success." );
			retCode = RETCODE__SUCCESS;
			return;
		end
		if ( trialCount >= trialCountLimit )
			msg_main( verbLev, thisFile, __LINE__, "Giving up." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		if ( omega < omegaTarget )
			muLo = mu;
		else
			muHi = mu;
			haveHi = true;
		end
		%
		vecDeltaPrime = -(matH + mu*matR)\(matR*vecDelta);
		omegaPrime = vecDeltaPrime'*vecG + vecDeltaPrime'*matH*vecDelta;
		mu_next = mu + (omega0-omega) * (omegaTarget-omega) / ( (omega0-omegaTarget) * omegaPrime );
		%
		if (haveHi)
			muBracketLo = muLo + bracketBumper*(muHi-muLo);
			muBracketHi = muHi - bracketBumper*(muHi-muLo);
			mu = cap( mu_next, muBracketLo, muBracketHi );
		else
			assert( mu_next > mu );
			mu = mu_next;
		end
	end
end

%!test
%!	assert(0);
