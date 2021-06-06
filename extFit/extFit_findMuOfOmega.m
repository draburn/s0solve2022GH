function [ mu, retCode, datOut ] = extFit_findMuOfOmega( omegaTarget, omega0, vecG, matH, matR, prm=[] )
	%
	commondefs;
	thisFile = "extFit_findMuOfOmega";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	retCode = RETCODE__NOT_SET;
	%
	%
	% Validate input.
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
	% Require (and enforce) R and H sym.
	hScale = sqrt( sum(sum(matH.^2)) / (size1*size2) );
	rScale = sqrt( sum(sum(matR.^2)) / (size1*size2) );
	assert( rScale > 0.0 );
	muScale = hScale / rScale;
	datOut.muScale = muScale;
	%
	symTol = mygetfield( prm, "symTol", sqrt(eps) );
	assert( isrealscalar(symTol) );
	assert( symTol > 0.0 );
	%
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
	% So, the input looks good.
	% Let's start looking for our mu.
	omegaTol = mygetfield( prm, "omegaTol", sqrt(eps)*(abs(omega0)+abs(omegaTarget)) );
	assert( isrealscalar(omegaTol) );
	assert( omegaTol > 0.0 );
	assert( omegaTol < abs(omega0-omegaTarget) );
	%
	%
	% See if we can hit below omegaTarget.
	% Find min eigenvalue, which limits how low mu can go.
	rvecLamda = eig(matH,matR); % Returns eig( matR^-1 * matH ), right?
	lambdaMin = min(rvecLamda);
	datOut.lambdaMin = lambdaMin;
	msg_progress( verbLev, thisFile, __LINE__, sprintf( "lambdaMin = %f.", lambdaMin ) );
	%
	epsMu = sqrt(eps)*muScale + sqrt(eps)*max([-lambdaMin,0.0]);
	epsMu = mygetfield( prm, "epsMu", epsMu );
	assert( isrealscalar(epsMu) );
	assert( epsMu > 0.0 );
	%
	mu = max([ -lambdaMin + epsMu, 0.0 ]);
	vecDelta = -(matH + mu*matR)\vecG;
	omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
	msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
	  " %3d, %10.3e, (%10.3e), %11.3e, %11.3e.", 0, mu, norm(vecDelta), omega, omega-omegaTarget ) );
	%
	if ( abs(omega-omegaTarget) <= omegaTol )
		msg_main( verbLev, thisFile, __LINE__, sprintf( "Success. Returning mu = %10.3e.", mu ) );
		retCode = RETCODE__SUCCESS;
		return;
	elseif ( omega > omegaTarget )
		msg_notify( verbLev, thisFile, __LINE__, "Target appears to be unreachable." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	%
	%
	% Prep for iterative solver.
	trialCountLimit = mygetfield( prm, "trialCountLimit", 50 );
	assert( isrealscalar(trialCountLimit) );
	assert( trialCountLimit > 0 );
	%
	bracketBumper = mygetfield( prm, "bracketBumper", 0.05 );
	assert( isrealscalar(bracketBumper) );
	assert( bracketBumper > 0.0 );
	assert( bracketBumper <= 0.5 );
	%
	muLo = mu;
	muHi = [];
	haveHi = false;
	%
	% Initial guess, based on omega in lim mu -> inf.
	mu = ( vecG' * ( matR \ vecG ) ) / ( omega0 - omegaTarget );
	if ( mu < muLo )
		mu = muLo + muScale;
	end
	trialCount = 0;
	while (1)
		assert( mu > muLo );
		if ( haveHi )
			assert( mu < muHi );
		end
		trialCount++;
		vecDelta = -(matH + mu*matR)\vecG;
		omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
		msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
		  "%3d, %10.3e, (%10.3e), %11.3e, %11.3e.", trialCount, mu, norm(vecDelta), omega, omega-omegaTarget ) );
		%
		if ( abs(omega-omegaTarget) <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( "Success. Returning mu = %10.3e.", mu ) );
			retCode = RETCODE__SUCCESS;
			return;
		end
		if ( trialCount >= trialCountLimit )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached trialCountLimit (%d).", trialCountLimit )  );
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
		% Next guess based on omega = omega0 + A / ( mu + B ),
		% set to match current mu, omega, and dOmega/dMu.
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
