function [ muLim, retCode, datOut ] = extFit_findMuLim( omegaLim, omega0, vecG, matH, matR, prm=[] )
	%
	% delta(mu) = -( matH + mu*matR ) \ vecG.
	% omega(mu) = omega0 + vecG' * delta + 0.5 * delta' * matH * delta.
	% If the curve would hit omegaLim before mu reaches zero, then,
	%  return the value of mu at which that happens;
	%  otherwise, return zero.
	% "R" for "regularization".
	%
	thisFile = "extFit_findMuLim";
	commondefs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	msg( thisFile, __LINE__, "DEPRECATED." );
	%
	sizeK = size(matH,1);
	sizeK = size(matH,2);
	assert( isrealscalar(omegaLim) );
	assert( isrealscalar(omega0) );
	assert( omegaLim < omega0 );
	assert( isrealarray(vecG,[sizeK,1]) );
	assert( isrealarray(matH,[sizeK,sizeK]) );
	assert( isrealarray(matR,[sizeK,sizeK]) );
	%
	%
	omegaTol = mygetfield( prm, "omegaTol", sqrt(eps)*(abs(omega0)+abs(omegaLim)) );
	%omegaTol = mygetfield( prm, "omegaTol", 0.01*(abs(omega0)+abs(omegaLim))/1.5 );
	%omegaTol = mygetfield( prm, "omegaTol", 0.1*(abs(omega0)+abs(omegaLim))/1.5 );
	%omegaTol/omega0
	assert( isrealscalar(omegaTol) );
	assert( omegaTol > 0.0 );
	assert( omegaTol < abs(omega0-omegaLim) );
	symTol = mygetfield( prm, "symTol", sqrt(eps) );
	assert( isrealscalar(symTol) );
	assert( symTol > 0.0 );
	%
	hScale = sqrt(sum(sum(matH.^2))) / sizeK;
	rScale = sqrt(sum(sum(matR.^2))) / sizeK;
	assert( rScale > 0.0 );
	muScale = hScale / rScale;
	epsMu = mygetfield( prm, "epsMu", sqrt(eps) * muScale );
	assert( isrealscalar(epsMu) );
	assert( epsMu > 0.0 );
	%
	%msg_main( verbLev, thisFile, __LINE__, "BUMPER AND TOLERANCE SHOULD BE IN TERMS OF STEP LENGTH, NOT MU!" );
	% Isn't using strict tolerances (/small bumper) good enough?
	%
	% Default return values.
	muLim = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	% Require (and enforce) R and H sym.
	assert( sum(sum((matH'-matH).^2)) < (symTol*hScale*sizeK)^2 );
	assert( sum(sum((matR'-matR).^2)) < (symTol*rScale*sizeK)^2 );
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
	%
	%
	msg( thisFile, __LINE__, "Need proper handling for case(s?) where borderline pos-def!" );
	% Get muCrit.
	matRInvH = matR\matH;
	rvecEigRInvH = eig(matRInvH);
	minEig = min(rvecEigRInvH);
	if ( minEig > epsMu )
		msg_main( verbLev, thisFile, __LINE__, "omega appears to be pos-def." );
		muCrit = 0.0;
		mu = 0.0;
		vecDelta = -(matH + mu*matR)\vecG;
		omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
		msg_main( verbLev, thisFile, __LINE__, sprintf( "Calculation %d, %g, (%g), %g.", 0, mu, norm(vecDelta), omega ) );
		%
		if ( omega > omegaLim + omegaTol )
			msg_main( verbLev, thisFile, __LINE__, "omegaLim appears to be unreachable." );
			muLim = mu;
			msg_main( verbLev, thisFile, __LINE__, sprintf( "Returning %d, %g, (%g), %g.", 0, mu, norm(vecDelta), omega ) );
			datOut.muCrit = muCrit;
			datOut.muScale = muScale;
			return;
		elseif ( abs(omega - omegaLim) <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, "Min of omega appears to be omegaLim." );
			muLim = mu;
			msg_main( verbLev, thisFile, __LINE__, sprintf( "Returning %d, %g, (%g), %g.", 0, mu, norm(vecDelta), omega ) );
			datOut.muCrit = muCrit;
			datOut.muScale = muScale;
			return;
		end
	else
		muCrit = -minEig;
	end
	msg_main( verbLev, thisFile, __LINE__, sprintf( "muCrit = %g.", muCrit ) );
	%
	muLo = muCrit;
	muHi = [];
	haveHi = false;
	%
	gTRInvG = vecG' * ( matR \ vecG );
	assert( isrealscalar(gTRInvG) );
	assert( gTRInvG > 0.0 );
	mu = gTRInvG / ( omega0 - omegaLim );
	if ( mu < muCrit )
		mu = muCrit + muScale;
	end
	%
	bracketBumper = mygetfield( prm, "bracketBumper", 0.1 );
	assert( isrealscalar(bracketBumper) );
	assert( bracketBumper > 0.0 );
	assert( bracketBumper <= 0.5 );
	%
	trialCountLimit = mygetfield( prm, "trialCountLimit", 100 );
	assert( isrealscalar(trialCountLimit) );
	assert( trialCountLimit > 0 );
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
		if ( abs(omega-omegaLim) < omegaTol )
			muLim = mu;
			msg_main( verbLev, thisFile, __LINE__, sprintf( "Returning %d, %g, (%g), %g.", trialCount, mu, norm(vecDelta), omega ) );
			datOut.muCrit = muCrit;
			datOut.muScale = muScale;
			return;
		end
		assert( trialCount < trialCountLimit );
		%
		if ( omega < omegaLim )
			muLo = mu;
		else
			muHi = mu;
			haveHi = true;
		end
		%
		vecDeltaPrime = -(matH + mu*matR)\(matR*vecDelta);
		omegaPrime = vecDeltaPrime'*vecG + vecDeltaPrime'*matH*vecDelta;
		mu_next = mu + (omega0-omega) * (omegaLim-omega) / ( (omega0-omegaLim) * omegaPrime );
		%
		if (haveHi)
			muBracketLo = (1.0-bracketBumper)*muLo + bracketBumper*muHi;
			muBracketHi = (1.0-bracketBumper)*muHi + bracketBumper*muLo;
			mu = cap( mu_next, muBracketLo, muBracketHi );
		else
			assert( mu_next > mu );
			mu = mu_next;
		end
	end
return;
end

%!test
%!	assert(0);
