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
	%omegaTol = mygetfield( prm, "omegaTol", sqrt(eps)*(abs(omega0)+abs(omegaLim)) );
	omegaTol = mygetfield( prm, "omegaTol", 0.01*(abs(omega0)+abs(omegaLim))/1.5 );
	%omegaTol = mygetfield( prm, "omegaTol", 0.25*(abs(omega0)+abs(omegaLim))/1.5 );
	omegaTol/omega0
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
	muBumper = mygetfield( prm, "muBumper", sqrt(eps) * muScale );
	%muBumper = mygetfield( prm, "muBumper", muScale/100.0 );
	%muBumper = mygetfield( prm, "muBumper", 0.1 * muScale );
	assert( isrealscalar(muBumper) );
	assert( muBumper > 0.0 );
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
	% Get muCrit.
	matRInvH = matR\matH;
	rvecEigRInvH = eig(matRInvH);
	minEig = min(rvecEigRInvH);
	if ( minEig > muBumper )
		muCrit = 0.0;
	else
		muCrit = -minEig + muBumper;
	end
	msg_main( verbLev, thisFile, __LINE__, sprintf( "muCrit = %g.", muCrit ) );
	%
	mu = muCrit;
	vecDelta = -(matH + mu*matR)\vecG;
	omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
	msg_main( verbLev, thisFile, __LINE__, sprintf( "Calculation %d, %g, (%g), %g.", 0, mu, norm(vecDelta), omega ) );
	omegaCrit = omega;
	if ( omegaCrit > omegaLim - omegaTol )
		msg_main( verbLev, thisFile, __LINE__, "omega appears to be pos-semi-def (or very nearly so.)" );
		muLim = muCrit;
		msg_main( verbLev, thisFile, __LINE__, sprintf( "Returning %d, %g, (%g), %g.", 0, mu, norm(vecDelta), omega ) );
		datOut.omegaCrit = omegaCrit;
		datOut.muCrit = muCrit;
		datOut.muScale = muScale;
		return;
	end
	%
	gTRInvG = vecG' * ( matR \ vecG );
	assert( isrealscalar(gTRInvG) );
	assert( gTRInvG > 0.0 );
	mu = gTRInvG / ( omega0 - omegaLim );
	if ( mu < muCrit )
		mu = muCrit + muScale;
	end
	%
	trialCountLimit = mygetfield( prm, "trialCountLimit", 100 );
	assert( isrealscalar(trialCountLimit) );
	assert( trialCountLimit > 0 );
	%
	muLo = muCrit;
	omegaLo = omegaCrit;
	muHi = [];
	omegaHi = [];
	haveHi = false;
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
			datOut.omegaCrit = omegaCrit;
			datOut.muCrit = muCrit;
			datOut.muScale = muScale;
			return;
		end
		assert( trialCount < trialCountLimit );
		%
		if ( omega < omegaLim )
			muLo = mu;
			omegaLo = omega;
		else
			muHi = mu;
			omegaHi = mu;
			haveHi = true;
		end
		%
		vecDeltaPrime = -(matH + mu*matR)\(matR*vecDelta);
		omegaPrime = vecDeltaPrime'*vecG + vecDeltaPrime'*matH*vecDelta;
		mu_next = mu + (omega0-omega) * (omegaLim-omega) / ( (omega0-omegaLim) * omegaPrime );
		%
		muBondLo = ( mu + 3.0*muLo ) / 4.0;
		if ( mu_next < muBondLo )
			mu_next = muBondLo;
		elseif ( haveHi )
			muBondHi = ( mu + 3.0*muHi ) / 4.0;
			if ( mu_next > muBondHi )
				mu_next = muBondHi;
			end
		end
		mu = mu_next;
	end
	assert(0);
return;
%end

%!test
%!	assert(0);
