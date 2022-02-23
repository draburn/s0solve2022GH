	deltaNormTrgt = ( deltaNormMax + deltaNormMin ) / 2.0; % For first iteration.
	iterLimit = 10; % Arbitrary-ish.
	%
	iterCount = 0;
	while ( norm(vecDelta) > deltaNormMax )
		iterCount++;
		if ( iterCount > iterLimit )
			error( "Failed to satisfy trust region (cross deltaNormMax)." );
		endif
		%
		muLo = mu;
		%matR_muLo = matR;
		vecDelta_muLo = vecDelta;
		%
		vecDeltaPrime = -( matR \ ( matR' \ vecDelta ) );
		dsq = sumsq(vecDelta,1);
		ddsqdmu = 2.0*(vecDelta'*vecDeltaPrime);
		assert( 0.0 > ddsqdmu );
		b = 2.0*dsq/(-ddsqdmu) - mu;
		a = 2.0*(dsq^1.5)/(-ddsqdmu);
		%
		mu = (a/deltaNormTrgt) - b;
		matM = matH + mu*matI;
		matR = chol( matM );
		vecDelta = -( matR \ ( matR' \ vecG ) );
		%
		deltaNormTrgt = deltaNormMin;
	endwhile
	%
	assert( norm(vecDelta) <= deltaNormMax );
	if ( norm(vecDelta) >= deltaNormMin )
		return;
	endif
	%
	muHi = mu;
	%matR_muHi = matR;
	vecDelta_muHi = vecDelta;
	%
	iterLimit = 100; % Arbitrary-ish.
	iterCount = 0;
	while ( norm(vecDelta) > deltaNormMax || norm(vecDelta) < deltaNormMin )
		iterCount++;
		if ( iterCount > iterLimit )
			error( "Failed to satisfy trust region (cross deltaNormMax)." );
		endif
		%
		% Use a crude bisection. POITROME.
		mu = ( muLo + muHi ) / 2.0;
		matM = matH + mu*matI;
		matR = chol( matM );
		vecDelta = -( matR \ ( matR' \ vecG ) );
	endwhile
