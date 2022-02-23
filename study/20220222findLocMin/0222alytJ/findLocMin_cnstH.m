% Function...
%  [ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm=[], datIn=[] )
% Returns vecDelta corresponding to the local min of the omega model,
%  possibly subject to a trust region and a "reasonableness" constraint (such as omega >= 0).

function [ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm=[], datIn=[] )
	%
	%
	% Parse input.
	sizeX = size(vecG,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if ( debugMode )
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG,[sizeX,1]) );
		assert( isrealarray(matH,[sizeX,sizeX]) );
		assert( issymmetric(matH) );
	endif
	gNormSq = sumsq(vecG);
	gNorm = sqrt(gNormSq);
	hNorm = sqrt(sum(sumsq(matH)));
	matI = eye(sizeX,sizeX);
	%
	muReguCoeff = muReguCoeff( prm, "muReguCoeff", 1.0e-5 );
	deltaNormMax = mygetfield( prm, "deltaNormMax", 100.0*gNorm/hNorm );
	deltaNormMaxRelTol = mygetfield( prm, "deltaNormMaxRelTol", 0.4 );
	omegaModelMin = mygetfield( prm, "omegaModelMin", [] );
	omegaModelMinRelTol = mygetifled( prm, "omegaModelMinRelTol", 0.4 );
	useDeltaNormMax = ~isempty(deltaNormMax);
	useOmegaModelMin = ~isempty(omegaModelMin);
	if ( debugMode )
		assert( isrealscalar(muReguCoeff) );
		assert( 0.0 < muReguCoeff );
		if ( useDeltaNormMax )
			assert( isrealscalar(deltaNormMax) );
			assert( 0.0 < deltaNormMax );
			assert( isrealscalar(deltaNormMaxRelTol) );
			assert( 0.0 < deltaNormMaxRelTol );
			assert( deltaNormMaxRelTol < 1.0 );
		endif
		if ( useOmegaModelMin )
			assert( isrealscalar(omegaModelMin) );
			assert( isrealscalar(omegaModelMinRelTol) );
			assert( 0.0 < omegaModelMinRelTol );
			assert( omegaModelMinRelTol < 1.0 );
		elseif
	endif
	%
	%
	% Set default return values.
	vecDelta = zeros(sizeX,1);
	if (nargout>=2)
		datOut = [];
	endif
	%
	%
	% Handle "corner zero" cases.
	if ( useOmegaModelMin )
	if ( omegaModelMin >= omega0 )
		msgif( debugMode, __FILE__, __LINE__, "Initial objective value is already below specified minimum." );
		return;
	endif
	endif
	if ( 0.0 == gNorm )
		msgif( debugMode, __FILE__, __LINE__, "Initial gradient is zero." );
		return;
	endif
	if ( 0.0 == hNorm )
		if ( useDeltaNormMax && useOmegaModelMin  )
			vecDelta = -min([ deltaNormMax/gNorm, (omega0-omegaModelMin)/gNormSq ]) * vecG;
		elseif ( useDeltaNormMax )
			vecDelta = -(deltaNormMax/gNorm) * vecG;
		elseif ( useOmegaModelMin )
			vecDelta = -((omega0-omegaModelMin)/gNormSq) * vecG;
		else
			msgif( debugMode, __FILE__, __LINE__, "Hessian is zero, gradient is not, and there is no constraint on step." );
		endif
		return;
	endif
	%
	%
	% Find the largest valid step we can (or, at least, something close).
	% Consider letting an inital guess for mu be passed in.
	mu = 0.0;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 ~= cholFlag )
		msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with mu zero failed; trying with a small mu." );
		mu = muReguCoeff*hMaxAbs;
		matM = matH + mu*matI;
		[ matR, cholFlag ] = chol( matM );
		if ( 0 ~= cholFlag )
			msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with small mu failed; finding muCrit." );
			msgif( debugMode, __FILE__, __LINE__, "Calling eig(). This may be slow. Faster approaches may be possible, such as:" );
			msgif( debugMode, __FILE__, __LINE__, "  start with mu = upper bound for eigenvalue of H, and target omega = 0.0; or," );
			msgif( debugMode, __FILE__, __LINE__, "  increase mu exponentially until chol() works." );
			[ matPsi_eig, matLambda_eig ] = eig( matH );
			muCrit = -min(diag(matLambda_eig));
			mu = muCrit + muReguCoeff * ( muCrit + hNorm );
			matM = matH + mu*matI;
			matR = chol( matM );
			%
			clear matPsi_eig;
			clear matLambda_eig;
			clear muCrit;
		endif
	endif
	vecDelta = -( matR \ (matR'\vecG) );
	%
	%
	%
	% Apply constraints...
	% As we apply a constraint and backtrack, keep track of previous value of mu (which doesn't satisfy the constraint)
	% so we can forwardtrack later in case we overshot the constraint.
	% We'll keep track of quantities that require some effort to calculate:
	%  mu, matR, vecDelta, and (late) omegaModel.
	% We could track more, but, POITROME.
	%
	%
	% If necessary, backtrack to satisfy deltaNormMax (trust region).
	haveBTedForDeltaNormMax = false;
	if ( isrealscalar(deltaNormMax) )
	if ( norm(vecDelta) > deltaNormMax )
		assert( deltaNormMax > 0.0 );
		deltaNormMin = (1.0-deltaNormMaxRelTol)*deltaNormMax;
		%
		deltaNormTrgt = ( deltaNormMax + deltaNormMin ) / 2.0; % For first iteration.
		iterLimit = 10; % Arbitrary.
		iterCount = 0;
		while ( norm(vecDelta) > deltaNormMax )
			iterCount++;
			if ( iterCount > iterLimit )
				msg( __FILE__, __LINE__, "Failed to satisfy deltaNormMax (trust region constraint)." );
				return;
			endif
			%
			muLo = mu;
			matR_muLo = matR;
			vecDelta_muLo = vecDelta;
			%
			% Model: ||delta||^2 = ( a / (b+mu) )^2,
			%  match ||delta||^2 and d/dmu(||delta||^2) at previous mu.
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
			assert( mu > muLo );
			%
			haveBTedForDeltaNormMax = true;
			deltaNormTrgt = deltaNormMin; % Be more aggressive for subsequent iterations.
		endwhile
	endif
	endif
	%
	%
	%
	omegaModel = omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
	haveBTedForOmegaModelMin = false;
	if ( useOmegaModelMin )
	if ( omegaModel < omegaModelMin )
		% Since we apply omegaModelMin after deltaNormMax,
		%  and, since we've gotten here,
		% This means the deltaNormMax constraint is superfluous;
		%  only the omegaModelMin constraint matters.
		assert( omegaModelMin < omega0 );
		omegaModelMax = omegaModelMin + omegaModelMinRelTol*(omega0-omegaModelMin);
		%
		omegaTrgt = ( omegaModelMax + omegaModelMin ) / 2.0; % For first iteration.
		iterLimit = 10; % Arbitrary.
		iterCount = 0;
		while ( norm(vecDelta) > deltaNormMax )
			iterCount++;
			if ( iterCount > iterLimit )
				msg( __FILE__, __LINE__, "Failed to satisfy omegaModelMin (model reasonableness constraint)." );
				return;
			endif
			%
			muLo = mu;
			matR_muLo = matR;
			vecDetla_muLo = vecDelta;
			omegaModel_muLo = omegaModel;
			%
			% Model: omegaModel = omega0 - ( g^2 / (c+mu) ).
			%  match omega at previous mu.
			% Not a great model, but, should be good enough.
			mu = mu + normGSq*(omegaTrgt-omegaModel)/((omega0-omegaTrgt)*(omega0-omegaModel));
			matM = matH + mu*matI;
			matR = chol( matM );
			vecDelta = -( matR \ ( matR' \ vecG ) );
			omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
			assert( mu > muLo );
			%
			haveBTedForOmegaModelMin = true;
			omegaTrgt = omegaModelMax; % Be more aggressive for subsequent iterations.
		endwhile
	endif
	endif
	%
	%
	%
	% We possibly may have over backtracked.
	if ( haveBTedForOmegaModelMin )
		% Since omegaModelMin was applied after deltaNormMax,
		%  the deltaNormMax constraint is irrelevant.
		% mu should now satisfy omegaModelMin, and muLo should not.
		assert( omegaModel >= omegaModelMin );
		assert( omegaModel_muLo < omegaModelMin );
		% If mu also satisfies omegaModelMax, then we're done.
		% Otherwise, we'll forward track to get a mu that does satisfy omegaModelMax.
		if ( omegaModel > omegaModelMax )
			muHi = mu;
			matR_muHi = matR
			vecDelta_muHi = vecDelta;
			omegaModel_muHi = omegaModel;
			error( "Not implemented." );
		endif
	elseif (haveBTedForDeltaNormMax)
		% mu should now satisfy deltaNormMax, and muLo should not.
		assert( norm(vecDelta) >= deltaNormMax );
		assert( norm(vecDelta_muLo) < deltaNormMax );
		% If mu also satisfies deltaNormMin, then we're done.
		% Otherwise, we'll forward track to get a mu that does satisfy deltaNormMin.
		if ( norm(vecDelta) < deltaNormMin )
			muHi = mu;
			matR_muHi = matR
			vecDelta_muHi = vecDelta;
			omegaModel_muHi = omegaModel;
			error( "Not implemented." );
		endif
	endif
return;
endfunction
