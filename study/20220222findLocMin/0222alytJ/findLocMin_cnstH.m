% Function...
%  [ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm=[] )
% Returns vecDelta corresponding to the local min of the omega model,
%  possibly subject to a trust region and a "reasonableness" constraint (such as omega >= 0).

function [ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm=[] )
	%
	%
	% Parse input.
	sizeX = size(vecG,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if ( debugMode )
		msg( __FILE__, __LINE__, "Using debugMode." );
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
	cdmlPrm = [];
	deltaNormMax = [];
	deltaNormMaxRelTol = 0.4;
	omegaModelMin = 0.0;
	omegaModelMinRelTol = 0.4;
	if ( ~isempty(prm) )
		cdmlPrm = mygetfield( prm, "cdmlPrm", cdmlPrm );
		deltaNormMax = mygetfield( prm, "deltaNormMax", deltaNormMax );
		deltaNormMaxRelTol = mygetfield( prm, "deltaNormMaxRelTol", deltaNormMaxRelTol );
		omegaModelMin = mygetfield( prm, "omegaModelMin", omegaModelMin );
		omegaModelMinRelTol = mygetfield( prm, "omegaModelMinRelTol", omegaModelMinRelTol );
	endif
	useDeltaNormMax = ~isempty(deltaNormMax);
	useOmegaModelMin = ~isempty(omegaModelMin);
	if ( debugMode )
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
		endif
	endif
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
		msgif( __FILE__, __LINE__, "WARNING: Initial objective value is already below specified minimum." );
		return;
	endif
	endif
	if ( 0.0 == gNorm )
		msgif( __FILE__, __LINE__, "WARNING: Initial gradient is zero." );
		return;
	endif
	if ( 0.0 == hNorm )
		stepConstrainedByOmegaModelMin = false;
		if ( useDeltaNormMax && useOmegaModelMin  )
			assert( deltaNormMax > 0.0 );
			assert( omega0 > omegaModelMin );
			muD = gNorm/deltaNormMax;
			muO = gNormSq/(omega0-omegaModelMin);
			if (muO>muD)
				mu = muO;
				stepConstrainedByOmegaModelMin = true;
			else
				mu = muD;
			endif
		elseif ( useDeltaNormMax )
			assert( deltaNormMax > 0.0 );
			mu = gNorm/deltaNormMax;
		elseif ( useOmegaModelMin )
			assert( omega0 > omegaModelMin );
			mu = gNormSq/(omega0-omegaModelMin);
			stepConstrainedByOmegaModelMin = true;
		else
			msgif( __FILE__, __LINE__, "WARNING: Hessian is zero, gradient is not, and there is no constraint on step." );
			return;
		endif
		vecDelta = vecG/(-mu);
		datOut.mu = mu;
		datOut.matR = sqrt(mu)*matI;
		datOut.omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
		datOut.trustRegionShouldBeUpdated = stepConstrainedByOmegaModelMin;
		return;
	endif
	%
	%
	%
	% Find the largest valid step we can (or, at least, something close).
	[ vecDelta, cdmlDatOut ] = calcDeltaMaxLev( vecG, matH, cdmlPrm );
	mu = cdmlDatOut.mu;
	matR = cdmlDatOut.matR;
	%
	%
	%
	% Apply constraints...
	% As we apply a constraint and backtrack, keep track of previous value of mu (which doesn't satisfy the constraint)
	%  so we can forwardtrack later in case we overshot the constraint.
	% We'll keep track of quantities that are useful: mu, matR, vecDelta, and (later) omegaModel.
	%
	% If necessary, backtrack to satisfy deltaNormMax (trust region).
	haveBTedForDeltaNormMax = false;
	if ( isrealscalar(deltaNormMax) )
	if ( norm(vecDelta) > deltaNormMax )
		deltaNormMin = (1.0-deltaNormMaxRelTol)*deltaNormMax;
		deltaNormTrgt = ( deltaNormMax + deltaNormMin ) / 2.0; % For first iteration.
		iterLimit = 10; % Arbitrary.
		iterCount = 0;
		while ( norm(vecDelta) > deltaNormMax )
			iterCount++;
			if ( iterCount > iterLimit )
				msg( __FILE__, __LINE__, "WARNING: Failed to satisfy deltaNormMax (trust region constraint)." );
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
			assert( mu > muLo );
			matM = matH + mu*matI;
			matR = chol( matM );
			vecDelta = -( matR \ ( matR' \ vecG ) );
			%
			haveBTedForDeltaNormMax = true;
			deltaNormTrgt = deltaNormMin; % Be more aggressive for subsequent iterations.
		endwhile
		if ( debugMode )
			clear deltaNormTrgt;
			clear iterLimit;
			clear iterCount;
			clear vecDeltaPrime;
			clear dsq;
			clear ddsqdmu;
			clear b;
			clear a;
		endif
	endif
	endif
	%
	%
	%
	omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
	haveBTedForOmegaModelMin = false;
	if ( useOmegaModelMin )
	if ( omegaModel < omegaModelMin )
		omegaModelMax = omegaModelMin + omegaModelMinRelTol*(omega0-omegaModelMin);
		omegaTrgt = ( omegaModelMax + omegaModelMin ) / 2.0; % For first iteration.
		iterLimit = 10; % Arbitrary.
		iterCount = 0;
		while ( omegaModel < omegaModelMin )
			iterCount++;
			if ( iterCount > iterLimit )
				msg( __FILE__, __LINE__, "WARNING: Failed to satisfy omegaModelMin (model reasonableness constraint)." );
				return;
			endif
			%
			muLo = mu;
			matR_muLo = matR;
			vecDelta_muLo = vecDelta;
			omegaModel_muLo = omegaModel;
			%
			% Model: omegaModel = omega0 - ( g^2 / (c+mu) ).
			%  match omega at previous mu.
			% Omega is harder to model than delta. This is not a great model, but, should be good enough.
			mu = mu + gNormSq*(omegaTrgt-omegaModel)/((omega0-omegaTrgt)*(omega0-omegaModel));
			assert( mu > muLo );
			matM = matH + mu*matI;
			matR = chol( matM );
			vecDelta = -( matR \ ( matR' \ vecG ) );
			omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
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
		% Since omegaModelMin was applied after deltaNormMax, the deltaNormMax constraint is now irrelevant.
		% muLo should not satisfy the constraint but mu should.
		assert( omegaModel_muLo < omegaModelMin );
		assert( omegaModel >= omegaModelMin );
		% If mu also satisfies omegaModelMax, then we're done.
		% Otherwise, we'll forward track to get a mu that does satisfy omegaModelMax.
		if ( omegaModel > omegaModelMax )
			muHi = mu;
			matR_muHi = matR;
			vecDelta_muHi = vecDelta;
			omegaModel_muHi = omegaModel;
			omegaModelTrgt = ( omegaModelMax + omegaModelMin ) / 2.0;
			iterLimit = 10;
			iterCount = 0;
			while ( 1 )
				iterCount++;
				if ( iterCount > iterLimit )
					msg( __FILE__, __LINE__, "WARNING: Failed to satisfy omegaModelMinRelTol." );
					mu = muHi;
					vecDelta = vecDelta_muHi;
					matR = matR_muHi;
					omegaModel = omegaModel_muHi;
					break;
				endif
				%
				mu = muLo + (muHi-muLo)*(omegaModelTrgt-omegaModel_muLo)/(omegaModel_muHi-omegaModel_muLo);
				assert( muLo < mu );
				assert( mu < muHi );
				matM = matH + mu*matI;
				matR = chol( matM );
				vecDelta = -( matR \ ( matR' \ vecG ) );
				omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
				%
				if ( omegaModel < omegaModelMin )
					assert( mu > muLo );
					muLo = mu;
					matR_muLo = matR;
					vecDelta_muLo = vecDelta;
					omegaModel_muLo = omegaModel;
					continue;
				elseif ( omegaModel > omegaModelMax )
					assert( mu < muHi );
					muHi = mu;
					matR_muHi = matR;
					vecDelta_muHi = vecDelta;
					omegaModel_muHi = omegaModel;
					continue;
				endif
				break;
			endwhile
		endif
	elseif ( haveBTedForDeltaNormMax )
		assert( ~haveBTedForOmegaModelMin );
		% muLo should not satisfy deltaNormMax but mu should.
		assert( norm(vecDelta_muLo) > deltaNormMax );
		assert( norm(vecDelta) <= deltaNormMax );
		% If mu also satisfies deltaNormMin, then we're done.
		% Otherwise, we'll forward track to get a mu that does satisfy deltaNormMin.
		if ( norm(vecDelta) < deltaNormMin )
			muHi = mu;
			matR_muHi = matR;
			vecDelta_muHi = vecDelta;
			%
			deltaNorm_muHi = norm(vecDelta_muHi);
			deltaNorm_muLo = norm(vecDelta_muLo);
			% Use secant method.
			deltaNormTrgt = ( deltaNormMax + deltaNormMin ) / 2.0
			iterLimit = 10;
			iterCount = 0;
			while ( 1 )
				iterCount++;
				if ( iterCount > iterLimit )
					msg( __FILE__, __LINE__, "WARNING: Failed to satisfy deltaNormMaxRelTol." );
					mu = muHi;
					vecDelta = vecDelta_muHi;
					matR = matR_muHi;
					omegaModel = omegaModel_muHi;
					break;
				endif
				%
				mu = muLo + (muHi-muLo)*(deltaNormTrgt-deltaNorm_muLo)/(deltaNorm_muHi-deltaNorm_muLo);
				%
				assert( muLo < mu );
				assert( mu < muHi );
				matM = matH + mu*matI;
				matR = chol( matM );
				vecDelta = -( matR \ ( matR' \ vecG ) );
				%
				deltaNorm = norm(vecDelta);
				if ( deltaNorm < deltaNormMin )
					assert( mu < muHi );
					muHi = mu;
					matR_muHi = matR;
					vecDelta_muHi = vecDelta;
					deltaNorm_muHi = deltaNorm;
					continue;
				elseif ( deltaNorm > deltaNormMax )
					assert( mu > muLo );
					muLo = mu;
					matR_muLo = matR;
					vecDelta_muLo = vecDelta;
					omegaModel_muLo = omegaModel;
					deltaNorm_muLo = deltaNorm;
					continue;
				endif
				break;
			endwhile
			omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
		endif
	endif
	%
	if ( nargout >= 2 )
		datOut.mu = mu;
		datOut.matR = matR;
		datOut.omegaModel = omegaModel;
		datOut.trustRegionShouldBeUpdated = haveBTedForOmegaModelMin;
	endif
return;
endfunction


%!test
%!	msg( __FILE__, __LINE__, "~~~ Positive Definite Test ~~~ " );
%!	omega0 = 1.0
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ OmegaModelMin Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	prm.omegaModelMinRelTol = 0.01;
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Sans OmegaModelMin Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	prm.omegaModelMin = [];
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ DeltaNormMax Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	prm.deltaNormMax = 0.3;
%!	prm.deltaNormMaxRelTol = 0.01;
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ DeltaNormMax Sans OmegaModelMin Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	prm.omegaModelMin = [];
%!	prm.deltaNormMax = 0.3;
%!	prm.deltaNormMaxRelTol = 0.01;
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Weak DeltaNormMax Sans OmegaModelMin Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	prm.omegaModelMin = [];
%!	prm.deltaNormMax = 0.5;
%!	prm.deltaNormMaxRelTol = 0.01;
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Weak DeltaNormMax Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = eye(2,2)
%!	prm = [];
%!	prm.omegaModelMin = 0.0;
%!	prm.omegaModelMinRelTol = 0.01;
%!	prm.deltaNormMax = 0.5;
%!	prm.deltaNormMaxRelTol = 0.01;
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Zero Hessian Test ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = zeros(2,2)
%!	prm = [];
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Zero Hessian Test With DeltaNormMax ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = zeros(2,2)
%!	prm = [];
%!	prm.deltaNormMax = 0.2;
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )


%!test
%!	msg( __FILE__, __LINE__, "~~~ Zero Hessian Test With DeltaNormMax Sans OmegaModelMin ~~~ " );
%!	omega0 = 0.3
%!	vecG = [ 1.0; 0.0 ]
%!	matH = zeros(2,2)
%!	prm = [];
%!	prm.deltaNormMax = 0.2;
%!	prm.omegaModelMin = [];
%!	echo__prm = prm
%!	[ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm )
