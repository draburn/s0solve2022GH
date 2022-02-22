%msg( __FILE__, __LINE__, "WIP RETURN!" ); vecX_next = vecX - 0.0001*vecG; return;
%
%echo__vecX = vecX
%echo__vecG = vecG
matK = zeros(sizeX,sizeX);
trialCount = 0;
haveBestSoFar = false;
vecX_next = []; % Keep track of best so far.
vecF_next = []; % Keep track of best so far.
omega_next = []; % Keep track of best so far.
dGenMin = dTreg;
while (1)
	% Check pre-iter imposed limits.
	trialCount++;
	trialLimit = 100;
	if ( trialCount > trialLimit )
		msgif( debugMode, __FILE__, __LINE__, "Reached trialLimit." );
		return;
	endif
	%
	%
	%
	% Generate step
	genCount = 0;
	mu = 0.0;
	matH = matJTJ + matK;
	matH = (matH'+matH)/2.0;
	while (1)
		genCount++;
		genLimit = 100; % MAKE PARAM
		if ( genCount > genLimit )
			msgif( debugMode, __FILE__, __LINE__, "Reached genLimit." );
			return;
		endif
		%
		matM = matH + mu*eye(sizeX);
		[ matR, cholFlag ] = chol( matM );
		hScale = max(abs(diag(matH)));
		if ( 0~=cholFlag || min(abs(diag(matR))) <= 1e-4*max(abs(diag(matR))) )
			msgif( false, __FILE__, __LINE__, "Chol failed." );
			msgif( false, __FILE__, __LINE__, "We need a better way to set mu when chol fails." );
			msgif( false, __FILE__, __LINE__, "And, we need to check to ensure the calculated delta is good." );
			if ( 0.0 == mu )
				mu = 1e-2*hScale;
				% Using mu = sqrt(eps)*hScale causes numerical issues, caseNum 100.
				continue;
			endif
			msg( __FILE__, __LINE__, "ERROR: matH appears to have a negative eigenvalue; this should be impossible." );
			mu *= 2.0;
			mu += 0.001*hScale;
			continue;
		endif
		vecDelta_trial = -( matR \ ( matR' \ vecG ) );
		if ( false && reldiff(vecDelta_trial(1),vecDelta_trial(2)) > 0.1 )
			msg( __FILE__, __LINE__, "Looks like matH is poorly behaved!" );
			echo__vecX = vecX
			echo__vecF = vecF
			echo__matJ = matJ
			echo__omega = omega
			echo__vecG = vecG
			echo__matJTJ = matJTJ
			echo__matK = matK
			echo__matH = matH
			echo__mu = mu
			echo__matM = matM
			echo__matR = matR
			echo__vecDelta_trial = vecDelta_trial
			echo__Hdelta = matH*vecDelta_trial
			error( "Looks like matH is poorly behaved." );
		end
		omegaModelMin = omega + vecG'*vecDelta_trial + 0.5*(vecDelta_trial'*matH*vecDelta_trial);
		%
		% Note that dGenMin starts at dTreg.
		if ( norm(vecDelta_trial) >= dGenMin )
			% Increase mu to be under trust region / dGenMin.
			muPrev = mu;
			% Model: ||delta||^2 = ( a / (b+mu) )^2,
			%  match ||delta||^2 and d/dmu(||delta||^2) at muPrev.
			vecDeltaPrime_trial = -( matR \ ( matR' \ vecDelta_trial ) );
			dsq = sumsq(vecDelta_trial,1);
			ddsqdmu = 2.0*(vecDelta_trial'*vecDeltaPrime_trial);
			assert( 0.0 > ddsqdmu );
			b = 2.0*dsq/(-ddsqdmu) - muPrev;
			a = 2.0*(dsq^1.5)/(-ddsqdmu);
			% We'll target about half of dGenMin.
			dTrgt = 0.5*dGenMin;
			mu = (a/dTrgt) - b;
			assert( mu > muPrev );
			% This could possibly over backtrack, but, we won't worry about that case.
			% A linear model (||delta||^2 = a + b*mu) would avoid over backtracking,
			%  but, the convergence there may be very slow.
			continue;
		endif
		dGenMin = norm(vecDelta_trial);
		%
		omegaModel_trial = omega + vecG'*vecDelta_trial + 0.5*(vecDelta_trial'*matH*vecDelta_trial);
		if ( omegaModel_trial < -eps*abs(omega) )
			% Increase mu so that omegaModel is not negative.
			muPrev = mu
			% Model: omegaModel = omega - ( g^2 / (c+mu) ).
			%  match omega at muPrev.
			omegaTrgt = 0.5*omega;
			mu = muPrev + sumsq(vecG)*(omegaTrgt-omegaModel_trial)/((omega-omegaTrgt)*(omega-omegaModel_trial));
			assert( mu > muPrev );
			% Could this overbacktrack? I don't know.
			continue;
		endif
		%
		%
		% Check if step is too small.
		% Should we also check omegaFallTol and stepSizeTol???
		wolfeC = 1.0; % MAKE PARAM
		if ( abs(vecDelta_trial'*(vecG+matH*vecDelta_trial)) >= wolfeC * abs(vecDelta_trial'*vecG) )
			msgif( true, __FILE__, __LINE__, "Step is too small per Wolfe curvature condition. Giving up." );
			msgif( true, __FILE__, __LINE__, "Is this happens frequently, consider modifying how mu is calculated" );
			msgif( true, __FILE__, __LINE__, " while maintaining other constraints." );
			fooLHS = vecDelta_trial'*(vecG+matH*vecDelta_trial)
			fooRHS = vecDelta_trial'*vecG
			return;
		endif
		%
		break;
	endwhile
	vecFModel_trial = vecF + matJ*vecDelta_trial; % No quad info here.
	omegaModel_trial = omega + vecG'*vecDelta_trial + 0.5*(vecDelta_trial'*matH*vecDelta_trial); % Repeated. POITROME.
	vecGModel_trial = vecG + matH*vecDelta_trial; % Some quad info here.
	%
	%
	%
	vecX_trial = vecX + vecDelta_trial;
	vecF_trial = funchFJ( vecX_trial );
	fevalCount++;
	if ( ~isrealarray(vecF_trial,[sizeF,1]) )
		msgif( debugMode, __FILE__, __LINE__, "Function evaluation failed." );
		if (haveBestSoFar)
			msgif( __FILE__, __LINE__, "Function evaluation failed even though it succeeded for an earlier trial step!" );
			msgif( __FILE__, __LINE__, "This suggests the invalid space is not simply at distant values of x." );
		endif
		dGenFailureBTFactor = 0.1; % MAKE PARAM.
		dGenMin *= dGenFailureBTFactor;
		msgif( debugMode, __FILE__, __LINE__, "Updating dTreg." );
		dTreg = dGenMin;
		continue;
	endif
	omega_trial = sumsq(vecF_trial,1)/2.0;
	%
	%
	%
	% Find how far vecF_trial is from vecF, vecFModel.
	% Since ||vecFModel-vecF|| <= 2.0*||vecF|| is guaranteed,
	%  we can simply compare against vecF.
	c_hugeChangeToF = 10.0; % MAKE PARAM.
	assert( c_hugeChangeToF >= 2.0 );
	if ( norm(vecF_trial-vecF) >= c_hugeChangeToF*norm(vecF) )
		msgif( debugMode, __FILE__, __LINE__, "Function was radically different at trial point." );
		if (haveBestSoFar)
			msgif( __FILE__, __LINE__, "Function was radically different even though it was reasonable for an earlier trial step!" );
			msgif( __FILE__, __LINE__, "This suggests the presence of a sepratrix." );
		endif
		dGenRadicallyWrongBTFactor = 0.1; % MAKE PARAM.
		dGenMin *= dGenRadicallyWrongBTFactor;
		dTreg = dGenMin;
		continue;
	endif
	%
	%
	% Having gotten here, the step is reasonable enough that we can use the results to update matK.
	if ( omega_trial >= omega )
		msgif( debugMode, __FILE__, __LINE__, "omega_trial >= omega" );
		msgif( debugMode, __FILE__, __LINE__, "Updating K!" );
		matK += (2.0*(omega_trial-omegaModel_trial)/sumsq(vecDelta_trial)^2)*(vecDelta_trial*(vecDelta_trial'));
		continue;
	endif
	%
	%
	%
	% Check outright acceptance criterion.
	c_acceptOutright = 0.3; % MAKE PARAM
	omegaAcceptOutright = (1.0-c_acceptOutright)*omega + c_acceptOutright*omegaModelMin;
	if ( omega_trial <= omegaAcceptOutright )
		msgif( debugMode, __FILE__, __LINE__, "Step is very good!" );
		msgif( omega_trial < omegaModel_trial, __FILE__, __LINE__, "Trial is better than model!" );
		msgif( omega_trial < omegaModelMin, __FILE__, __LINE__, "Trial is better than model min!!!" );
		if (haveBestSoFar)
			assert( ~isempty(omega_trial) );
			if ( omega_trial >= omega_next )
				error( "Somehow this very good step is no better than an earlier one. This should be impossible." );
			endif
		endif
		% Update to matK won'd be used, but, let's do it anyway.
		matK += (2.0*(omega_trial-omegaModel_trial)/sumsq(vecDelta_trial)^2)*(vecDelta_trial*(vecDelta_trial'));
		haveBestSoFar = true; % Though, this result won't be used anywhere.
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		omega_next = omega_trial;
		return;
	endif
	%
	%msg( __FILE__, __LINE__, "WIP RETURN!" ); vecX_next = vecX + vecDelta_trial; return;
	%
	% Is this the best result so far?
	if ( ~haveBestSoFar )
		haveBestSoFar = true;
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		omega_next = omega_trial;
	elseif ( omega_trial > omega_next )
		msgif( debugMode, __FILE__, __LINE__, "omega_trial > omega_next. This suggests a sepratrix." );
		% We'll accept the best so far (==next) if it's below a loose threshold...
		c_acceptReluctant = 1e-3; % MAKE PARAM
		vecDelta_next = vecX_next - vecX;
		omegaModel_next = omega + vecG'*vecDelta_next + 0.5*(vecDelta_next'*matH*vecDelta_next); % Repeated. POITROME.
		% Note that omegaModel_next uses the current K, not the one used to calc _next.
		if ( omega_next <= (1.0-c_acceptReluctant)*omega + c_acceptReluctant*omegaModel_next )
			msg( __FILE__, __LINE__, "Omega increased while backtracking and previous(?) trial is good enough." );
			return;
		end
	elseif ( omega_trial < omega_next )
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		omega_next = omega_trial;
	endif
	%
	msgif( debugMode, __FILE__, __LINE__, "Updating K!" );
	matK += (2.0*(omega_trial-omegaModel_trial)/sumsq(vecDelta_trial)^2)*(vecDelta_trial*(vecDelta_trial'));
endwhile
