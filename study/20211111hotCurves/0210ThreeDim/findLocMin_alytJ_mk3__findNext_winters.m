% This is (to be) with inter-Step data handling, and a revised algorithm.

function [ vecX_next, datOut ] = findLocMin_alytJ_mk3__findNext_winters( vecX0, vecF0, matJ0, funchF, prm=[], datIn=[] )
	%
	% Unpack input.
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if ( debugMode )
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		%
		vecF0_test = funchF(vecX0);
		assert( isrealarray(vecF0_test,[sizeF,1]) );
		assert( reldiff(vecF0,vecF0_test) < sqrt(eps) );
		clear vecF0_test;
		%
		epsX_test = 0.001;
		vecDelta_test = epsX_test*ones(sizeX,1);
		vecJF_testA = ( funchF( vecX0 + vecDelta_test ) - funchF( vecX0 - vecDelta_test ) ) / 2.0;
		vecJF_testB = matJ0*vecDelta_test;
		if ( sumsq(vecJF_testA-vecJF_testB) > sqrt(eps)*(sumsq(vecJF_testA)+sumsq(vecJF_testA)) )
			msg( __FILE__, __LINE__, "*** WARNING: Jacobian calculated by funchFJ appears to be incorrect. ***" );
		endif
		clear vecDelta_test;
		clear vecJF_testA;
		clear vecJF_testB;
	endif
	%
	omega0 = sumsq(vecF0,1)/2.0;
	if ( 0.0 == omega0 )
		msg( __FILE__, __LINE__, "omega0 is already zero." );
		return;
	end
	vecG0 = matJ0'*vecF0;
	if ( 0.0 == norm(vecG0) )
		msg( __FILE__, __LINE__, "vecG0 is already zero." );
		return;
	end
	matI = eye(sizeX,sizeX);
	matJTJ0 = matJ0'*matJ0;
	matK0 = mygetfield( prm, "matK0", zeros(sizeX,sizeX) );
	dTreg0 = mygetfield( prm, "dTreg", 100.0 * norm(vecF0) / sqrt(max(abs(diag(matJTJ0)))));
	if ( debugMode )
		assert( isrealarray(matK0,[sizeX,sizeX]) );
		assert( issymmetric(matK0) );
		assert( isrealscalar(dTreg0) );
	end
	%
	haveBSF = false;
	vecX_bsf = [];
	vecF_bsf = [];
	omega_bsf = [];
	havePrev = false;
	vecX_prev = [];
	vecF_prev = [];
	omega_prev = [];
	%
	dTreg = dTreg0; % Trust region size *IS* modified in this code.
	matJTJ = matJTJ0; % matJTJ is *NOT* modified in this code; J is assumed to be exact.
	matK = matK0; % matK *IS* modified in this code.
	doMainLoop = true;
	while ( doMainLoop )
		%
		% Find the vecX_trial corresponding to the min omega of the current model,
		% subject to the trust region.
		% But first, find a value of mu that works...
		% We could try to make use of a previous mu, if we have one, but, POITROME.
		matH = matJTJ + matK;
		hAbsMax = max(max(abs(matH)));
		if ( 0.0 == hAbsMax )
			msg( __FILE__, __LINE__, "matH is zero. This is unexpected, but can be handled." );
			mu = norm(vecG0)/omega0;
			matM = matH + mu*matI;
			[ matR, cholFlag ] = chol( matM );
			if ( 0~=cholFlag )
				msg( __FILE__, __LINE__, "Cholesky factorization failed for matH = 0 case with mu for omegaModel = 0. This should be impossible!" );
				return;
			endif
		else
			mu = 0.0;
			matM = matH + mu*matI;
			[ matR, cholFlag ] = chol( matM );
			if ( 0~=cholFlag )
				muReguCoeff = 1e-5;
				mu = muReguCoeff * hAbsMax;
				matM = matH + mu*matI;
				[ matR, cholFlag ] = chol( matM );
				if ( 0~=cholFlag )
					msgif( debugMode, __FILE__, __LINE__, "Hessian appears to have a negative eigenvalue." );
					msgif( debugMode, __FILE__, __LINE__, "Calling eig(). This may be slow. Faster approaches may be possible..." );
					msgif( debugMode, __FILE__, __LINE__, "  start with mu = upper bound for eigenvalue of H, and target omega = 0.0." );
					msgif( debugMode, __FILE__, __LINE__, "  increase mu exponentially until chol() works." );
					msgif( debugMode, __FILE__, __LINE__, "But, POITROME." );
					[ matPsi_eig, matLambad_eig ] = eig( matH );
					muCrit = -min(diag(matLambda_eig));
					assert( 0.0 < muCrit );
					mu = muCrit + muReguCoeff * ( muCrit + hAbsMax );
					matM = matH + mu*matI;
					matR = chol( matM ); % Can but shouldn't throw an error.
					cholFlag = 0; % Since we got here.
				endif
			endif
		endif
		vecDelta = -( matR \ ( matR' \ vecG0 ) )
		deltaNorm = norm(vecDelta);
		%
		%
		if ( norm(vecDelta) > dTreg )
			treg_relTol = 0.3;
			treg_iterLimit = 10;
			haveMuHi = false;
			muHi = [];
			muLo = mu;
			treg_iterCount = 0;
			while ( norm(vecDelta) > dTreg || norm(vecDelta) < (1.0-treg_relTol)*dTreg )
				treg_iterCount++;
				if ( treg_iterCount > treg_iterLimit )
					msg( __FILE__, __LINE__, "Failed to satisfy trust region constraint in allowed number of iterations." );
					break;
				end
				if ( norm(vecDelta) > dTreg )
					muLo = mu;
				else
					muHi = mu;
					haveMuHi = true;
				endif
				% Model: ||delta||^2 = ( a / (b+mu) )^2,
				%  match ||delta||^2 and d/dmu(||delta||^2) at previous mu.
				vecDeltaPrime = -( matR \ ( matR' \ vecDelta ) );
				dsq = sumsq(vecDelta,1);
				ddsqdmu = 2.0*(vecDelta'*vecDeltaPrime);
				assert( 0.0 > ddsqdmu );
				b = 2.0*dsq/(-ddsqdmu) - mu;
				a = 2.0*(dsq^1.5)/(-ddsqdmu);
				dTrgt = (1.0-0.5*treg_relTol)*dTreg;
				mu = (a/dTrgt) - b;
				if ( haveMuHi )
					% Limit change if have muHi.
					if ( mu > muHi )
						mu = 0.9*muHi + 0.1*muLo;
					elseif ( mu < muLo )
						mu = 0.9*muLo + 0.1*muHi;
					endif
				endif
				assert( mu > muLo );
				% Update matM, matR, and vecDelta.
				matM = matH + mu*matI;
				matR = chol( matM ); % Can but shouldn't throw an error.
				vecDelta = -( matR \ ( matR' \ vecG ) );
			endwhile
		endif
		%
		%
		%
		omegaModel = omega0 + vecG0'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
		if ( omegaModel < -eps*omega0 )
			% Make omegaModel non-negative.
			omnn_tol = 0.3;
			omnn_iterLimit = 10;
			haveMuHi = false;
			muHi = [];
			omegaModelOfMuHi = [];
			muLo = mu;
			omegaModelOfMuLo = omegaModel;
			omnn_iterCount = 0;
			while ( omegaModel < -eps*omega0 || omegaModel > omegaModelPositiveTol*omega0 )
				omnn_iterCount++;
				if ( omnn_iterCount > omnn_iterLimit )
					msg( __FILE__, __LINE__, "Failed to satisfy omegaModel non-negative constraint in allowed number of iterations." );
					break;
				end
				if ( omegaModel < -eps*omega0 )
					muLo = mu;
					omgeaModelOfMuLo = omegaModel;
				else
					muHi = mu;
					omegaModelOfMuHi = omegaModel;
					haveMuHi = true;
				endif
				omegaTrgt = 0.5*omnn_tol*omega0;
				if ( ~haveMuHi )
					% Model: omegaModel = omega - ( g^2 / (c+mu) ).
					%  match omega at muPrev.
					mu = mu + sumsq(vecG0)*(omegaTrgt-omegaModel)/((omega0-omegaTrgt)*(omega0-omegaModel));
					assert( mu > muLo );
				else
					% Just use a linear model.
					mu = muLo + (muHi-muLo)*(omegaTrgt-omegaModelOfMuLo)/(omegaModelOfMuHi-omegaModelOfMuLo);
					assert( mu > muLo );
					assert( mu < muHi );
				end
				% Update matM, matR, and vecDelta.
				matM = matH + mu*matI;
				matR = chol( matM ); % Can but shouldn't throw an error.
				vecDelta = -( matR \ ( matR' \ vecG ) );
			endwhile
			dTreg = norm(vecDelta);
		endif
		%
		%
		%
		msg( __FILE__, __LINE__, "WIP RETURN!" ); vecX_next = vecX0 - 0.0001*vecG0; return;
		
		
		
		% Now we have our suggested step.
		% Is it worth doing a feval?
		
		
		error( "FOLLOWING CODE IS PRE-WIP." );
		
		doStepGenLoop = true;
		while ( doStepGenLoop )
			% If ||deltaX|| is very large compared to ||g||/||H||, issue a warning.
			%
			% If the newly generated step is larger than the previous step, issue a warning.
			%
			doStepGenLoop = false;
		endwhile
		%
		%
		% Determine whether or not its worth doing a feval.
		% If our new vecX_trial is sufficientlly close to our best trial,
		% or omegaModel_trial wouldn't be much better than omegaBest,
		% then, just accept the best and return.
		% TO-DO!
		% if we have a BSF and Model_trial is worse, flag a warning and return.
		% if we have no BSF, compare omegaModel_trial to omega0
		%   ( omegaModel_trial > 0.999 * omega0 ) % Reduction in omega is too small to be worthwhile.
		% if we have a BSF, compare omegaModel_trial to omega_bsf.
		%
		% Do a feval.
		vecF_trial = funchF( vecX_trial );
		%
		% If feval failed, cut trust region size.
		if ( ~isrealarray(vecF_trial,[sizeF,1]) )
			msgif( debugMode, __FILE__, __LINE__, "Function evaluation failed." );
			if (havePrevTrial)
				msgif( __FILE__, __LINE__, "Function evaluation failed even though it succeeded for an earlier trial step!" );
				msgif( __FILE__, __LINE__, "This suggests the invalid space is not simply at distant values of x." );
			endif
			if ( 0.0 >= coeff_reduceTregOnFevalFail )
				msgif( debugMode, __FILE__, __LINE__, "Bailing because feval failed and 0.0 >= coeff_reduceTregOnFevalFail." );
				return;
			endif
			dTreg = coeff_reduceTregOnFevalFail * norm(vecDelta);
			continue;
		endif
		%
		% If vecF_trial is radically different from vecFModel_trial, bail.
		% We could check to see how close vecF_trial is to (1-b)*vecF0 + b*vecFModel_trial,
		%  for the value of b in [0,1] that minimizes the difference, but,
		%  assuming norm( vecFModel_trial - vecF0 ) <= norm( vecF0 ),
		%  and coeff_declareModelIsRadicallyWrong is > 2 ish(?), this isn't necessary.
		if ( norm( vecF_trial - vecFModel_trial ) >= coeff_declareModelIsRadicallyWrong * norm(vecF0) )
			msgif( debugMode, __FILE__, __LINE__, "Function was radically different at trial point." );
			if (havePrevTrial)
				msgif( __FILE__, __LINE__, "Function was radically different even though it was reasonable for an earlier trial step!" );
				msgif( __FILE__, __LINE__, "This suggests the presence of a sepratrix." );
			endif
			if ( 0.0 >= coeff_reduceTregOnRadicallyWrong )
				msgif( debugMode, __FILE__, __LINE__, "Bailing because feval failed and 0.0 >= coeff_reduceTregOnFevalFail." );
				return;
			endif
			dTreg = coeff_reduceTregOnRadicallyWrong * norm(vecDelta);
			continue;
		endif
		%
		omega_trial = sumsq(vecF_trial,1)/2.0;
		if ( havePrevTrial )
		if ( omega_trial > omega_prevTrial )
			msg( __FILE__, __LINE__, "omega_trial > omega_prevTrial; this suggests a sepratrix." );
		endif
		endif
		%
		%
		%
		% Let's update our model with this new info.
		msgif( debugMode, __FILE__, __LINE__, "Updating K!" );
		matK += (2.0*(omega_trial-omegaModel_trial)/sumsq(vecDeltaX_trial)^2)*(vecDeltaX_trial*(vecDeltaX_trial'));
		%
		if ( ~haveBSF && omega_trial < omega0 )
			haveBSF = true;
			vecX_bsf = vecX_trial;
			vecF_bsf = vecF_trial;
			omega_bsf = omega_trial;
		elseif ( omega_trial < omega_bsf )
			haveBSF = true; % Redundant.
			vecX_bsf = vecX_trial;
			vecF_bsf = vecF_trial;
			omega_bsf = omega_trial;
		endif
		%
		% Under various circumstances, we could know that our trial point would be the next trial point,
		% so there's no need to calculate a new trial.
		% But, POITROME.
		havePrevTrial = true;
		vecX_prevTrial = vecX_trial;
		vecF_prevTrial = vecF_trial;
		omega_prevTrial = omega_trial;
		continue;
	endwhile
	
	
	
	error( "Reached junk code." );
	
	
	
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
		if ( omega_trial < omegaModelMin )
			msgif( omega_trial < omegaModelMin, __FILE__, __LINE__, "Trial is better than model min!!!" );
		elseif ( omega_trial < omegaModel_trial )
			msgif( omega_trial < omegaModel_trial, __FILE__, __LINE__, "Trial is better than model!" );
		endif
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

endfunction
