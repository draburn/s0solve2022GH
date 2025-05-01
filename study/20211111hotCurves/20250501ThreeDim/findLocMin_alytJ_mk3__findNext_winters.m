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
	matI = eye(sizeX,sizeX);
	omega0 = sumsq(vecF0,1)/2.0;
	vecG0 = matJ0'*vecF0;
	matJTJ0 = matJ0'*matJ0;
	if ( 0.0 == norm(vecG0) )
		msg( __FILE__, __LINE__, "vecG0 is zero." );
		return;
	elseif ( 0.0 == omega0 )
		msg( __FILE__, __LINE__, "omega0 is below tolerance." );
		return;
	endif
	matK0 = mygetfield( prm, "matK0", zeros(sizeX,sizeX) );
	dTreg0 = mygetfield( prm, "dTreg", 100.0 * norm(vecG0) / max(abs(diag(matJTJ0))) );
	omegaTol = mygetfield( prm, "omegaRelTol", 1e-4*omega0 );
	if ( debugMode )
		assert( isrealarray(matK0,[sizeX,sizeX]) );
		assert( issymmetric(matK0) );
		assert( isrealscalar(dTreg0) );
		assert( isrealscalar(omegaTol) );
		assert( 0.0 <= omegaTol );
	end
	%
	if (debugMode)
		msg( __FILE__, __LINE__, "Input values..." );
		echo__matK0 = matK0
		echo__vecX0 = vecX0
		echo__vecF0 = vecF0
		echo__matJTJ0 = matJTJ0
		echo__omega0 = omega0
		echo__vecG0 = vecG0
		echo__dTreg0 = dTreg0
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
	matK = matK0;
	dTreg = dTreg0;
	doMainLoop = true;
	iterCount = 0;
	datOut.fevalCount = 0;
	while ( doMainLoop )
		%
		iterCount++;
		iterLimit = 100;
		if ( iterCount > iterLimit )
			msg( __FILE__, __LINE__, "Reached iterLimit." );
			doMainLoop = false;
			break;
		end
		%
		%
		%matH = matJTJ0 + matK;
		%%% [ vecDelta, ftdatOut ] = findLocMinJ__findTrial( omega0, vecG0, matH, dTreg, ftprm );
		%%% assert( norm(vecDelta) < dTreg );
		%%% omegaModel = omega0 + vecG0'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
		%%% assert( omegaModel >= 0.0 );
		%%% vecFModel = vecF0 + matJ0*vecDelta;
		%%% dTreg = ftdatOut.dTreg;
		%
		% Find the vecX_trial corresponding to the min omega of the current model,
		% subject to the trust region.
		% But first, find a value of mu that works...
		% We could try to make use of a previous mu, if we have one, but, POITROME.
		matH = matJTJ0 + matK;
		hAbsMax = max(max(abs(matH)));
		if ( 0.0 == hAbsMax )
			msg( __FILE__, __LINE__, "matH is zero. This is unexpected, but can be handled." );
			mu = norm(vecG0)/omega0;
			matM = matH + mu*matI;
			matR = chol( matM );
			cholFlag = 0; % Since chol() didn't error away.
		else
			mu = 0.0;
			matM = matH + mu*matI;
			msgif( debugMode, __FILE__, __LINE__, "Trying chol with mu = 0.0." );
			[ matR, cholFlag ] = chol( matM );
			if ( 0~=cholFlag )
				muReguCoeff = 1e-5;
				mu = muReguCoeff * hAbsMax;
				matM = matH + mu*matI;
				msgif( debugMode, __FILE__, __LINE__, "Trying chol with small mu for regularization." );
				[ matR, cholFlag ] = chol( matM );
				if ( 0~=cholFlag )
					msgif( debugMode, __FILE__, __LINE__, "Hessian appears to have a negative eigenvalue." );
					msgif( debugMode, __FILE__, __LINE__, "Calling eig(). This may be slow. Faster approaches may be possible..." );
					msgif( debugMode, __FILE__, __LINE__, "  start with mu = upper bound for eigenvalue of H, and target omega = 0.0." );
					msgif( debugMode, __FILE__, __LINE__, "  increase mu exponentially until chol() works." );
					msgif( debugMode, __FILE__, __LINE__, "But, POITROME." );
					[ matPsi_eig, matLambda_eig ] = eig( matH );
					muCrit = -min(diag(matLambda_eig));
					assert( 0.0 < muCrit );
					mu = muCrit + muReguCoeff * ( muCrit + hAbsMax );
					matM = matH + mu*matI;
					matR = chol( matM );
					cholFlag = 0; % Since chol() didn't error away.
				endif
			endif
		endif
		vecDelta = -( matR \ ( matR' \ vecG0 ) );
		deltaNorm = norm(vecDelta);
		%
		%
		%
		if ( norm(vecDelta) > dTreg )
			msgif( debugMode, __FILE__, __LINE__, sprintf( "mu = %g, norm(vecDelta) = %g, dTreg = %g.", mu, norm(vecDelta), dTreg ) );
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
					if ( ~haveMuHi )
						mu = muHi;
					end
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
				vecDelta = -( matR \ ( matR' \ vecG0 ) );
				msgif( debugMode, __FILE__, __LINE__, sprintf( "mu = %g, norm(vecDelta) = %g, dTreg = %g.", mu, norm(vecDelta), dTreg ) );
			endwhile
		endif
		%
		%
		%
		omegaModel = omega0 + vecG0'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
		if ( omegaModel < -eps*omega0 )
			msgif( debugMode, __FILE__, __LINE__, sprintf( "mu = %g, norm(vecDelta) = %g, omegaModel = %g.", mu, norm(vecDelta), omegaModel ) );
			% Make omegaModel non-negative.
			omnn_tol = 0.3;
			omnn_iterLimit = 10;
			haveMuHi = false;
			muHi = [];
			omegaModelOfMuHi = [];
			muLo = mu;
			omegaModelOfMuLo = omegaModel;
			omnn_iterCount = 0;
			while ( omegaModel < -eps*omega0 || omegaModel > omnn_tol*omega0 )
				omnn_iterCount++;
				if ( omnn_iterCount > omnn_iterLimit )
					msg( __FILE__, __LINE__, "Failed to satisfy omegaModel non-negative constraint in allowed number of iterations." );
					mu = muLo; % Omega model there is negative, but, okay.
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
					%%%% Just use a linear model.
					mu = muLo + (muHi-muLo)*(omegaTrgt-omegaModelOfMuLo)/(omegaModelOfMuHi-omegaModelOfMuLo);
					assert( mu > muLo );
					assert( mu < muHi );
					% Actually, bisection looks faster.
					%%%mu = (muLo+muHi)/2.0;
				end
				% Update matM, matR, and vecDelta.
				matM = matH + mu*matI;
				matR = chol( matM ); % Can but shouldn't throw an error.
				vecDelta = -( matR \ ( matR' \ vecG0 ) );
				omegaModel = omega0 + vecG0'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
				msgif( debugMode, __FILE__, __LINE__, sprintf( "mu = %g, norm(vecDelta) = %g, omegaModel = %g.", mu, norm(vecDelta), omegaModel ) );
			endwhile
			dTreg = norm(vecDelta); % This isn't strictly necessary?
		endif
		if ( omegaModel < 0.0 )
			omegaModel = 0.0;
		end
		vecFModel = vecF0 + matJ0*vecDelta;
		%
		%
		%
		% Decide whether or not to do a feval.
		if (haveBSF)
			if ( omegaModel >= omega_bsf )
				msgif( debugMode, __FILE__, __LINE__, sprintf( "Model is worse than best so far ( %g >= %g ).", omegaModel, omega_bsf ) );
				doMainLoop = false;
				break;
			end
			%%% THIS CRITERIA IS NOT VERY SENSIBLE!
			if ( abs(omegaModel-omega_bsf) <= omegaTol )
				msgif( debugMode, __FILE__, __LINE__, sprintf( "Model offers insufficient improvement over best so far ( %g vs %g ).", omegaModel, omega_bsf ) );
				doMainLoop = false;
				break;
			end
		end
		if ( abs(omegaModel-omega0) <= omegaTol )
			msgif( debugMode, __FILE__, __LINE__, sprintf( "Model offers insufficient improvement ( %g vs %g ).", omegaModel, omega0 ) );
			doMainLoop = false;
			break;
		end
		%
		%
		%
		vecX = vecX0 + vecDelta;
		if (debugMode)
			msgif( debugMode, __FILE__, __LINE__, "Performing function evaluation..." );
			echo__vecX = vecX
			echo__vecFModel = vecFModel
			echo__omegaModel = omegaModel
		end
		vecF = funchF( vecX );
		datOut.fevalCount++;
		%
		% If feval failed, cut trust region size.
		if ( ~isrealarray(vecF,[sizeF,1]) )
			msgif( debugMode, __FILE__, __LINE__, "Function evaluation failed." );
			if (havePrev)
				msgif( __FILE__, __LINE__, "Function evaluation failed even though it succeeded for an earlier trial step!" );
				msgif( __FILE__, __LINE__, "This suggests the invalid space is not simply at distant values of x." );
			endif
			coeff_reduceTregOnFevalFail = 0.1;
			dTreg = coeff_reduceTregOnFevalFail * norm(vecDelta);
			continue;
		endif
		if (debugMode)
			echo__vecF = vecF
		end
		%
		%
		%
		% If vecF_trial is radically different from vecFModel_trial, bail.
		% We could check to see how close vecF_trial is to (1-b)*vecF0 + b*vecFModel_trial,
		%  for the value of b in [0,1] that minimizes the difference, but,
		%  assuming norm( vecFModel_trial - vecF0 ) <= norm( vecF0 ),
		%  and coeff_declareModelIsRadicallyWrong is > 2 ish(?), this isn't necessary.
		coeff_declareModelIsRadicallyWrong = 2.0;
		if ( norm( vecF - vecFModel ) >= coeff_declareModelIsRadicallyWrong * norm(vecF0) )
			msgif( debugMode, __FILE__, __LINE__, "Function was radically different at trial point." );
			if (havePrev)
				msgif( __FILE__, __LINE__, "Function was radically different even though it was reasonable for an earlier trial step!" );
				msgif( __FILE__, __LINE__, "This suggests the presence of a sepratrix." );
			endif
			coeff_reduceTregOnRadicallyWrong = 0.1;
			dTreg = coeff_reduceTregOnRadicallyWrong * norm(vecDelta);
			msgif( debugMode, __FILE__, __LINE__, sprintf( "Set dTreg to %g.", dTreg ) );
			continue;
		endif
		%
		omega = sumsq(vecF,1)/2.0
		if ( havePrev )
		if ( omega > omega_prev )
			msg( __FILE__, __LINE__, "omega > omega_prev; this suggests a sepratrix." );
		endif
		endif
		%
		%
		%
		% Let's update our model with this new info.
		msgif( debugMode, __FILE__, __LINE__, "Updating K!" );
		matK += (2.0*(omega-omegaModel)/sumsq(vecDelta)^2)*(vecDelta*(vecDelta'));
		%
		if ( ~haveBSF && omega < omega0 )
			haveBSF = true;
			vecX_bsf = vecX;
			vecF_bsf = vecF;
			omega_bsf = omega;
		elseif ( omega < omega_bsf )
			haveBSF = true; % Redundant.
			vecX_bsf = vecX;
			vecF_bsf = vecF;
			omega_bsf = omega;
		endif
		%
		% Under various circumstances, we could know that our trial point would be the next trial point,
		% so there's no need to calculate a new trial.
		% But, POITROME.
		havePrev = true;
		vecX_prev = vecX;
		vecF_prev = vecF;
		omega_prev = omega;
		continue;
	endwhile
	%
	if (nargout>=2)
		datOut.dTreg = dTreg;
		datOut.matK = matK;
	end
	if ( haveBSF )
		vecX_next = vecX_bsf
	else
		vecX_next = vecX0
		msg( __FILE__, __LINE__, "Failed." );
	end
	return;
endfunction
