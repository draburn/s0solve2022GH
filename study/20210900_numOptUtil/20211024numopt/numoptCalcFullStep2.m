function [ vecDeltaFinal, retCode, datOut ] = numoptCalcFullStep( omega0, vecG, matH, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "numoptCalcFullStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	assert( isrealscalar(verbLev) );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	assert( isrealscalar(omega0) );
	assert( omega0 >= 0.0 );
	probSize = size( vecG,1 );
	assert( 1 <= probSize );
	assert( isrealarray(vecG,[probSize,1]) );
	assert( isrealarray(matH,[probSize,probSize]) );
	assert( issymmetric(matH) );
	%
	msg( thisFile, __LINE__, "TODO: Add LUT." );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK
	%
	%
	if ( 0.0 == omega0 )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: omega0 = 0.0." );
		muFinal = +Inf;
		vecDeltaFinal = zeros(sizeX,1);
		omegaFinal = omega0;
		return;
	end
	%
	gSq = vecG'*vecG;
	if ( 0.0 == gSq )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: vecG = 0.0." );
		muFinal = +Inf;
		vecDeltaFinal = zeros(sizeX,1);
		omegaFinal = omega0;
		return;
	end
	%
	hAbsMax = max(max(abs(matH)));
	if ( 0.0 == hAbsMax )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: matH = 0.0." );
		muFinal = gSq / omega0;
		vecDeltaFinal = -vecG * omega0 / gSq;
		omegaFinal = 0.0;
		return;
	end
	%
	%
	%
	funchOmega = @(vecDummy)( omega0 + (vecG'*vecDummy) + (0.5*(vecDummy' * matH * vecDummy)) );
	matI = eye(probSize,probSize);
	hFrobNorm = sqrt(sum(sum(matH.^2)));
	gthg = vecG'*matH*vecG;
	%
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag ) % Nominally pos-def.
	if ( min(diag(matR)) > eps050*max(abs(diag(matR))) ) % Pos-def within tolerance.
		muFinal = 0.0;
		vecDeltaFinal = -( matR \ (matR'\vecG) );
		omegaFinal = funchOmega( vecDeltaFinal );
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Apparent Result: Posdef (%e, %e).", muFinal, omegaFinal ) );
		return;
	end
	end
	clear matR;
	clear cholFlag;
	%
	%
	%
	% Try positive-semi-definite case.
	% Omega may go to negative infinity, but in case it doesn't,
	%  we want to extrapolate to get the pseudo-inverse point.
	muExtrap = eps075*hAbsMax;
	[ matR1, cholFlag ] = chol( matH + (muExtrap * matI) );
	if ( 0 == cholFlag )
	if ( min(diag(matR1)) > eps050*max(abs(diag(matR1))) )
		vecDelta1 = -( matR1 \ ( matR1'\vecG ) );
		omega1 = funchOmega( vecDelta1 );
		if ( omega1 < -eps075*omega0 )
			muFinal = muExtrap;
			vecDeltaFinal = vecDelta1;
			omegaFinal = funchOmega( vecDeltaFinal );
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Apparent Result: Divergent PSD (%e, %e).", muFinal, omegaFinal ) );
			return;
		end
		%
		matR2 = chol( matH + (2.0 * muExtrap * matI) );
		assert( 0 == cholFlag );
		vecDelta2 = -( matR2 \ ( matR2'\vecG ) );
		%
		muFinal = 0.0;
		vecDeltaFinal = (2.0*vecDelta1) - vecDelta2;
		omegaFinal = funchOmega( vecDeltaFinal );
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Apparent Result: Finite PSD (%e, %e).", muFinal, omegaFinal ) );
		return;
	end
	end
	%
	%
	% We've gotten here, so matH has at least one negative eigenvalue.
	%
	% First, do a simple method that would work if we don't need to get close to
	% the singularity (as measured by mu).
	%
	omegaThresh = mygetfield( prm, "omegaThresh", 0.0 );
	assert( isrealscalar(omegaThresh) );
	if ( omegaThresh >= omega0 )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "WARNING: omegaThresh >= omega0 (%e, %e).", omegaThresh, omega0 ) );
		muFinal = +Inf;
		vecDeltaFinal = zeros(sizeX,1);
		omegaFinal = omega0;
		return;
	end
	%
	deltaNormThresh = mygetfield( prm, "deltaNormThresh", +Inf );
	%
	mu1 = hFrobNorm*(1.0+eps025);
	%
	[ matR1, cholFlag ] = chol( matH + mu1*matI );
	assert( 0 == cholFlag );
	vecDelta1 = -(matR1 \ (matR1'\vecG) );
	omega1 = funchOmega(vecDelta1);
	if ( omega1 <= omegaThresh )
		msg_flagged( verbLev, thisFile, __LINE__, "THIS SITUATION IS VERY UNLIKELY! (SEE CODE.)" );
		muFinal = mu1;
		vecDeltaFinal = vecDelta1;
		omegaFinal = omega1;
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Apparent Result: RARE, ABSURDLY divergent HN (%e, %e).", muFinal, omegaFinal ) );
		return;
	end
	% Model: kappa = A / ( mu + B )^2, iota = C / ( mu + B )^D.
	kappa1 = 0.5*(vecDelta'*vecDelta);
	iota1 = omega0 - omega1;
	vecDeltaPrime1 = -(matR1 \ (matR1'\vecDelta1) );
	kappaPrime1 = vecDelta'*vecDeltaPrime;
	iotaPrime1 = -( (vecG'*vecDeltaPrime) + (vecDelta'*matH*vecDeltaPrime) );
	cnstX1 = -2.0*kappa1/kappaPrime1;
	cnstA = kappa1 * (cnstX1^2);
	cnstB = cnstX1 - mu1;
	cnstD = -cnstX1*iotaPrime1/iota1;
	cnstC = iota1*(cnstX1^cnstD);
	%
	iotaTarget = 1.2*(omega0 - omegaThresh); % Coeff is free param.
	mu2 = pow( cnstC / iotaTarget, 1.0 / cnstD ) - cnstB;
	if (isrealscalar(deltaNormThresh))
		kappaThresh = deltaNormThresh^2;
		kappaTarget = 1.2 * kappaThresh; % Coeff is free param.
		mu2Kappa = sqrt( cnstA / kappaTarget ) - cnstB;
		mu2 = max([ mu2, mu2Kappa ]); % Larger mu = smaller delta.
	end
	% Note that mu2 might even be negative.
	if ( mu2 > muExtrap )
		[ matR2, cholFlag ] = chol( matH + mu2*matI );
		if ( 0 == cholFlag )
			vecDelta2 = -(matR2 \ (matR2'\vecG) );
			omega2 = funchOmega(vecDelta2);
			if ( omega2 <= omegaThresh )
				muFinal = mu2;
				vecDeltaFinal = vecDelta2;
				omegaFinal = omega2;
				msg_main( verbLev, thisFile, __LINE__, sprintf( ...
				  "Apparent Result: Easy divergent HN (%e, %e).", muFinal, omegaFinal ) );
				return;
			end
			%
			% If nothing else, "2" serves as a new upper-bound.
			mu1 = mu2;
			matR1 = matR2;
			vecDelta1 = vecDelta2;
			omega1 = omega2;
		end
	end
	%
	%
	% That simple method didn't work.
	% We could try an iterative modification of that simple method, but,
	% analysis has suggested that, when the above simple method doesn't work,
	% we probably need to find a good approximation to the most negative eigenvalue anyway.
	% There seem to be three meaningful approaches now:
	%  1. Do a direct eigenvale calculation;
	%  2. Use a crude non-infinite mu (perhaps mu1 or mu2 from above); or,
	%  3. Start a linear solver such as GMRes or CG.
	%   At least for GMRes, a breakdown in the conventional algorithm would
	%   correspond to information about a direction for which matH is negative,
	%   and would be a good thing rather than a bad thing.
	%   Also, the solution space would be different from the Levenberg(-Marquardt) curve.
	%
	doFullEig = mygetfield( prm, "doFullEig", true );
	if (doFullEig)
		% Octave-wise, in cases, eig(matH) works better than eigs(matH,1,'sa').
		[ matPsi, matLambda ] = eig(matH);
		muCrit = -min(diag(matLambda));
		assert( muCrit >= 0.0 );
		%
		mu1 = (1.0+eps050)*muCrit + eps075*hAbsMax;
		[ matR1, cholFlag ] = chol( matH + mu1*matI );
		assert( 0 == cholFlag );
		vecDelta1 = -(matR1 \ (matR1'\vecG) );
		omega1 = funchOmega(vecDelta1);
		%
		if ( omega1 < 0.0 )
			muFinal = mu1;
			vecDeltaFinal = vecDelta1;
			omegaFinal = omega1;
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Apparent Result: Non-easy divergent HN (%e, %e).", muFinal, omegaFinal ) );
			return;
		end
		%
		mu2 = mu1 + eps075*hAbsMax;
		[ matR2, cholFlag ] = chol( matH + mu2*matI );
		assert( 0 == cholFlag );
		vecDelta2 = -(matR2 \ (matR2'\vecG) );
		%
		muFinal = muCrit;
		vecDeltaFinal = (  (vecDelta1 * ( mu2 - muCrit )) - (vecDelta2 * ( mu1 - muCrit )) ) / ( mu2 - mu1 );
		omegaFinal = funchOmega( vecDeltaFinal );
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Apparent Result: Finite HN (which is never easy) (%e, %e).", muFinal, omegaFinal ) );
		return;
	end
	%
	% Here, for simplicity, let's use method 2.
	muFinal = mu1;
	vecDeltaFinal = vecDelta1;
	omegaFinal = omega1;
	%
	%
return;
end
