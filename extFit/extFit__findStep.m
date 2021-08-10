function [ s, p, retCode, datOut ] = extFit__findStep( ...
  s0, p0, xVals, fVals, wVals, prm=[] )
	commondefs;
	thisFile = "extFit__findStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	%
	% Default return values.
	s = s0;
	p = p0;
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
	end
	%
	%
	prm_calcAboutPt = mygetfield( prm, "prm_calcAboutPt", [] );
	[ rhoVals, bigF0, bigF1, omega0, vecG, matH, retCode ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, wVals, prm_calcAboutPt );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_warn( verbLev, thisFile, __LINE__, "Calculation about initial point failed." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	retCode = RETCODE__NOT_SET;
	if ( verbLev >= VERBLEV__COPIOUS )
		msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
		echo__bigF0 = bigF0
		echo__bigF1 = bigF1
		echo__omega0 = omega0
		echo__matH = matH
		echo__vecG = vecG
		msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	end
	assert( omega0 >= 0.0 );
	assert( matH(1,1) >= 0.0 );
	assert( matH(2,2) >= 0.0 );
	assert( matH(1,2).^2 <= matH(1,1)*matH(2,2) );
	assert( fleq(matH(1,2),matH(2,1)) );
	assert(~(  0.0 == matH(1,1)  &&  0.0 ~= vecG(1)  ));
	assert(~(  0.0 == matH(2,2)  &&  0.0 ~= vecG(2)  ));
	if ( omega0 == 0.0 )
		msg_notify( verbLev, thisFile, __LINE__, "Initial omega is already zero." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	if ( norm(vecG) == 0.0 )
		msg_notify( verbLev, thisFile, __LINE__, "Initial gradient is zero." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	if (  0.0 == matH(1,1)  ||  0.0 == matH(2,2)  )
		msg_notify( verbLev, thisFile, __LINE__, ...
		  "Hessian has a zero on-diagonal; this case is not supported." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	rvecHTHDiag = sum(matH.^2,1);
	hScale = sqrt(sum(rvecHTHDiag)/2.0);
	funch_deltaOmegaModel = @(vecD)( vecD'*vecG + 0.5*vecD'*matH*vecD );
	if ( mygetfield(prm,"useLevMarq",true) )
		matD = diag(abs(diag(matH)));
	else
		matD = hScale*eye(2,2);
	end
	%
	%
	%
	sMin = mygetfield( prm, "sMin", [] );
	sMax = mygetfield( prm, "sMax", [] );
	pMin = mygetfield( prm, "pMin", 0.0 );
	pMax = mygetfield( prm, "pMax", [] );
	mu0 = mygetfield( prm, "mu0", eps050 );
	mu1 = mygetfield( prm, "mu1", 1.0/eps075 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	sufficientDecreaseCoeff = mygetfield( prm, "sufficientDecreaseCoeff", 0.01 );
	if ( doChecks )
		if (~isempty(sMin))
			assert( isrealscalar(sMin) );
			assert( sMin <= s0 );
		end
		if (~isempty(sMax))
			assert( isrealscalar(sMax) );
			assert( s0 <= sMax );
		end
		if ( ~isempty(sMin) && ~isempty(sMax) )
			assert( sMin < sMax );
		end
		if (~isempty(pMin))
			assert( isrealscalar(pMin) );
			assert( 0.0 <= pMin );
			assert( pMin <= p0 );
		end
		if (~isempty(pMax))
			assert( isrealscalar(pMax) );
			assert( p0 <= pMax );
		end
		if ( ~isempty(pMin) && ~isempty(pMax) )
			assert( pMin < pMax );
		end
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
		assert( isrealscalar(sufficientDecreaseCoeff) );
		assert( 0.0 <= sufficientDecreaseCoeff );
		assert( sufficientDecreaseCoeff < 1.0 );
	end
	%
	%
	%
	% DO WORK.
	haveTriedLoLo = false;
	haveTriedLoHi = false;
	haveTriedHiLo = false;
	haveTriedHiHi = false;
	iterCount = 0;
	while (1)
		iterCount++;
		switch (iterCount)
		case 1
			mu_trial = 0.0;
		case 2
			mu_trial = mu0;
		otherwise
			mu_trial *= muStep;
		end
		if ( mu_trial > mu1 )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached maximum backtracking (%g).", mu1 ) );
			retCode = RETCODE__IMPOSED_STOP; % Unless changed below.
			break;
		end
		%
		vecDelta_trial = -(matH + (mu_trial*matD))\vecG; % Unless modified below.
		%
		%
		sRegion = 0;
		pRegion = 0;
		if ( ~isempty(sMin) )
		if ( s + vecDelta_trial(1) < sMin )
			sRegion = -1;
		end
		end
		if ( ~isempty(sMax) )
		if ( s + vecDelta_trial(1) > sMax )
			sRegion = +1;
		end
		end
		if ( ~isempty(pMin) )
		if ( p + vecDelta_trial(2) < pMin )
			pRegion = -1;
		end
		end
		if ( ~isempty(pMax) )
		if ( p - vecDelta_trial(2) > pMax )
			pRegion = +1;
		end
		end
		%
		if ( sRegion > 0 )
			vecDelta_trial(1) = (sMax - s)*(1.0-eps075);
			if ( pRegion > 0 )
				vecDelta_trial(2) = (pMax - p)*(1.0-eps075);
				if ( haveTriedHiHi )
					continue;
				end
				haveTriedHiHi = true;
			elseif ( pRegion < 0 )
				vecDelta_trial(2) = (pMin - p)*(1.0-eps075);
				if ( haveTriedHiLo )
					continue;
				end
				haveTriedHiLo = true;
			else
				vecDelta_trial(2) = -( vecG(2) + vecDelta_trial(1)*(matH(1,2)+mu_trial*matD(1,2)) ) ...
				  / ( matH(2,2) + mu_trial*matD(2,2) );
			end
		elseif ( sRegion < 0 )
			vecDelta_trial(1) = (sMin - s)*(1.0-eps075);
			if ( pRegion > 0 )
				vecDelta_trial(2) = (pMax - p)*(1.0-eps075);
				if ( haveTriedLoHi )
					continue;
				end
				haveTriedLoHi = true;
			elseif ( pRegion < 0 )
				vecDelta_trial(2) = (pMin - p)*(1.0-eps075);
				if ( haveTriedLoLo )
					continue;
				end
				haveTriedLoLo = true;
			else
				vecDelta_trial(2) = -( vecG(2) + vecDelta_trial(1)*(matH(1,2)+mu_trial*matD(1,2)) ) ...
				  / ( matH(2,2) + mu_trial*matD(2,2) );
			end
		else
			if ( pRegion > 0 )
				vecDelta_trial(2) = (pMax - p)*(1.0-eps075);
				vecDelta_trial(1) = -( vecG(1) + vecDelta_trial(2)*(matH(1,2)+mu_trial*matD(1,2)) ) ...
				  / ( matH(1,1) + mu_trial*matD(1,1) );
			elseif ( pRegion < 0 )
				vecDelta_trial(2) = (pMin - p)*(1.0-eps075);
				vecDelta_trial(1) = -( vecG(1) + vecDelta_trial(2)*(matH(1,2)+mu_trial*matD(1,2)) ) ...
				  / ( matH(1,1) + mu_trial*matD(1,1) );
			else
				% Nothing to do!
			end
		end
		s_trial = s0 + vecDelta_trial(1);
		p_trial = p0 + vecDelta_trial(2);
		%
		if ( ~isempty(sMin) )
			assert( s_trial >= sMin );
		end
		if ( ~isempty(sMax) )
			assert( s_trial <= sMax );
		end
		if ( ~isempty(pMin) )
			assert( p_trial >= pMin);
		end
		if ( ~isempty(pMax) )
			assert( p_trial <= pMax );
		end
		%
		deltaOmegaModel_trial = funch_deltaOmegaModel( vecDelta_trial );
		assert( 0.0 > deltaOmegaModel_trial );
		if ( abs(deltaOmegaModel_trial) < omega0*eps075 )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Hit local minimum (%g from %g).", -deltaOmegaModel_trial, omega0 ) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN; % Unless changed below.
			break;
		end
		%
		[ rhoVals_trial, bigF0_trial, bigF1_trial, omega_trial, retCode ] = extFit__calcAtPt( ...
		  s_trial, p_trial, xVals, fVals, wVals, prm_calcAboutPt );
		if ( RETCODE__SUCCESS ~= retCode )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( "Calculation failed for mu = %g.", mu_trial ) );
			if ( verbLev >= VERBLEV__COPIOUS )
				msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
				echo__mu_trial = mu_trial
				echo__s_trial = s_trial
				echo__p_trial = p_trial
				msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
			end
			retCode = RETCODE__NOT_SET;
			continue;
		end
		retCode = RETCODE__NOT_SET;
		assert( 0.0 <= omega_trial );
		if ( omega_trial >= omega0 )
			continue;
		end
		if ( omega_trial > omega0 - sufficientDecreaseCoeff*abs(deltaOmegaModel_trial) )
			continue;
		end
		%
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Decreased omega from %g to %g.", omega0, omega_trial ) );
		s = s_trial;
		p = p_trial;
		retCode = RETCODE__SUCCESS;
		return;
	end
	%
	msg_main( verbLev, thisFile, __LINE__, sprintf( ...
	  "Failed to decrease omega from %g.", omega0 )  );
return;
end
