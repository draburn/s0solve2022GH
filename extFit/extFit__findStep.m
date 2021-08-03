function [ s, p, retCode, datOut ] = extFit__findStep( ...
  s0, p0, xVals, fVals, wVals, prm=[] )
	commondefs;
	thisFile = "extFit__findStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	doChecks = mygetfield( prm, "doChecks", true );
	datOut = [];
	%
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
	end
	%
	prm_calcAboutPt = mygetfield( prm, "prm_calcAboutPt", [] );
	[ rhoVals, bigF0, bigF1, omega0, vecG, matH ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, wVals, prm_calcAboutPt );
	assert( omega0 >= 0.0 );
	assert( matH(1,1) >= 0.0 );
	assert( matH(2,2) >= 0.0 );
	assert( matH(1,2).^2 <= matH(1,1)*matH(2,2) );
	assert( fleq(matH(1,2),matH(2,1)) );
	assert(~(  0.0 == matH(1,1)  &&  0.0 ~= vecG(1)  ));
	assert(~(  0.0 == matH(2,2)  &&  0.0 ~= vecG(2)  ));
	if ( omega0 == 0.0 )
		msg_notify( verbLev, thisFile, __LINE__, "Initial omega is already zero." );
		s = s0;
		p = p0;
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	if ( norm(vecG) == 0.0 )
		msg_notify( verbLev, thisFile, __LINE__, "Initial gradient is zero." );
		s = s0;
		p = p0;
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	if (  0.0 == matH(1,1)  ||  0.0 == matH(2,2)  )
		msg_notify( verbLev, thisFile, __LINE__, ...
		  "Hessian has a zero on-diagonal; this case is not supported." );
		s = s0;
		p = p0;
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
		if (~isempty(pMin))
			assert( isrealscalar(pMin) );
			assert( 0.0 <= pMin );
			assert( pMin <= p0 );
		end
		if (~isempty(pMax))
			assert( isrealscalar(pMax) );
			assert( p0 <= pMax );
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
			  "Only excessive backtracking (%g > %g) might reduce residual (%g).", ...
			   mu_trial, mu1, omega0 ) );
			break;
		end
		%
		vecDelta_trial = -(matH + (mu_trial*matD))\vecG;
		s_trial = s0 + vecDelta_trial(1);
		p_trial = p0 + vecDelta_trial(2);
		if ( !isempty(sMin) )
		if ( s_trial < sMin )
			continue;
		end
		end
		if ( !isempty(sMax) )
		if ( s_trial > sMax )
			continue;
		end
		end
		if ( !isempty(pMin) )
		if ( p_trial < pMin )
			continue;
		end
		end
		if ( !isempty(pMax) )
		if ( p_trial > pMax )
			continue;
		end
		end
		%
		deltaOmegaModel_trial = funch_deltaOmegaModel( vecDelta_trial );
		assert( 0.0 > deltaOmegaModel_trial );
		if ( abs(deltaOmegaModel_trial) < omega0*eps075 )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Model suggests only insignificant decrease (%g) from omega (%g) is possible.", ...
			   deltaOmegaModel_trial, omega0 ) );
			break;
		end
		%
		[ rhoVals_trial, bigF0_trial, bigF1_trial, omega_trial ] = extFit__calcAtPt( ...
		  s_trial, p_trial, xVals, fVals, wVals, prm_calcAboutPt );
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
	s = s0;
	p = p0;
	retCode = RETCODE__IMPOSED_STOP;
return;
end
