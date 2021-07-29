function [ s, p, retCode, datOut ] = extFit__findStep( s0, p0, xVals, fVals, nPtWiseExt, wVals=[], prm=[] )
	commondefs;
	thisFile = "extFit__findStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	doChecks = mygetfield( prm, "doChecks", true );
	datOut = [];
	%
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
		numPts = size(xVals,2);
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isposintscalar(nPtWiseExt) );
		assert( 1 <= nPtWiseExt );
		assert( nPtWiseExt <= numPts );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	datOut.s0 = s0;
	datOut.p0 = p0;
	prm_calcAboutPt = mygetfield( prm, "prm_calcAboutPt", [] );
	[ rhoVals, bigF0, bigF1, omega0, vecG, matH ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, nPtWiseExt, wVals, prm_calcAboutPt );
	datOut.bigF0 = bigF0;
	datOut.bigF1 = bigF1;
	datOut.omega0 = omega0;
	datOut.vecG = vecG;
	datOut.matH = matH;
	assert( omega0 >= 0.0 );
	assert( matH(1,1) >= 0.0 );
	assert( matH(2,2) >= 0.0 );
	assert( matH(1,2).^2 <= matH(1,1)*matH(2,2) );
	assert( fleq(matH(1,2),matH(2,1)) );
	assert(~(  0.0 == matH(1,1)  &&  0.0 ~= vecG(1)  ));
	assert(~(  0.0 == matH(2,2)  &&  0.0 ~= vecG(2)  ));
	funch_omegaModel = @(vecD)( omega0 + vecD'*vecG + 0.5*vecD'*matH*vecD );
	%
	%
	%
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
	%
	hScale = sqrt(sum(sum(matH.^2)));
	if ( mygetfield( prm, "useLevMarq", false ) )
		epsRelLevMarq = mygetfield( prm, "epsRelLevMarq", eps );
		assert( isrealscalar(epsRelLevMarq) );
		assert( 0.0 < epsRelLevMarq );
		matD = sqrt( diag(diag(matH'*matH)) + hScale^2*epsRelLevMarq*eye(2,2) );
	else
		matD = hScale*eye(2,2);
	end
	datOut.matD = matD;
	rcondTol = mygetfield( prm, "rcondTol", eps^0.75 );
	assert( rcond(matD) > rcondTol );
	%
	if ( rcond(matH) > rcondTol )
		vecDeltaN = -matH\vecG;
	else
		vecDelta1 = -(matH+(1.0*epsMu*matD))\vecG;
		vecDelta2 = -(matH+(2.0*epsMu*matD))\vecG;
		vecDelta3 = -(matH+(3.0*epsMu*matD))\vecG;
		%
		%vecDeltaN0 = vecDelta1; % Cnst
		%vecDeltaN1 = (2.0*vecDelta1) - vecDelta2; % Linear.
		vecDeltaN2 = 3.0*(vecDelta1-vecDelta2) + vecDelta3; %Quadratic.
		vecDeltaN = vecDeltaN2;
	end
	omegaModelN = funch_omegaModel(vecDeltaN);
	if ( omegaModelN >= omega0*(1.0-eps^0.75) )
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Local model indicates we can't decrease from %g.", omega0 )  );
		s = s0;
		p = p0;
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	%
	%
	%
	% Set defaults for pMin and pMax.
	if ( 1 == nPtWiseExt )
		sMin = min(xVals) - 1.0*(max(xVals)-min(xVals));
		sMax = xVals(2)
	elseif ( numPts == nPtWiseExt )
		sMin = xVals(numPts-1);
		sMax = max(xVals) + 1.0*(max(xVals)-min(xVals));
	else
		if (doChecks)
		if ( isempty(mygetfield( prm, "sMin", [] ))  || isempty(mygetfield( prm, "sMax", [] )) )
		if ( 0.0 < ...
		    ( fVals(nPtWiseExt+1) - fVals(nPtWiseExt) ) ...
		  * ( fVals(nPtWiseExt) - fVals(nPtWiseExt-1) )  )
			msg_warn( verbLev, thisFile, __LINE__, "WARNING: nPtWiseExt is not an exterma." );
			msg_warn( verbLev, thisFile, __LINE__, "  But, default sMin/Max are set on this basis." );
		end
		end
		sMin = xVals(nPtWiseExt-1);
		sMax = xVals(nPtWiseExt+1);
	end
	pMin = min([ p0, 1.0 ]);
	pMax = max([ p0, 20.0 ]);
	%
	sMin = mygetfield( prm, "sMin", sMin );
	sMax = mygetfield( prm, "sMax", sMax );
	pMin = mygetfield( prm, "pMin", pMin );
	pMax = mygetfield( prm, "pMax", pMax );
	mu0 = mygetfield( prm, "mu0", (eps) );
	mu1 = mygetfield( prm, "mu1", 1.0/(eps) );
	muStep = mygetfield( prm, "muStep", 10.0 );
	epsMu = mygetfield( prm, "epsMu", eps^0.5 );
	sufficientDecreaseCoeff = mygetfield( prm, "sufficientDecreaseCoeff", 0.01 );
	if ( doChecks )
		assert( isrealscalar(sMin) );
		assert( isrealscalar(sMax) );
		assert( isrealscalar(pMin) );
		assert( isrealscalar(pMax) );
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( isrealscalar(epsMu) );
		assert( sMin <= s0 );
		assert( s0 <= sMax );
		assert( 0.0 < pMin );
		assert( pMin <= p0 );
		assert( p0 <= pMax );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
		assert( 0.0 < epsMu );
	end
	%
	%
	%
	% DO WORK.
	%
	%
	%
	iterCount = 0;
	while (1)
		iterCount++;
		switch (iterCount)
		case 1
			mu_trial = 0.0;
			vecDelta_trial = vecDeltaN;
		case 2
			mu_trial = mu0;
			vecDelta_trial = -(matH + (mu_trial*matD))\vecG;
		otherwise
			mu_trial *= muStep;
			vecDelta_trial = -(matH + (mu_trial*matD))\vecG;
		end
		if ( mu_trial > mu1 )
			break;
		end
		%
		s_trial = s0 + vecDelta_trial(1);
		p_trial = p0 + vecDelta_trial(2);
		omegaModel_trial = funch_omegaModel( vecDelta_trial );
		if ( omegaModel_trial >= omega0*(1.0-eps^0.75) )
			break;
		end
		if ( s_trial < sMin || s_trial > sMax ...
		  || p_trial < pMin || p_trial > pMax )
			continue;
		end
		[ rhoVals_trial, bigF0_trial, bigF1_trial, omega_trial ] = extFit__calcAtPt( ...
		  s_trial, p_trial, xVals, fVals, nPtWiseExt, wVals, prm_calcAboutPt );
		if ( omega_trial >= omega0 )
			continue;
		end
		if ( omega_trial > omega0 - sufficientDecreaseCoeff*abs(omegaModel_trial-omega0) )
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