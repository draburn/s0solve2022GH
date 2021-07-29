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
	prm_calcAboutPt = mygetfield( prm, "prm_calcAboutPt", [] );
	[ rhoVals, bigF0, bigF1, omega0, vecG, matH ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, nPtWiseExt, wVals, prm_calcAboutPt );
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
	funch_deltaOmegaModel = @(vecD)( vecD'*vecG + 0.5*vecD'*matH*vecD );
	%
	%
	%
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
	reguCoeff = mygetfield( prm, "reguCoeff", eps075 );
	rcondTol = mygetfield( prm, "rcondTol", eps075 );
	epsMu = mygetfield( prm, "epsMu", eps075 );
	sMin = mygetfield( prm, "sMin", sMin );
	sMax = mygetfield( prm, "sMax", sMax );
	pMin = mygetfield( prm, "pMin", pMin );
	pMax = mygetfield( prm, "pMax", pMax );
	mu0 = mygetfield( prm, "mu0", eps050 );
	mu1 = mygetfield( prm, "mu1", 1.0/eps075 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	sufficientDecreaseCoeff = mygetfield( prm, "sufficientDecreaseCoeff", 0.01 );
	if ( doChecks )
		assert( isrealscalar(reguCoeff) );
		assert( isrealscalar(rcondTol) );
		assert( isrealscalar(epsMu) );
		assert( isrealscalar(sMin) );
		assert( isrealscalar(sMax) );
		assert( isrealscalar(pMin) );
		assert( isrealscalar(pMax) );
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( 0.0 < reguCoeff );
		assert( 0.0 < rcondTol );
		assert( 0.0 < epsMu );
		assert( sMin <= s0 );
		assert( s0 <= sMax );
		assert( 0.0 < pMin );
		assert( pMin <= p0 );
		assert( p0 <= pMax );
		assert( epsMu < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
	end
	%
	%
	%
	% DO WORK.
	rvecHTHDiag = sum(matH.^2,1);
	hScale = sqrt(sum(rvecHTHDiag)/2.0);
	matD = sqrt( diag(rvecHTHDiag) + (hScale^2)*reguCoeff*eye(2,2) );
	%
	if ( rcond(matH) > rcondTol )
		vecDeltaN = -matH\vecG;
	else
		vecDelta1 = -(matH+(1.0*epsMu*matD))\vecG;
		vecDelta2 = -(matH+(2.0*epsMu*matD))\vecG;
		vecDelta3 = -(matH+(3.0*epsMu*matD))\vecG;
		%
		vecDeltaN0 = vecDelta1; % Cnst
		vecDeltaN1 = (2.0*vecDelta1) - vecDelta2; % Linear.
		vecDeltaN2 = 3.0*(vecDelta1-vecDelta2) + vecDelta3; %Quadratic.
		vecDeltaN = vecDeltaN2;
	end
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
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Only excessive backtracking (%g > %g) might reduce residual (%g).", ...
			   mu_trial, mu1, omega0 ) );
			break;
		end
		%
		s_trial = s0 + vecDelta_trial(1);
		p_trial = p0 + vecDelta_trial(2);
		deltaOmegaModel_trial = funch_deltaOmegaModel( vecDelta_trial );
		assert( 0.0 > deltaOmegaModel_trial );
		if ( abs(deltaOmegaModel_trial) < omega0*eps075 )
			% No hope.
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Model suggests only insignificant decrease (%g) from omega (%g) is possible.", ...
			   deltaOmegaModel_trial, omega0 ) );
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
