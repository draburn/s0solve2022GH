function [ s, p, retCode, datOut ] = extFit__findStep( s0, p0, xVals, fVals, nPtWiseExt, wVals=[], prm=[] )
	commondefs;
	thisFile = "extFit__findStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
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
	[ rhoVals, bigF0, bigF1, omega, vecG, matH ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, nPtWiseExt, wVals, prm_calcAboutPt );
	datOut.bigF0 = bigF0;
	datOut.bigF1 = bigF1;
	datOut.omega = omega;
	datOut.vecG = vecG;
	datOut.matH = matH;
	if ( omega == 0.0 )
		msg_notify( verbLev, thisFile, __LINE__, "Input omega is already zero." );
		retCode = RETCODE__IMPOSED_STOP;
		return;
	end
	assert( matH(1,1) >= 0.0 );
	assert( matH(2,2) >= 0.0 );
	%
	if ( mygetfield( prm, "useLevMarq", false ) )
		matD = diag(diag(matH));
	else
		matD = eye(2,2);
	end
	%
	%
	%
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
	pMin = 1.0;
	pMax = 20.0;
	%
	sMin = mygetfield( prm, "sMin", sMin );
	sMax = mygetfield( prm, "sMax", sMax );
	pMin = mygetfield( prm, "pMin", pMin );
	pMax = mygetfield( prm, "pMax", pMax );
	mu0 = mygetfield( prm, "mu0", sqrt(eps) );
	mu1 = mygetfield( prm, "mu1", 1.0/sqrt(eps) );
	muStep = mygetfield( prm, "muStep", 10.0 );
	if ( doChecks )
		assert( isrealscalar(sMin) );
		assert( isrealscalar(sMax) );
		assert( isrealscalar(pMin) );
		assert( isrealscalar(pMax) );
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( sMin <= s0 );
		assert( s0 <= sMax );
		assert( 0.0 < pMin );
		assert( pMin <= p0 );
		assert( p0 <= pMax );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
	end
	%
	iterNum = 0;
	while (1)
		switch (iterNum)
		case 0
			mu = 0.0;
		case 1
			mu = mu0;
		otherwise
			mu *= muStep;
		end
		iterNum++;
		%echo__mu = mu
		if ( mu > mu1 )
			break;
		end
		%
		matA = matH + (mu*matD);
		vecDelta_trial = -matA\vecG;
		s_trial = s0 + vecDelta_trial(1);
		p_trial = p0 + vecDelta_trial(2);
		if ( p_trial < 1.0 )
			continue;
		elseif ( s_trial < sMin )
			continue;
		elseif ( s_trial > sMax )
			continue;
		end
		%
		datOut.s_trial = s_trial;
		datOut.p_trial = p_trial;
		[ rhoVals_trial, bigF0_trial, bigF1_trial, omega_trial ] = extFit__calcAtPt( ...
		  s_trial, p_trial, xVals, fVals, nPtWiseExt, wVals, prm_calcAboutPt );
		datOut.rhoVals_trial = rhoVals_trial;
		datOut.bigF0_trial = bigF0_trial;
		datOut.bigF1_trial = bigF1_trial;
		datOut.omega_trial = omega_trial;
		if ( omega_trial >= omega )
			continue;
		end
		%
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Decreased omega from %g to %g.", omega, omega_trial ) );
		s = s_trial;
		p = p_trial;
		retCode = RETCODE__SUCCESS;
		return;
	end
	msg_main( verbLev, thisFile, __LINE__, sprintf( ...
	  "Failed to decrease omega from %g.", omega )  );
	s = s0;
	p = p0;
	retCode = RETCODE__IMPOSED_STOP;
return;
end
