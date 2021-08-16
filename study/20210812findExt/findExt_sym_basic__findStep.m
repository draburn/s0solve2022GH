function [ xExt_new, p_new, retCode, datOut ] = findExt_sym_basic__findStep( ...
  xExt, p, xVals, fVals, ...
  xExtMin, xExtMax, pMin, pMax, dVals, ...
  prm );
	%
	% Init
	commondefs;
	thisFile = "findExt_sym_basic__findStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	%
	%
	numPts = size(xVals,2);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(xExt) );
		assert( isrealscalar(p) );
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealorinfscalar(xExtMin) );
		assert( isrealorinfscalar(xExtMax) );
		assert( xExtMin <= xExt );
		assert( xExt <= xExtMax );
		assert( isrealorinfscalar(pMin) );
		assert( isrealorinfscalar(pMax) );
		assert( pMin <= p );
		assert( p <= pMax );
		assert( isrealarray(dVals,[1,numPts]) );
		assert( sum(dVals) > 0.0 );
	end
	epsXExt = mygetfield( prm, "epsXExt", eps025*(max(xVals)-min(xVals)) );
	epsP = mygetfield( prm, "epsP", eps025 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(epsXExt) );
		assert( isrealscalar(epsP) );
		assert( epsXExt > 0.0 );
		assert( epsP > 0.0 );
	end
	%
	xExt_new = [];
	p_new = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	% 
	% Calculate Hessian.
	[ flag_00, rhoVals, omega, fExt, f1 ] = findExt_sym__calcRes( xExt, p, xVals, fVals, dVals );
	[ flag_p0, rhoVals_p0, omega_p0 ] = findExt_sym__calcRes( xExt+epsXExt, p, xVals, fVals, dVals );
	[ flag_m0, rhoVals_m0, omega_m0 ] = findExt_sym__calcRes( xExt-epsXExt, p, xVals, fVals, dVals );
	[ flag_0p, rhoVals_0p, omega_0p ] = findExt_sym__calcRes( xExt,    p+epsP, xVals, fVals, dVals );
	[ flag_0m, rhoVals_0m, omega_0m ] = findExt_sym__calcRes( xExt,    p-epsP, xVals, fVals, dVals );
	if ( flag_00 || flag_p0 || flag_m0 || flag_0m || flag_0p )
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	%
	rhoD0Vals = dVals .* rhoVals;
	rhoDXVals = dVals .* ( rhoVals_p0 - rhoVals_m0 ) / ( 2.0*epsXExt );
	rhoDPVals = dVals .* ( rhoVals_0p - rhoVals_0m ) / ( 2.0*epsP );
	%
	vecG = [ ...
	  sum( rhoD0Vals.*rhoDXVals ); ...
	  sum( rhoD0Vals.*rhoDPVals) ];
	matH = [ ...
	  sum(rhoDXVals.*rhoDXVals), sum(rhoDXVals.*rhoDPVals); ...
	  sum(rhoDXVals.*rhoDPVals), sum(rhoDPVals.*rhoDPVals) ];
	if ( verbLev >= VERBLEV__COPIOUS )
		msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
		echo__omega = omega
		echo__fExt = fExt
		echo__f1 = f1
		echo__matH = matH
		echo__vecG = vecG
		%echo__vecG1 = vecG(1)
		%echo__vecG2 = vecG(2)
		msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( omega >= 0.0 );
		assert( matH(1,1) >= 0.0 );
		assert( matH(2,2) >= 0.0 );
		assert( matH(1,2).^2 <= matH(1,1)*matH(2,2) );
		assert( fleq(matH(1,2),matH(2,1)) );
		assert(~(  0.0 == matH(1,1)  &&  0.0 ~= vecG(1)  ));
		assert(~(  0.0 == matH(2,2)  &&  0.0 ~= vecG(2)  ));
	end
	if ( 0.0 == omega )
		msg_notify( verbLev, thisFile, __LINE__, "Initial omega is already zero." );
		xExt_new = xExt;
		p_new = p;
		retCode = RETCODE__SUCCESS;
		return;
	end
	if ( 0.0 == norm(vecG) )
		msg_notify( verbLev, thisFile, __LINE__, "Initial gradient is zero." );
		xExt_new = xExt;
		p_new = p;
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		return;
	end
	if ( matH(1,1) <= eps075*sum(sum(abs(matH))) )
		msg_warn( verbLev, thisFile, __LINE__, "Hessian (1,1) is near-singular; using a hack." );
		matH(1,1) = norm(vecG);
	end
	if ( matH(2,2) <= eps075*sum(sum(abs(matH))) )
		msg_warn( verbLev, thisFile, __LINE__, "Hessian (2,2) is near-singular; using a hack." );
		matH(2,2) = norm(vecG);
	end
	matD = diag(abs(diag(matH)));
	%funch_deltaOmegaModel = @(vecDelta)( vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta );
	%
	%
	%
	mu0 = mygetfield( prm, "mu0", 1.0E-8 );
	mu1 = mygetfield( prm, "mu1", 1.0E8 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
	end
	%
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
			msg_main( verbLev, thisFile, __LINE__, sprintf( "Reached maximum backtracking (%g).", mu1 ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		%
		vecDelta_trial = -(matH + (mu_trial*matD))\vecG;
		xExt_trial = xExt + vecDelta_trial(1);
		p_trial = p + vecDelta_trial(2);
		if ( (xExt_trial < xExtMin) || (xExt_trial > xExtMax) || (p_trial < pMin) || (p_trial > pMax) )
			continue;
		end
		[ flag_trial, rhoVals_trial, omega_trial ] = findExt_sym__calcRes( xExt_trial, p_trial, xVals, fVals, dVals );
		if ( omega_trial < omega )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reduced omega from %g to %g with mu %g.", omega, omega_trial, mu_trial ) );
			xExt_new = xExt_trial;
			p_new = p_trial;
			retCode = RETCODE__SUCCESS;
			return;
		end
	end
end
