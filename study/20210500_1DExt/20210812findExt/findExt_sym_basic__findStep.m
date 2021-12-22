function [ xExt_new, p_new, retCode, datOut ] = findExt_sym_basic__findStep( ...
  xExt, p, xVals, fVals, ...
  xExtMin, xExtMax, pMin, pMax, dVals, ...
  prm );
	%
	% Init
	commondefs;
	thisFile = "findExt_sym_basic__findStep";
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
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
		msg_warn( verbLev, thisFile, __LINE__, "Initial guess produced invalid results." );
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
		echo__xExt = xExt
		echo__p = p
		echo__omega = omega
		echo__fExt = fExt
		echo__f1 = f1
		echo__matH = matH
		echo__vecG = vecG
		msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	end
	
	msg( thisFile, __LINE__, "TODO: Fix regularization." );
	%omegaScale = omega;
	%xScale = max(xVals)-min(xVals);
	%pScale = 1.0;
	%vecG += eps025*omegaScale*[ xExt/(xScale^2); (p-2.0)/(pScale^2) ];
	%matH(1,1) += eps025*omegaScale/(xScale^2);
	%matH(2,2) += eps025*omegaScale/(pScale^2);
	
	
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
	
	
	
	mu0 = mygetfield( prm, "mu0", 1.0E-8 );
	mu1 = mygetfield( prm, "mu1", 1.0E8 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	sufficientDecreaseCoeff = mygetfield( prm, "sufficientDecreaseCoeff", 0.01 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( isrealscalar(sufficientDecreaseCoeff) );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
		assert( 0.0 < sufficientDecreaseCoeff );
		assert( sufficientDecreaseCoeff < 1.0 );
	end
	
	
	if (1)
	vecG0 = vecG;
	matH0 = matH;
	vecR = [ ...
	  eps050*(xExt-0.5*(max(xVals)+min(xVals)))/((max(xVals)-min(xVals))^2); ...
	  p-2.0 ];
	matR = [ ...
	  eps050/((max(xVals)-min(xVals))^2), 0.0; ...
	  0.0, 1.0 ];
	iterCount = 0;
	while (1)
		iterCount++;
		switch (iterCount)
		case 1
			nu = 0.0;
		case 2
			nu = mu0;
		otherwise
			nu *= muStep;
		end
		if ( nu > mu1 )
			msg_main( verbLev, thisFile, __LINE__, sprintf( "Failed to normalize Hessian (%g).", mu1 ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		%
		clear matG;
		clear matH;
		echo__nu = nu
		vecG = vecG0 + nu*vecR
		matH = matH0 + nu*matR
		if ( rcond(matH)>eps075 )
			break;
		end
	end
	clear vecG0;
	clear matH0;
	clear matR;
	clear vecR;
	end
	
	
	%if ( matH(1,1) <= eps075*sum(sum(abs(matH))) )
	%	msg_warn( verbLev, thisFile, __LINE__, "Hessian (1,1) is near-singular; using a hack." );
	%	matH(1,1) += norm(vecG);
	%end
	%if ( matH(2,2) <= eps075*sum(sum(abs(matH))) )
	%	msg_warn( verbLev, thisFile, __LINE__, "Hessian (2,2) is near-singular; using a hack." );
	%	matH(2,2) += norm(vecG);
	%end
	
	rch = rcond(matH);
	if ( rch <= eps075 )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "rcond(matH) = %g; this case is not supported.", rch ) );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	matD = diag(abs(diag(matH)));
	funch_deltaOmegaModel = @(vecDelta)( vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta );
	%
	%
	%
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
		if ( omega_trial >= omega )
			continue;
		end
		if ( abs(omega_trial - omega) > sufficientDecreaseCoeff*funch_deltaOmegaModel(vecDelta_trial) );
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reduced omega from %g to %g with mu %g.", omega, omega_trial, mu_trial ) );
			xExt_new = xExt_trial;
			p_new = p_trial;
			retCode = RETCODE__SUCCESS;
			return;
		end
	end
end
