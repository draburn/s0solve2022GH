function [ s, p, bigF0, bigF1, retCode, datOut ] = extFit__findFit( ...
  s0, p0, xVals, fVals, wVals, prm=[] )
	commondefs;
	thisFile = "extFit__findFit";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", false );
	datOut = [];
	%
	numPts = size(xVals,2);
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( isrealarray(xVals,[1,numPts]) );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealarray(wVals,[1,numPts]) );
	end
	%
	%
	%
	% Process parameters...
	sMin = mygetfield( prm, "sMin", [] );
	sMax = mygetfield( prm, "sMax", [] );
	pMin = mygetfield( prm, "pMin", 0.0 );
	pMax = mygetfield( prm, "pMax", [] );
	%
	% Defaults...
	omegaTolCnvg = eps100*0.5*sum(wVals.*(fVals.^2));
	%omegaTolCnvg = eps150*0.5*sum(wVals.*(fVals.^2)); % If you want to be very accurate.
	omegaTolSucc = eps025*0.5*sum(wVals.*(fVals.^2));
	deltaSThresh = eps075*(max(xVals)-min(xVals));
	if ( ~isempty(sMin)  &&  ~isempty(sMax)  )
		deltaSThresh = max([ deltaSThresh, eps050*(sMax-sMin) ]);
	end
	deltaPThresh = eps050;
	if ( ~isempty(pMin)  &&  ~isempty(pMax)  )
		deltaPThresh = max([ deltaPThres, eps050*(pMax-pMin) ]);
	end
	omegaRelThresh = eps050;
	numIterLimit = 100;
	%
	% More params...
	omegaTolCnvg = mygetfield( prm, "omegaTolCnvg", omegaTolCnvg );
	omegaTolSucc = mygetfield( prm, "omegaTolSucc", omegaTolSucc );
	deltaSThresh = mygetfield( prm, "deltaSThresh", deltaSThresh );
	deltaPThresh = mygetfield( prm, "deltaPThresh", deltaPThresh );
	omegaRelThresh = mygetfield( prm, "omegaRelThresh", omegaRelThresh );
	numIterLimit = mygetfield( prm, "numIterLimit", numIterLimit );
	if (doChecks)
		assert( isrealscalar(omegaTolCnvg) );
		assert( isrealscalar(omegaTolSucc) );
		assert( isrealscalar(deltaSThresh) );
		assert( isrealscalar(deltaPThresh) );
		assert( isrealscalar(omegaRelThresh) );
		assert( 0.0 < omegaTolCnvg );
		assert( omegaTolCnvg <= omegaTolSucc );
		assert( 0.0 < deltaSThresh );
		assert( 0.0 < deltaPThresh );
		assert( 0.0 < omegaRelThresh );
		assert( omegaRelThresh < 1.0 );
		assert( isposintscalar(numIterLimit) );
	end
	%
	prm_calcAtPt = mygetfield( prm, "prm_calcAtPt", [] );
	%
	prm_findStep = mygetfield( prm, "prm_findStep", [] );
	prm_findStep.sMin = sMin;
	prm_findStep.sMax = sMax;
	prm_findStep.pMin = pMin;
	prm_findStep.pMax = pMax;
	%
	%
	%
	% DO WORK
	%
	[ rhoVals0, bigF00, bigF10, omega0, retCode ] = extFit__calcAtPt( ...
	  s0, p0, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__calcAtPt() failed for initial guess (%d).", retCode ) );
		return;
	end
	retCode = RETCODE__NOT_SET;
	if (doChecks)
		assert( isrealscalar(omega0) );
		assert( 0.0 <= omega0 );
	end
	%
	%
	%
	s = s0;
	p = p0;
	bigF0 = bigF00;
	bigF1 = bigF10;
	omega = omega0;
	numIter = 0;
	omegaPrev = -1.0;
	while (1)
		msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
		  "   %3d,   %11.3e, %11.3e,   %11.3e", ...
		  numIter, s, p, omega ) );
		% Check success condition first.
		if ( omega <= omegaTolCnvg )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", omega, omegaTolCnvg ) );
			retCode = RETCODE__SUCCESS;
			break;
		end
		%
		% Check imposed stop conditions.
		if ( numIter >= numIterLimit )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached numIterLimit (%d).", numIterLimit ) );
			retCode = RETCODE__IMPOSED_STOP; % Unless changed below.
			break;
		end
		%
		%
		%
		[ sNew, pNew, retCode, datOut_findStep ] = extFit__findStep( ...
		  s, p, xVals, fVals, wVals, prm_findStep );
		if ( RETCODE__SUCCESS != retCode )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "__findStep() was unsuccessful (%d).", retCode ) );
			% Leave retCode as-is, unless changed below.
			break;
		end
		if (doChecks)
			assert( isrealscalar(sNew) );
			assert( isrealscalar(pNew) );
			if (!isempty(sMin))
				assert( sMin <= sNew );
			end
			if (!isempty(sMax))
				assert( sNew <= sMax );
			end
			if (!isempty(pMin))
				assert( pMin <= pNew );
			end
			if (!isempty(pMax))
				assert( pNew <= pMax );
			end
		end
		numIter++;
		[ rhoValsNew, bigF0New, bigF1New, omegaNew, retCode ] = extFit__calcAtPt( ...
		  sNew, pNew, xVals, fVals, wVals, prm_calcAtPt );
		if ( RETCODE__SUCCESS ~= retCode )
			% DRaburn 2021.08.08:
			% This should actually never happen, since __findStep should have already checked this point.
			msg_error( verbLev, thisFile, __LINE__, "__calcAtPt() failed." );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			break;
		end
		retCode = RETCODE__NOT_SET;
		if (doChecks)
			assert( isrealscalar(omegaNew) );
			assert( 0.0 <= omegaNew );
		end
		if ( omegaNew >= omega )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to decrease omega (%g >= %g).", omegaNew, omega ) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN; % Unless changed below.
			break;
		end
		%
		%
		% Accept new values.
		sPrev = s;
		pPrev = p;
		bigF0Prev = bigF0;
		bigF1Prev = bigF1;
		omegaPrev = omega;
		s = sNew;
		p = pNew;
		bigF0 = bigF0New;
		bigF1 = bigF1New;
		omega = omegaNew;
		%
		if ( abs(omegaPrev-omega) <= abs(omegaPrev*omegaRelThresh) )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to decrease omega sufficiently (%g, %g).", ...
			  abs(omegaPrev-omega), abs(omegaPrev*omegaRelThresh) ) );
			retCode = RETCODE__IMPOSED_STOP; % Unless changed below.
			break;
		end
		if (  abs(s-sPrev)<deltaSThresh  &&  abs(p-pPrev)<deltaPThresh  )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Step was too small (%e, %e) vs (%e, %e).", ...
			  s-sPrev, p-pPrev, deltaSThresh, deltaPThresh ) );
			retCode = RETCODE__IMPOSED_STOP; % Unless changed below.
			break;
		end
	end
	msg_main( verbLev, thisFile, __LINE__, sprintf( ...
	  "  %3d,   %11.3e, %11.3e,   %11.3e", numIter, s, p, omega ) );
	if ( omega <= omegaTolSucc )
		retCode = RETCODE__SUCCESS;
	end
return;
end
