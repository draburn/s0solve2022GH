function [ xExt, fExt, retCode, datOut ] = findExt_basic( xVals=[1,2,3], fVals=[1,4,9], ...
  xExt_initial=[], pL_initial=2.0, pR_initial=2.0, prm=[], datIn=[] )
	%
	findExt_basic__init;
	thisFile = "findExt_basicAsym";
	%
	return;
	
	
	%
	% DO WORK
	%
	[ rhoVals0, bigF00, bigF10, omega0, retCode ] = extFit__calcAtPt( ...
	  s0, p0, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__calcAtPt() failed for initial guess (%d).", retCode ) );
		s = s0;
		p = p0;
		bigF0 = 0.0;
		bigF1 = 0.0;
		% Return retCode as-is.
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
