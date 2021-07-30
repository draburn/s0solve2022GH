function [ s, p, retCode, datOut ] = extFit__mainLoop( ...
  s0, p0, xVals, fVals, nExactFit, wVals=[], prm=[] )
	commondefs;
	thisFile = "extFit__mainLoop";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	doChecks = mygetfield( prm, "doChecks", true );
	datOut = [];
	%
	numPts = size(xVals,2);
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isposintscalar(nExactFit) );
		assert( 1 <= nExactFit );
		assert( nExactFit <= numPts );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	%
	%
	% Process parameters...
	%
	% findStep params...
	if ( 1 == nExactFit )
		sMin = min(xVals) - 1.0*(max(xVals)-min(xVals));
		sMax = xVals(2)
	elseif ( numPts == nExactFit )
		sMin = xVals(numPts-1);
		sMax = max(xVals) + 1.0*(max(xVals)-min(xVals));
	else
		if (doChecks)
		if ( isempty(mygetfield( prm, "sMin", [] ))  || isempty(mygetfield( prm, "sMax", [] )) )
		if ( 0.0 < ...
		    ( fVals(nExactFit+1) - fVals(nExactFit) ) ...
		  * ( fVals(nExactFit) - fVals(nExactFit-1) )  )
			msg_warn( verbLev, thisFile, __LINE__, "WARNING: nExactFit is not an exterma." );
			msg_warn( verbLev, thisFile, __LINE__, "  But, default sMin/Max are set on this basis." );
		end
		end
		sMin = xVals(nExactFit-1);
		sMax = xVals(nExactFit+1);
	end
	pMin = min([ p0, 1.0 ]);
	pMax = max([ p0, 20.0 ]);
	prm_findStep = mygetfield( prm, "prm_findStep", [] );
	sMin = mygetfield( prm_findStep, "sMin", sMin );
	sMax = mygetfield( prm_findStep, "sMax", sMax );
	pMin = mygetfield( prm_findStep, "pMin", pMin );
	pMax = mygetfield( prm_findStep, "pMax", pMax );
	prm_findStep.sMin = sMin;
	prm_findStep.sMax = sMax;
	prm_findStep.pMin = pMin;
	prm_findStep.pMax = pMax;
	%prm_findStep.useLevMarq = false;
	%
	% mainLoop params...
	omegaTol = eps075*sqrt(sum((fVals.*wVals).^2));
	deltaSThresh = max([ eps050*(sMax-sMin), eps075*(max(xVals)-min(xVals)) ]);
	deltaPThresh = eps025;
	omegaRelThresh = eps025;
	numIterLimit = 100;
	omegaTol = mygetfield( prm, "omegaTol", omegaTol );
	deltaSThresh = mygetfield( prm, "deltaSThresh", deltaSThresh );
	deltaPThresh = mygetfield( prm, "deltaPThresh", deltaPThresh );
	omegaRelThresh = mygetfield( prm, "omegaRelThresh", omegaRelThresh );
	numIterLimit = mygetfield( prm, "numIterLimit", numIterLimit );
	if (doChecks)
		assert( isrealscalar(omegaTol) );
		assert( isrealscalar(deltaSThresh) );
		assert( isrealscalar(deltaPThresh) );
		assert( isrealscalar(omegaRelThresh) );
		assert( 0.0 < omegaTol );
		assert( 0.0 < deltaSThresh );
		assert( 0.0 < deltaPThresh );
		assert( 0.0 < omegaRelThresh );
		assert( omegaRelThresh < 1.0 );
		assert( isposintscalar(numIterLimit) );
	end
	%
	% calcAtPt params...
	prm_calcAtPt = mygetfield( prm, "prm_calcAtPt", [] );
	[ rhoVals0, bigF00, bigF10, omega0 ] = extFit__calcAtPt( ...
	  s0, p0, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
	if (doChecks)
		assert( isrealscalar(omega0) );
		assert( 0.0 <= omega0 );
	end
	%
	%
	%
	s = s0;
	p = p0;
	omega = omega0;
	numIter = 0;
	omegaPrev = -1.0;
	while (1)
		msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
		  "  %3d,   %11.3e, %11.3e,   %11.3e", ...
		  numIter, s, p, omega ) );
		% Check success condition first.
		if ( omega <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", omega, omegaTol ) );
			retCode = RETCODE__SUCCESS;
			return;
		end
		%
		% Check imposed stop conditions.
		if ( numIter >= numIterLimit )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached numIterLimit (%d).", numIterLimit ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		%
		%
		%
		[ sNew, pNew, retCode, datOut_findStep ] = extFit__findStep( ...
		  s, p, xVals, fVals, nExactFit, wVals, prm_findStep );
		if ( RETCODE__SUCCESS != retCode )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "__findStep() was unsuccessful (%d).", retCode ) );
			retCode = RETCODE__BAD_DEPENDENCY;
			return; 
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
		[ rhoValsNew, bigF0New, bigF1New, omegaNew ] = extFit__calcAtPt( ...
		  sNew, pNew, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
		if (doChecks)
			assert( isrealscalar(omegaNew) );
			assert( 0.0 <= omegaNew );
		end
		%
		if ( omegaPrev > 0.0 )
		if ( omega >= omegaPrev )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to decrease omega (%g >= %g).", omega, omegaPrev ) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			return;
		end
		if ( abs(omegaPrev-omega) <= abs(omegaPrev*omegaRelThresh) )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to decrease omega sufficiently (%g, %g).", omega, omegaPrev ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		end
		%
		%
		numIter++;
		sPrev = s;
		pPrev = p;
		omegaPrev = omega;
		s = sNew;
		p = pNew;
		omega = omegaNew;
		if (  abs(s-sPrev)<deltaSThresh  &&  abs(p-pPrev)<deltaPThresh  )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Step was too small (%e, %e) vs (%e, %e).", ...
			  s-sPrev, p-pPrev, deltaSThresh, deltaPThresh ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
	end
return;
end
