function [ sOmegaMin, retCode, datOut ] = getSOmegaMin( ...
  funchFOfX, vecX0, funchDeltaOfS, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "getSOmegaMin";
	startTime = time();
	retCode = RETCODE__NOT_SET;
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchFOfX(vecX0);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]));
	%
	sLo = mygetfield( prm, "sLo", 0.0 );
	sHi = mygetfield( prm, "sHi", 1.0 );
	numPts1 = mygetfield( prm, "numPts1", 201 );
	numPtsX = mygetfield( prm, "numPtsX", 21 );
	numIterLimit = mygetfield( prm, "numIterLimit", 10 );
	sTol = mygetfield( prm, "sTol", eps^0.5 );
	omegaRelTol = mygetfield( prm, "omegaRelTol", eps^0.5 );
	omegaTol = mygetfield( prm, "omegaTol", eps^1.5 );
	assert( isrealscalar(sLo) );
	assert( isrealscalar(sHi) );
	assert( isrealscalar(numPts1) );
	assert( isrealscalar(numPtsX) );
	assert( isrealscalar(numIterLimit) );
	assert( isrealscalar(sTol) );
	assert( isrealscalar(omegaRelTol) );
	assert( sLo < sHi );
	assert( 4 <= numPts1 );
	assert( 4 <= numPtsX );
	assert( 1 <= numIterLimit );
	%
	funchOmegaOfF_default = @(f)( 0.5*sum(f.^2,1) );
	funchOmegaOfF = mygetfield( prm, "funchOmegaOfF", funchOmegaOfF_default );
	%funchNormOfDelta_default = @(delta)( sqrt(sum(delta.^2,1)) );
	%funchNormOfDelta = mygetfield( prm, "funchOmegaOfF", funchNormOfDelta_default );
	%
	funchDeltaSupportsMultiArg = mygetfield( prm, "funchDeltaSupportsMultiArg", false );
	funchFSupportsMultiArg = mygetfield( prm, "funchFSupportsMultiArg", true );
	funchOmegaSupportsMultiArg = mygetfield( prm, "funchOmegaSupportsMultiArg", true );
	%
	omega0 = funchOmegaOfF( vecF0 );
	assert( isrealscalar(omega0) );
	msg_progress( verbLev, thisFile, __LINE__, sprintf( "  omega0 = %g", omega0 ) );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PREP WORK.
	%
	%
	numIter = 0;
	while (true)
		if (0==numIter)
			numPtsDesired = numPts1;
		else
			numPtsDesired = numPtsX;
		end
		sVals = linspace( sLo, sHi, numPtsDesired );
		numPts = size(sVals,2);
		%
		if (funchDeltaSupportsMultiArg)
			matDelta = funchDeltaOfS(sVals);
		else
			clear matDelta;
			for n=1:numPts
				matDelta(:,n) = funchDeltaOfS(sVals(n));
			end
		end
		assert( isrealarray(matDelta,[sizeX,numPts]) );
		%
		if (funchFSupportsMultiArg)
			matF = funchFOfX( repmat(vecX0,[1,numPts]) + matDelta );
		else
			clear matF;
			for n=1:numPts
				matF(:,n) = funchFOfX( vecX0 + matDelta(:,n) );
			end
		end
		assert( isrealarray(matF,[sizeF,numPts]) );
		%
		if (funchOmegaSupportsMultiArg)
			omegaVals = funchOmegaOfF( matF );
		else
			clear omegaVals;
			for n=1:numPts
				omegaVals(n) = funchOmegaOfF( matF(:,n) );
			end
		end
		assert( isrealarray(omegaVals,[1,numPts]) );
		%
		[ omegaMin, nOmegaMin ] = min( omegaVals );
		omegaLocalMax = max( omegaVals );
		sOmegaMin = sVals(nOmegaMin);
		%
		numIter++;
		%
		if ( 0.0 <= reportInterval )
		if ( time() > reportTimePrev + reportInterval )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "  %4d  %6.2f  %15.12f  %15.12f  %15.12f  %19.12e", ...
			   numIter, ...
			   time()-startTime, ...
			   sLo, sOmegaMin, sHi, ...
			   omegaMin )  );
			reportTimePrev = time();
		end
		end
		%
		if ( omegaMin <= omegaTol )
			retCode = RETCODE__SUCCESS;
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged to absolute omega tolerance (%g <= %g). %s", ...
			  omegaMin, omegaTol, retcode2str(retCode) ) );
			break;
		elseif ( (omegaLocalMax-omegaMin) <= omegaRelTol*(omega0-omegaMin) )
			retCode = RETCODE__SUCCESS;
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged to relative omega tolerance ( %g <= %g * ( %g - %g ) ). %s", ...
			 omegaLocalMax - omegaMin, ...
			 omegaRelTol, ...
			 omega0, ...
			 omegaMin, ...
			 retcode2str(retCode) ) );
			break;
		elseif (abs(sHi-sLo)<=sTol)
			retCode = RETCODE__IMPOSED_STOP;
			% Arguably __SUCCESS or __ALGORITHM_BREAKDOWN.
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached step size interval tolerance ( abs(%g - %g) <= %g).", ...
			  sHi, sLo, sTol ) );
			break;
		elseif (numIter>=numIterLimit)
			retCode = RETCODE__IMPOSED_STOP;
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached limit on number of iterations (%d >= %d).", ...
			  numIter, numIterLimit ) );
			break;
		else
			if (1==nOmegaMin)
				sHi = sVals(3);
			elseif (numPts==nOmegaMin)
				sLo = sVals(numPts-2);
			else
				sLo = sVals(nOmegaMin-1);
				sHi = sVals(nOmegaMin+1);
			end
		end
	end
	%
	datOut = [];
return;
end

%!test
%!	test_getSOmegaMin;
