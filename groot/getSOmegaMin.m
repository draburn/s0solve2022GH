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
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
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
	sTol = mygetfield( prm, "sTol", 1.0E-8 );
	omegaVarTol = mygetfield( prm, "omegaVarTol", 1.0E-4 );
	omegaTol = mygetfield( prm, "omegaTol", eps^1.5 );
	assert( isrealscalar(sLo) );
	assert( isrealscalar(sHi) );
	assert( isrealscalar(numPts1) );
	assert( isrealscalar(numPtsX) );
	assert( isrealscalar(numIterLimit) );
	assert( isrealscalar(sTol) );
	assert( isrealscalar(omegaVarTol) );
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
		%echo__sRange = [ sLo, sHi ]
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
		if ( omegaMin <= omegaTol )
			retCode = RETCODE__SUCCESS;
			break;
		elseif ( (omegaLocalMax-omegaMin) <= omegaVarTol )
			retCode = RETCODE__SUCCESS;
			break;
		elseif (abs(sHi-sLo)<=sTol)
			retCode = RETCODE__IMPOSED_STOP;
			% Arguably __SUCCESS or __ALGORITHM_BREAKDOWN.
			break;
		elseif (numIter>=numIterLimit)
			retCode = RETCODE__IMPOSED_STOP;
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
