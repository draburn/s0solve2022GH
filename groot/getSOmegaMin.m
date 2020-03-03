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
	assert( isrealscalar(sLo) );
	assert( isrealscalar(sHi) );
	assert( isrealscalar(numPts1) );
	assert( isrealscalar(numPtsX) );
	assert( isrealscalar(numIterLimit) );
	assert( sLo < sHi );
	assert( 4 <= numPts1 );
	assert( 4 <= numPtsX );
	assert( 1 <= numIterLimit );
	%
	funchOmegaOfF_default = @(f)( 0.5*sum(f.^2,1) );
	funchNormOfDelta_default = @(delta)( sqrt(sum(delta.^2,1)) );
	funchOmegaOfF = mygetfield( prm, "funchOmegaOfF", funchOmegaOfF_default );
	funchNormOfDelta = mygetfield( prm, "funchOmegaOfF", funchNormOfDelta_default );
	%
	funchDeltaSupportsMultiArg = mygetfield( prm, "funchDeltaSupportsMultiArg", false );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PREP WORK.
	%
	%
	doMainLoop = true;
	numIter = 0;
	numPtsDesired = numPts1;
	while (doMainLoop)
		sVals = linspace( sLo, sHi, numPtsDesired );
		numPts = size(sVals,2);
		if (funchDeltaSupportsMultiArg)
			matDelta = funchDeltaOfS(sVals);
		else
			clear matDelta;
			for n=1:numPts1
				matDelta(:,n) = funchDeltaOfS(sVals(n));
			end
		end
		assert( isrealarray(matDelta,[sizeX,numPts]) );
		%
		normVals = funchNormOfDelta(matDelta);
		assert( isrealarray(normVals,[1,numPts]) );
		%
		numIter++;
		if (numIter>=numIterLimit)
			doMainLoop = false;
		end
	end
	%
	sOmegaMin = [];
	datOut = [];
return;
end

%!test
%!	test_getSOmegaMin;
