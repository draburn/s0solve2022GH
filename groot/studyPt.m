function [ retCode, datOut ] = studyPt( ...
  funchF, vecX0, matW, matV=[], vecXSecret=[], prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "studyPt";
	startTime = time();
	retCode = RETCODE__NOT_SET;
	datOut = [];
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
	vecF0 = funchF(vecX0);
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	sizeK = size(matW,2);
	assert( 1 <= sizeK );
	assert( sizeK <= sizeX );
	assert( sizeK <= sizeF );
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]));
	assert( isrealarray(matW,[sizeF,sizeK]) );
	%
	if ( ~isempty(matV) )
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( matV'*matV, eye(sizeK,sizeK), eps^0.75 );
	end
	%
	genStepFunchPrm = mygetfield( prm, "genStepFunchPrm", [] );
	genStepFunchDatIn = mygetfield( datIn, "genStepFunchDatIn", [] );
	numNuVals = mygetfield( prm, "numNuVals", 100 );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	[ retCodeTemp, stepFunchDat ] = genStepFunch( ...
	  funchF, vecX0, matW, matV, vecXSecret, genStepFunchPrm, genStepFunchDatIn );
	msg_retcode( verbLev, thisFile, __LINE__, retCodeTemp );
	vecG = stepFunchDat.vecG;
	matH = stepFunchDat.matH;
	%
	numCurves = size(stepFunchDat.curveDat,2);
	assert( issize(stepFunchDat.curveDat,[1,numCurves]) );
	%
	datOut.curveDat = stepFunchDat.curveDat;
	for n=1:numCurves
		%
		funchDNorm = @(nuDummy)(sqrt(sum( ...
		  (stepFunchDat.curveDat(n).matS * ...
		  stepFunchDat.curveDat(n).funchYOfNu(nuDummy)).^2, 1 )));
		[ rvecNuVals, retCodeTemp ] = flinspace( 0.0, 1.0, numNuVals, funchDNorm );
		msg_retcode( verbLev, thisFile, __LINE__, retCodeTemp );
		if ( RETCODE__SUCCESS == retCodeTemp )
			datOut.curveDat(n).rvecNuVals = rvecNuVals;
		else
			datOut.curveDat(n).rvecNuVals = linspace( 0.0, 1.0, numNuVals );
		end
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FINISH.
	%
	datOut.funchF = funchF;
	datOut.vecX0 = vecX0;
	datOut.vecF0 = vecF0;
	datOut.matW = matW;
	datOut.matV = matV;
	datOut.prm = prm;
	%
	retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt;
