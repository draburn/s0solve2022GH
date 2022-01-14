function [ retCode, datOut ] = studyPt( ...
  funchF, vecX0, matW, matV=[], vecXSecret=[], prm=[], datIn=[] )
	%
	error("Under renovation!");
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
	stepTypeList_default = [ ...
	  STEPTYPE__NEWTON, ...
	  STEPTYPE__GRADDIR, ...
	  STEPTYPE__GRADDIR_SCALED, ...
	  STEPTYPE__LEVCURVE, ...
	  STEPTYPE__GRADCURVE ];
	% When sizeK < sizeX, full-space delta is not a "step",
	% considering "steps" to be restricted to matV.
	% So, handle it separately.
	%
	genStepFunchPrm = mygetfield( prm, "genStepFunchPrm", [] );
	numNuVals = mygetfield( prm, "numNuVals", 100 );
	stepTypeList = mygetfield( prm, "stepTypeList", stepTypeList_default );
	%
	genStepFunchDatIn = mygetfield( datIn, "genStepFunchDatIn", [] );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	[ retCodeTemp, stepFunchDat ] = genStepFunch( ...
	  funchF, vecX0, matW, matV, vecXSecret, genStepFunchPrm, genStepFunchDatIn );
	msg_retcode( verbLev, thisFile, __LINE__, retCodeTemp );
	datOut.vecG = stepFunchDat.vecG;
	datOut.matH = stepFunchDat.matH;
	datOut.stepFunchDat_raw = stepFunchDat;
	numCurves_raw = size(stepFunchDat.curveDat,2);
	assert( issize(stepFunchDat.curveDat,[1,numCurves_raw]) );
	%
	numCurves = 0;
	for n=1:numCurves_raw
	if (ismember(stepFunchDat.curveDat(n).stepType,stepTypeList))
		numCurves++;
		datOut.curveDat(numCurves) = stepFunchDat.curveDat(n);
	end
	end
	%
	parfor n=1:numCurves
		datOut.curveDat(n).funchDNorm = @(nuDummy)(sqrt(sum( ...
		  (stepFunchDat.curveDat(n).matS * ...
		  stepFunchDat.curveDat(n).funchYOfNu(nuDummy)).^2, 1 )));
		if (stepFunchDat.curveDat(n).funchYIsLinear)
			datOut.curveDat(n).rvecNuVals = linspace( 0.0, 1.0, numNuVals );
		else
			datOut.curveDat(n).rvecNuVals = flinspace( ...
			  0.0, 1.0, numNuVals, datOut.curveDat(n).funchDNorm );
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
