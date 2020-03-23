function [ studyPtDat, retCode, datOut ] = studyPt( ...
  funchF, vecX0, matV, matW, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "studyPt";
	studyPtDat = [];
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	vecF0 = funchF(vecX0);
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	sizeK = size(matV,2);
	assert( 1 <= sizeK );
	assert( sizeK <= sizeX );
	assert( sizeK <= sizeF );
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]));
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( isrealarray(matW,[sizeF,sizeK]) );
	if (VALLEV__HIGH<= valLev)
		assert( matV'*matV, eye(sizeK,sizeK), eps^0.75 );
	end
	%
	stepTypeList_default = [ ...
	  STEPTYPE__NEWTON, ...
	  STEPTYPE__GRADDIR, ...
	  STEPTYPE__GRADDIR_SCALED, ...
	  STEPTYPE__LEVCURVE, ...
	  STEPTYPE__LEVCURVE_SCALED, ...
	  STEPTYPE__SPECIFIED_VECTOR ];
	%
	genStepFunchPrm = mygetfield( prm, "genStepFunchPrm", [] );
	numNuVals = mygetfield( prm, "numNuVals", 100 );
	stepTypeList = mygetfield( prm, "stepTypeList", stepTypeList_default );
	vecXSecret = mygetfield( prm, "vecXSecret", [] );
	assert( isrealarray(vecXSecret,[sizeX,1]) );
	%
	genStepFunchDatIn = mygetfield( datIn, "genStepFunchDatIn", [] );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	matH = matW' * matW;
	vecG = -matW' * vecF0;
	%
	numCurves = max(size(stepTypeList));
	for n=1:numCurves
		stepType = stepTypeList(n);
		if (STEPTYPE__SPECIFIED_VECTOR==stepType)
			genCurveDat_prm.vecY = matV' * ( vecXSecret - vecX0 );
		else
			genCurveDat_prm.vecY = [];
		end
		curveDat(n) = studyPt_genCurveDat( ...
		  funchF, vecX0, matV, matW, matH, vecG, stepType, genCurveDat_prm );
	end
	%
	%
	numFigs++; figure(numFigs);
	hold off;
	plot( curveDat(1).rvecDeltaNorm, curveDat(1).rvecOmega, "o-" );
	hold on
	for n=2:numCurves
		plot( curveDat(n).rvecDeltaNorm, curveDat(n).rvecOmega, "o-" );
	end
	grid on;
	hold off;
	%
	%
	numFigs++; figure(numFigs);
	hold off;
	plot( 0.0, 0.0, "k+", "linewidth", 3, "markersize", 20 );
	hold on;
	for n=1:numCurves
		stepType = stepTypeList(n);
		m = curveDat(n).indexOfMin;
		plot( ...
		  curveDat(n).matY(1,:), curveDat(n).matY(2,:), "o-", ...
		  curveDat(n).matY(1,end), curveDat(n).matY(2,end), "gx", ...
		  "linewidth", 3, "markersize", 20, ...
		  curveDat(n).matY(1,m), curveDat(n).matY(2,m), "gs", ...
		  "linewidth", 3, "markersize", 20 );
		if (STEPTYPE__SPECIFIED_VECTOR==stepType)
			plot( curveDat(n).matY(1,end), curveDat(n).matY(2,end), ...
			  "rs", "linewidth", 3, "markersize", 20 );
		end
	end
	grid on;
	hold off;
	return;
	
	
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
