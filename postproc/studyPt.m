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
	parfor n=1:numCurves
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SET OUTPUT.
	%
	studyPtDat.funchF = funchF;
	studyPtDat.vecX0 = vecX0;
	studyPtDat.matV = matV;
	studyPtDat.matW = matW;
	studyPtDat.prm = prm;
	%
	studyPtDat.matH = matH;
	studyPtDat.vecG = vecG;
	%
	studyPtDat.curveDat = curveDat;
	%
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt;
