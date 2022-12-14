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
	  STEPTYPE__GRADCURVE, ...
	  STEPTYPE__GRADCURVE_SCALED, ...
	  STEPTYPE__SPECIFIED_VECTOR ];
	%
	genStepFunchPrm = mygetfield( prm, "genStepFunchPrm", [] );
	numNuVals = mygetfield( prm, "numNuVals", 100 );
	stepTypeList = mygetfield( prm, "stepTypeList", stepTypeList_default );
	vecXSecret = mygetfield( prm, "vecXSecret", [] );
	if (~isempty(vecXSecret))
		assert( isrealarray(vecXSecret,[sizeX,1]) );
	end
	funchFSupportsMultiArg = mygetfield( prm, "funchFSupportsMultiArg", true );
	%
	genStepFunchDatIn = mygetfield( datIn, "genStepFunchDatIn", [] );
	includeFSOrcale = false;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	matH = matW' * matW;
	vecG = -matW' * vecF0;
	%
	numCurves = max(size(stepTypeList));
	curveColMap = jet(numCurves);
	%
	for n=1:numCurves
		stepType = stepTypeList(n);
		if (STEPTYPE__SPECIFIED_VECTOR==stepType)
			genCurveDat_prm.vecY = matV' * ( vecXSecret - vecX0 );
			curveName = "Oracle";
		else
			genCurveDat_prm.vecY = [];
			curveName = steptype2str(stepType);
		end
		% Do I *really* need "curveDat0" as an intermediate?
		% I *really* can't do curveDat(n).curveName = curveName???
		curveDat0 = studyPt_genCurveDat( ...
		  funchF, vecX0, matV, matW, matH, vecG, stepType, genCurveDat_prm );
		curveDat0.curveName = curveName;
		curveDat0.col = curveColMap(n,:);
		curveDat(n) = curveDat0;
	end
	%
	if (includeFSOrcale)
		% DRaburn 2020.03.25: This is hackish.
		numCurves++;
		genCurveDat_prm.vecY = vecXSecret - vecX0;
		curveName = "FSOracle";
		curveDat0 = studyPt_genCurveDat( funchF, vecX0, ...
		  eye(sizeX,sizeX), zeros(sizeF,sizeX), zeros(sizeX,sizeX), zeros(sizeX,1), ...
		  STEPTYPE__SPECIFIED_VECTOR, genCurveDat_prm );
		curveDat0.curveName = curveName;
		curveDat0.col = [ 0.6, 0.6, 0.6 ];
		curveDat(numCurves) = curveDat0;
	end
	%
	% Generate contour levels and color map;
	% for any given slice, select from among these when plotting.
	haveLevCurve = false;
	for n=1:numCurves
	if (curveDat(n).stepType == STEPTYPE__LEVCURVE )
		haveLevCurve = true;
		omega0 = curveDat(n).rvecOmega(1);
		omegaLo = min(curveDat(n).rvecOmega);
	end
	end
	assert( haveLevCurve );
	assert( omegaLo < omega0 );
	assert( omegaLo >= 0.0 );
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SET OUTPUT.
	%
	studyPtDat.funchFSupportsMultiArg = funchFSupportsMultiArg;
	studyPtDat.funchF = funchF;
	studyPtDat.vecX0 = vecX0;
	studyPtDat.vecF0 = vecF0;
	studyPtDat.matV = matV;
	studyPtDat.matW = matW;
	studyPtDat.prm = prm;
	%
	studyPtDat.matH = matH;
	studyPtDat.vecG = vecG;
	%
	studyPtDat.curveDat = curveDat;
	%
	studyPtDat.omegaLo = omegaLo;
	studyPtDat.omega0 = omega0;
	%
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt;
