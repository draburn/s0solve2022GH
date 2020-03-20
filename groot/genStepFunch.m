function [ retCode, datOut ] = genStepFunch( ...
  funchF, vecX0, matW, matV=[], prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "genStepFunch";
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
	vecF0 = funchF(vecX0);
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	sizeK = size(matW,2);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecF0,[sizeF,1]));
	assert( isrealarray(matW,[sizeF,sizeK]) );
	%
	if ( ~isempty(matV) )
		assert( isrealarray(matV,[sizeX,sizeK]) );
		matVTV = matV'*matV;
		assert( matVTV, eye(sizeK,sizeK), (eps^0.75)*max(max(abs(matVTV))) );
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	vecG = -matW'*vecF0;
	matH = matW'*matW;
	vecDiagH = diag(matH);
	assert( 0.0 < min(vecDiagH) );
	matD = diag(vecDiagH);
	matI = eye(sizeK,sizeK);
	%
	%
	curveIndex = 0;
	%
	%
	if (1)
		vecDN = matH \ vecG;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__NEWTON;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecDN * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		%
		clear vecDeltaN;
	end
	%
	%
	if (1)
		vecTemp = vecG;
		fTemp = vecTemp' * matH * vecTemp;
		assert( 0.0 < fTemp );
		vecDG = vecTemp * ((vecTemp'*vecG)/fTemp);
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__GRADDIR;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecDG * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		%
		clear vecTemp;
		clear fTemp,
		clear vecDG;
	end
	%
	%
	if (1)
		muScl = min(vecDiagH(vecDiagH>0.0));
		matL = muScl * matI;
		matA = matH + matL;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__LEVCURVE;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		  nuDummy*( (matA-(nuDummy*matL)) \ vecG )  );
		datOut.curveDat(curveIndex).funchYIsLinear = false;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = false;
		% But, could use eig() to support multiArg.
		%
		clear muScl;
		clear matL;
		clear matA;
	end
	%
	%
return;
end

%!test
%!	test_genCurveFunch;
