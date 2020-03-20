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
		vecYN = matH \ vecG;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__NEWTON;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYN * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		%
		clear vecYN;
	end
	%
	%
	if (1)
		assert( ~isempty(matV) );
		assert( sizeX == sizeF );
		vecTemp = matV' * vecF0;
		fTemp = vecTemp' * matH * vecTemp;
		if ( fTemp > 0.0 )
			vecYP = vecTemp * ((vecTemp'*vecG)/fTemp);
		else
			vecYP = 0.0*vecTemp;
		end
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__PICARD;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYP * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		%
		clear vecYP;
		clear fTemp;
		clear vecTemp;
	end
	%
	%
	if (1)
		vecTemp = vecG;
		fTemp = vecTemp' * matH * vecTemp;
		assert( 0.0 < fTemp );
		vecYG = vecTemp * ((vecTemp'*vecG)/fTemp);
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__GRADDIR;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYG * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		%
		clear vecYG;
		clear fTemp,
		clear vecTemp;
	end
	%
	%
	if (1)
		muSclA = max(vecDiagH);
		muSclB = min(vecDiagH(vecDiagH>0.0));
		powA = 1.0;
		powB = 1.0;
		muScl = ( (muSclA^powA) * (muSclB^powB) )^(1.0/(powA+powB));
		matL = muScl * matI;
		matA = matH - matL;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__LEVCURVE;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		  nuDummy*( (matL+(nuDummy*matA)) \ vecG )  );
		datOut.curveDat(curveIndex).funchYIsLinear = false;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = false;
		% But, could use eig() to support multiArg.
		%
		clear matA;
		clear matL;
		clear muScl;
	end
	%
	%
return;
end

%!test
%!	test_genCurveFunch;
