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
	assert( 1 <= sizeK );
	assert( sizeK <= sizeX );
	assert( sizeK <= sizeF );
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
	vecXSecret = mygetfield( prm, "vecXSecret", [] );
	if ( ~isempty(vecXSecret) )
		assert( isrealarray(vecXSecret,[sizeX,1]) );
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	vecG = -matW'*vecF0;
	matH = matW'*matW;
	vecDiagH = diag(matH);
	assert( 0.0 < max(vecDiagH) );
	if ( 0.0 == min(vecDiagH) )
		msg_warn( verbLev, thisFile, __LINE__, "Warning: matW has a zero vector." );
	end
	rcondH = rcond(matH);
	if ( eps^0.75 > rcondH )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "Warning: rcond(matH) is very small (%g).", rcondH ) );
	end
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
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).col = [ 0.7, 0.0, 0.0 ];
		%
		clear vecYN;
	end
	%
	%
	if (1)
		if (isempty(matV))
			vecTemp = eye(sizeK,sizeF)*vecF0;
		else
			vecTemp = matV' * (eye(sizeX,sizeF)*vecF0);
		end
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
		datOut.curveDat(curveIndex).col = [ 1.0, 0.0, 1.0 ];
		%
		clear vecYP;
		clear fTemp;
		clear vecTemp;
	end
	%
	%
	if (1)
		if (isempty(matV))
			vecTemp = eye(sizeK,sizeF)*vecF0;
		else
			vecTemp = matV' * (eye(sizeX,sizeF)*vecF0);
		end
		vecTemp = matD\vecTemp;
		fTemp = vecTemp' * matH * vecTemp;
		if ( fTemp > 0.0 )
			vecYP = vecTemp * ((vecTemp'*vecG)/fTemp);
		else
			vecYP = 0.0*vecTemp;
		end
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__PICARD_SCALED;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYP * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		datOut.curveDat(curveIndex).col = [ 0.5, 0.0, 0.5 ];
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
		datOut.curveDat(curveIndex).col = [ 0.0, 0.9, 0.0 ];
		%
		clear vecYG;
		clear fTemp,
		clear vecTemp;
	end
	%
	%
	if (1)
		vecTemp = matD \ vecG;
		fTemp = vecTemp' * matH * vecTemp;
		assert( 0.0 < fTemp );
		vecYGScl = vecTemp * ((vecTemp'*vecG)/fTemp);
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__GRADDIR_SCALED;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYGScl * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		datOut.curveDat(curveIndex).col = [ 0.0, 0.5, 0.0 ];
		%
		clear vecYGScl;
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
		datOut.curveDat(curveIndex).col = [ 0.9, 0.9, 0.0 ];
		% But, could use eig() to support multiArg.
		%
		clear matA;
		clear matL;
		clear muScl;
	end
	%
	%
	if (1)
		matA = matH - matD;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__LEVCURVE_SCALED;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		  nuDummy*( (matD+(nuDummy*matA)) \ vecG )  );
		datOut.curveDat(curveIndex).funchYIsLinear = false;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = false;
		datOut.curveDat(curveIndex).col = [ 0.5, 0.5, 0.0 ];
		% But, could use eig() to support multiArg.
		%
		clear matA;
	end
	%
	%
	if (1)
		vecYSSS = matV' * (vecXSecret-vecX0);
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__SECRET;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYSSS * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		datOut.curveDat(curveIndex).col = [ 0.5, 0.5, 0.5 ];
		%
		clear vecYSSS;
	end
	%
	%
	datOut.funchF = funchF;
	datOut.vecX0 = vecX0;
	datOut.vecF0 = vecF0;
	datOut.matW = matW;
	datOut.matV = matV;
	%
	datOut.vecG = vecG;
	datOut.matH = matH;
	datOut.vecDiagH = vecDiagH;
	datOut.matD = matD;
	datOut.matI = matI;
return;
end

%!test
%!	test_genCurveFunch;
