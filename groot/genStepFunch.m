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
		assert( matV'*matV, eye(sizeK,sizeK), eps^0.75 );
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
	if (0)
		vecYN = matH \ vecG;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__NEWTON;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYN * nuDummy );
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		%
		clear vecYN;
	end
	%
	%
	if (0)
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
		%
		clear vecYP;
		clear fTemp;
		clear vecTemp;
	end
	%
	%
	if (0)
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
		%
		clear vecYGScl;
		clear fTemp,
		clear vecTemp;
	end
	%
	%
	if (0)
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
	if (0)
		matA = matH - matD;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__LEVCURVE_SCALED;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		  nuDummy*( (matD+(nuDummy*matA)) \ vecG )  );
		datOut.curveDat(curveIndex).funchYIsLinear = false;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = false;
		% But, could use eig() to support multiArg.
		%
		clear matA;
	end
	%
	%
	if (0)
		[ matPsi, matLambda ] = eig( matH );
		assert( matPsi'*matPsi, matI, eps^0.75 );
		assert( matPsi*(matPsi'), matI, eps^0.75 );
		vecPsiTN = matPsi'*(matH\vecG);
		vecLambdaDiag = diag(matLambda);
		lambdaScale = min(vecLambdaDiag(vecLambdaDiag>0.0));
		assert( lambdaScale > 0.0 );
		vecSigma = vecLambdaDiag/lambdaScale;
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__GRADCURVE;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		  matPsi * ( vecPsiTN - (diag((1.0-nuDummy).^vecSigma)*vecPsiTN) )  );
		datOut.curveDat(curveIndex).funchYIsLinear = false;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = false;
		% Could probably modify to support multiArg, but, not worthwhile.
		%
		clear vecSigma;
		clear lambdaScale;
		clear vecLambdaDiag;
		clear vecPsiTN;
		clear matPsi;
		clear matLambda;
	end
	%
	%
	if (1)
		%vecDiagHMod = vecDiagH + ( (eps*max(vecDiagH)) * (vecDiagH==0.0) );
		%matDModInvSqrt = diag(1.0/sqrt(vecDiagHMod));
		matDModInvSqrt = diag(1.0/sqrt(vecDiagH));
		matHScl = matDModInvSqrt * matH * matDModInvSqrt;
		vecGScl = matDModInvSqrt * vecG;
		[ matPsi, matLambda ] = eig( matHScl );
		assert( matPsi'*matPsi, matI, eps^0.75 );
		assert( matPsi*(matPsi'), matI, eps^0.75 );
		vecPsiTN = matPsi'*(matHScl\vecGScl);
		vecLambdaDiag = diag(matLambda);
		lambdaScaleA = min(vecLambdaDiag(vecLambdaDiag>0.0));
		lambdaScaleB = max(vecLambdaDiag);
		powA = 1.0;
		powB = 0.0;
		lambdaScale = ((lambdaScaleA^powA)*(lambdaScaleB^powB))^(1.0/(powA+powB));
		assert( lambdaScale > 0.0 );
		vecSigma = vecLambdaDiag/lambdaScale;
		matDMISPsi = matDModInvSqrt * matPsi;
		%
		vecYN = (matH\vecG);
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__GRADCURVE_SCALED;
		%datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		%  vecYN - (matDMISPsi*( ((1.0-(nuDummy/1E3)).^vecSigma) .* vecPsiTN ))  );
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		  vecYN - (matDMISPsi*( ((1.0-(nuDummy/1E3)).^vecSigma) .* vecPsiTN ))  );
		%datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		%  matDMISPsi * ( vecPsiTN - (diag((1.0-nuDummy).^vecSigma)*vecPsiTN) )  );
		%datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		%  matDMISPsi * ( vecPsiTN * nuDummy )  );
		%datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( ...
		%  matDMISPsi * ( diag(vecSigma)*vecPsiTN * nuDummy )  );
		datOut.curveDat(curveIndex).funchYIsLinear = false;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = false;
		% Could probably modify to support multiArg, but, not worthwhile.
		%
		clear matDMISPsi;
		clear vecSigma;
		clear lambdaScale;
		clear vecLambdaDiag;
		clear vecPsiTN;
		clear matPsi;
		clear matLambda;
		clear vecGScl;
		clear matHScl;
		clear matDModInvSqrt;
		clear vecDiagHMod;
	end
	%
	%
	if (0)
		vecYSSS = matV' * (vecXSecret-vecX0);
		%
		curveIndex++;
		datOut.curveDat(curveIndex).stepType = STEPTYPE__SECRET;
		datOut.curveDat(curveIndex).funchYOfNu = @(nuDummy)( vecYSSS * nuDummy );
		datOut.curveDat(curveIndex).funchYIsLinear = true;
		datOut.curveDat(curveIndex).funchYSupportsMultiArg = true;
		%
		clear vecYSSS;
	end
	%
	numCurves = curveIndex;
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
	%
	datOut.numCurves = numCurves;
	%
	for n=1:numCurves
	switch (datOut.curveDat(n).stepType)
	case {STEPTYPE__NEWTON}
		datOut.curveDat(n).col = [ 0.7, 0.0, 0.0 ];
	case {STEPTYPE__PICARD}
		datOut.curveDat(n).col = [ 1.0, 0.0, 1.0 ];
	case {STEPTYPE__PICARD_SCALED}
		datOut.curveDat(n).col = [ 0.5, 0.0, 0.5 ];
	case {STEPTYPE__GRADDIR}
		datOut.curveDat(n).col = [ 0.0, 0.9, 0.0 ];
	case {STEPTYPE__GRADDIR_SCALED}
		datOut.curveDat(n).col = [ 0.0, 0.5, 0.0 ];
	case {STEPTYPE__LEVCURVE}
		datOut.curveDat(n).col = [ 0.9, 0.9, 0.0 ];
	case {STEPTYPE__LEVCURVE_SCALED}
		datOut.curveDat(n).col = [ 0.5, 0.5, 0.0 ];
	case {STEPTYPE__GRADCURVE}
		datOut.curveDat(n).col = [ 0.0, 0.0, 1.0 ];
	case {STEPTYPE__GRADCURVE_SCALED}
		datOut.curveDat(n).col = [ 0.0, 0.0, 0.5 ];
	case {STEPTYPE__SECRET}
		datOut.curveDat(n).col = [ 0.5, 0.5, 0.5 ];
	otherwise
		datOut.curveDat(n).col = [ 0.0, 0.0, 0.0 ];
	end
	end
	%
return;
end

%!test
%!	test_genCurveFunch;
