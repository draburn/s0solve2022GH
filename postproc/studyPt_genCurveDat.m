function [ curveDat, retCode, datOut ] = studyPt_genCurveDat( ...
  funchF, vecX0, matV, matW, matH, vecG, stepType, prm=[], datIn=[] )
	%
	% Expect...
	%   matV has orthonormal columns: matV' * matV = eye(sizeK,sizeK);
	%   matW = matJ * matV;
	%   matH = matW' * matW;
	%   vecG = -matW' * vecF0;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "studyPt_genCurveDat";
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
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	assert( isrealscalar(valLev) );
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
		%
		wSqScale = max(max(matW.^2));
		wScale = sqrt(wSqScale);
		fScale = sqrt(max(vecF0.^2));
		assert( ( (matW'*matW)-matH ) < (eps^0.75)*wSqScale );
		assert( ( (matW'*vecF0) + vecG ) < (eps^0.75)*wScale*fScale );
	end
	%
	numNuValsDesired = mygetfield( prm, "numNuValsDesired", 20 );
	%
	curveDat.vecX0 = vecX0;
	curveDat.vecF0 = vecF0;
	curveDat.matV = matV;
	curveDat.matH = matH;
	curveDat.vecG = vecG;
	curveDat.stepType = stepType;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PRE-WORK.
	%
	vecDiagH = diag(matH);
	assert( 0.0 < max(vecDiagH) );
	if ( 0.0 == min(vecDiagH) )
		msg_warn( verbLev, thisFile, __LINE__, "Warning: matH has a diagonal zero." );
	end
	if (VALLEV__HIGH<= valLev)
		rcondH = rcond(matH);
		if ( eps^0.75 > rcondH )
			msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
			  "Warning: rcond(matH) is very small (%g).", rcondH ) );
		end
	end
	matD = diag(vecDiagH);
	matI = eye(sizeK,sizeK);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO MAIN WORK.
	%
	%
	switch( stepType )
	case {STEPTYPE__LEVCURVE}
		muSclA = max(vecDiagH);
		muSclB = min(vecDiagH(vecDiagH>0.0));
		powA = 1.0;
		powB = 1.0;
		muScl = ( (muSclA^powA) * (muSclB^powB) )^(1.0/(powA+powB));
		matL = muScl * matI;
		matA = matH - matL;
		%
		curveDat.funchYOfNu = @(nuDummy)( nuDummy*( (matL+(nuDummy*matA)) \ vecG )  );
		curveDat.funchYIsLinear = false;
		curveDat.funchYSupportsMultiArg = false;
		curveDat.matS = matI;
		%
		%
	case {STEPTYPE__LEVCURVE_SCALED}
		matA = matH - matD;
		%
		curveDat.funchYOfNu = @(nuDummy)( nuDummy*( (matD+(nuDummy*matA)) \ vecG )  );
		curveDat.funchYIsLinear = false;
		curveDat.funchYSupportsMultiArg = false;
		curveDat.matS = matD;
		warning("curveDat.matS may be wrong!");
		%
		%
	case {STEPTYPE__GRADCURVE}
		error( "Not implemented!" );
		%
		%
	case {STEPTYPE__GRADCURVE_SCALED}
		error( "Not implemented!" );
		%
		%
	otherwise
		% Handle all linear cases here...
		switch (stepType)
		case {STEPTYPE__NEWTON}
			vecTemp = matH \ vecG;
			matS = matI;
		case {STEPTYPE__GRADDIR}
			vecTemp = vecG;
			matS = matI;
		case {STEPTYPE__GRADDIR_SCALED}
			vecTemp = matD \ vecG;
			matS = matD; % Is this right???
			warning("curveDat.matS may be wrong, not that it matters.");
		case {STEPTYPE__PICARD}
			vecTemp = matV' * (eye(sizeX,sizeF) * vecF0);
			matS = matI;
		case {STEPTYPE__PICARD_SCALED}
			vecTemp = matD \ (matV' * (eye(sizeX,sizeF) * vecF0));
			matS = matD; %Is this right???
			warning("curveDat.matS may be wrong, not that it matters.");
		otherwise
			error(sprintf( "Invalid value of stepType (%d).", stepType ));
		end
		%
		fTemp = vecTemp' * matH * vecTemp;
		if (STEPTYPE__NEWTON==stepType)
			assert( abs(fTemp-vecTemp'*vecG) < eps^0.75 );
		end
		if ( 0.0 < fTemp )
			vecY = vecTemp * ((vecTemp'*vecG)/fTemp);
		else
			vecY = vecTemp * 0.0;
		end
		%
		curveDat.funchYOfNu = @(nuDummy)( vecY * nuDummy );
		curveDat.funchYIsLinear = true;
		curveDat.funchYSupportsMultiArg = true;
		curveDat.matS = matS; % Doesn't matter because linear.
	end
	%
	%
	if (curveDat.funchYIsLinear)
		curveDat.rvecNuVals = linspace( 0.0, 1.0, numNuValsDesired );
		numNuVals = size(curveDat.rvecNuVals,2);
		if (curveDat.funchYSupportsMultiArg)
			curveDat.matY = curveDat.funchYOfNu(curveDat.rvecNuVals);
		else
			for n=1:numNuVals
				curveDat.matY(:,n) = curveDat.funchYOfNu(curveDat.rvecNuVals(n));
			end
		end
	else
		[ curveDat.rvecNuVals, retCode_dac, datOut_dac ] = daclinspace( ...
		  0.0, 1.0, numNuValsDesired, curveDat.funchYOfNu );
		numNuVals = size(curveDat.rvecNuVals,2);
		curveDat.matY = datOut_dac.matY;
	end
	assert( 1 <= numNuVals );
	assert(isrealarray(curveDat.rvecNuVals,[1,numNuVals]));
	assert(isrealarray(curveDat.matY,[sizeK,numNuVals]));
	%
	curveDat.matDelta = matV * curveDat.matY;
	curveDat.matX = repmat(vecX0,[1,numNuVals]) + curveDat.matDelta;
	assert(isrealarray(curveDat.matX,[sizeX,numNuVals]));
	%
	curveDat.matFLin = repmat(vecF0,[1,numNuVals]) + (matW * curveDat.matY);
	funchFSupportsMultiArg = mygetfield( prm, "funchFSupportsMultiArg", true );
	if (funchFSupportsMultiArg)
		curveDat.matF = funchF(curveDat.matX);
	else
		for n=1:numNuVals
			curveDat.matF(:,n) = funchF(curveDat.matX(:,n));
		end
	end
	assert(isrealarray(curveDat.matF,[sizeF,numNuVals]));
	%
	curveDat.rvecDeltaNorm = sqrt(sum(curveDat.matDelta.^2,1));
	curveDat.rvecOmegaLin = 0.5*sum(curveDat.matFLin.^2,1);
	curveDat.rvecOmega = 0.5*sum(curveDat.matF.^2,1);
	%
	% Find nuOfOmegaMin. ... minscan()?
	% Add this nu to rvecNu?
return;
end

%!test
%!	test_studyPt_genCurveDat;
