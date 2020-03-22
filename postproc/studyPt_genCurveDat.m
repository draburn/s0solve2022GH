function [ curveDat, retCode, datOut ] = studyPt_genCurveDat( ...
  funchF, vecX0, matV, matH, vecG, stepType, prm=[], datIn=[] )
	%
	warning("Work-in-progress!");
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
	if (VALLEV__HIGH<= valLev)
		assert( matV'*matV, eye(sizeK,sizeK), eps^0.75 );
	end
	%
	matW = mygetfield( prm, "matW", [] );
	if (~isempty(matW))
		assert( isrealarray(matW,[sizeF,sizeK]) );
		if (VALLEV__HIGH<= valLev)
			wSqScale = max(max(matW.^2));
			wScale = sqrt(wSqScale);
			fScale = sqrt(max(vecF0.^2));
			assert( ( (matW'*matW)-matH ) < (eps^0.75)*wSqScale );
			assert( ( (matW'*vecF0) + vecG ) < (eps^0.75)*wScale*fScale );
		end
	end
	%
	numNuValsDesired = mygetfield( prm, "numNuValsDesired", 100 );
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
	vecN = matH \ vecG;
	vecP = matV' * (eye(sizeX,sizeF) * vecF0);
	%
	curveDat.vecDiagH = vecDiagH;
	curveDat.matD = matD;
	curveDat.vecN = vecN;
	curveDat.vecP = vecP;
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
			vecTemp = vecN;
			matS = matI;
		case {STEPTYPE__GRADDIR}
			vecTemp = vecG;
			matS = matI;
		case {STEPTYPE__GRADDIR_SCALED}
			vecTemp = matD \ vecG;
			matS = matD; % Is this right???
		case {STEPTYPE__PICARD}
			vecTemp = vecP;
			matS = matI;
		case {STEPTYPE__PICARD_SCALED}
			vecTemp = matD \ vecP;
			matS = matD; %Is this right???
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
	else
		funchYNormOfNu = @(nuDummy)( sqrt(sum(curveDat.funchYOfNu(nuDummy).^2,1)) );
		curveDat.rvecNuVals = flinspace( 0.0, 1.0, numNuValsDesired, funchYNormOfNu );
	end
	curveDat.numNuVals = size(curveDat.rvecNuVals,2);
	assert( 1 <= curveDat.numNuVals );
	assert(isrealarray(curveDat.rvecNuVals,[1,curveDat.numNuVals]));
	%
	%
return;
end

%!test
%!	test_studyPt_genCurveDat;
