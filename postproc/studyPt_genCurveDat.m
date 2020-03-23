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
	funchFSupportsMultiArg = mygetfield( prm, "funchFSupportsMultiArg", true );
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Set funchY.
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
		funchYOfNu = @(nuDummy)( nuDummy*( (matL+(nuDummy*matA)) \ vecG )  );
		funchYIsLinear = false;
		funchYSupportsMultiArg = false;
		matS = matI;
		%
		%
	case {STEPTYPE__LEVCURVE_SCALED}
		matA = matH - matD;
		%
		funchYOfNu = @(nuDummy)( nuDummy*( (matD+(nuDummy*matA)) \ vecG )  );
		funchYIsLinear = false;
		funchYSupportsMultiArg = false;
		matS = matD;
		warning("matS may be wrong!");
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
			warning("matS may be wrong, not that it matters.");
		case {STEPTYPE__PICARD}
			vecTemp = matV' * (eye(sizeX,sizeF) * vecF0);
			matS = matI;
		case {STEPTYPE__PICARD_SCALED}
			vecTemp = matD \ (matV' * (eye(sizeX,sizeF) * vecF0));
			matS = matD; %Is this right???
			warning("matS may be wrong, not that it matters.");
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
		funchYOfNu = @(nuDummy)( vecY * nuDummy );
		funchYIsLinear = true;
		funchYSupportsMultiArg = true;
		matS = matS; % Doesn't matter because linear.
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get rvecNu and matY.
	%
	if (funchYIsLinear)
		rvecNu = linspace( 0.0, 1.0, numNuValsDesired );
		numNuVals = size(rvecNu,2);
		if (funchYSupportsMultiArg)
			matY = funchYOfNu(rvecNu);
		else
			for n=1:numNuVals
				matY(:,n) = funchYOfNu(rvecNu(n));
			end
		end
	else
		[ rvecNu, retCode_dac, datOut_dac ] = daclinspace( ...
		  0.0, 1.0, numNuValsDesired, funchYOfNu );
		numNuVals = size(rvecNu,2);
		matY = datOut_dac.matY;
	end
	assert( 1 <= numNuVals );
	assert(isrealarray(rvecNu,[1,numNuVals]));
	assert(isrealarray(matY,[sizeK,numNuVals]));
	%
	msg_notify( verbLev, thisFile, __LINE__, ...
	  "It would be nice to have a check that scaled Y is properly monotonic." );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get matDelta, matX, matF, and rvecOmega.
	%
	matDelta = matV * matY;
	matX = repmat(vecX0,[1,numNuVals]) + matDelta;
	if (funchFSupportsMultiArg)
		matF = funchF(matX);
	else
		for n=1:numNuVals
			matF(:,n) = funchF(matX(:,n));
		end
	end
	assert(isrealarray(matF,[sizeF,numNuVals]));
	rvecOmega = 0.5*sum(matF.^2,1);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get nuOfMin.
	%
	funchXOfNu = @(nu)( repmat(vecX0,size(nu)) + (matV*funchYOfNu(nu)) );
	funchOmegaOfNu = @(nu)( 0.5*sum((funchF(funchXOfNu(nu))).^2,1) );
	[ omegaTemp, indexTemp ] = min(rvecOmega);
	if (1==indexTemp)
		nu0 = rvecNu(1);
		nu1 = rvecNu(3);
	elseif (numNuVals==indexTemp)
		nu0 = rvecNu(end-2);
		nu1 = rvecNu(end);
	else
		nu0 = rvecNu(indexTemp-1);
		nu1 = rvecNu(indexTemp+1);
	end
	prm_minscan.funchOmegaSupportsMultiArg = ...
	  ( funchFSupportsMultiArg && funchYSupportsMultiArg );
	nuOfMin = minscan( nu0, nu1, funchOmegaOfNu, prm_minscan );
	clear nu0;
	clear nu1;
	clear omegaTemp;
	clear indexTemp;
	clear funchOmegaOfNu;
	clear funchXOfNu;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get other ofMin quantities and insert.
	%
	vecYOfMin = funchYOfNu(nuOfMin);
	vecDeltaOfMin = matV * vecYOfMin;
	vecXOfMin = vecX0 + vecDeltaOfMin;
	vecFOfMin = funchF(vecXOfMin);
	omegaOfMin = 0.5*sum(vecFOfMin.^2);
	%
	indexOfMin = 0;
	msg_notify( verbLev, thisFile, __LINE__, ...
	  "Add 'ofMin' to data here?" );
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate other quantiies.
	%
	rvecDeltaNorm = sqrt(sum((matDelta).^2,1));
	matFLin = repmat(vecF0,[1,numNuVals]) + (matW * matY);
	rvecOmegaLin = 0.5*sum(matFLin.^2,1);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COPY TO CURVEDAT.
	%
	%
	curveDat.vecX0 = vecX0;
	curveDat.vecF0 = vecF0;
	curveDat.matV = matV;
	curveDat.matH = matH;
	curveDat.vecG = vecG;
	curveDat.stepType = stepType;
	curveDat.prm = prm;
	%
	curveDat.funchYOfNu = funchYOfNu;
	curveDat.funchYIsLinear = funchYIsLinear;
	curveDat.funchYSupportsMultiArg = funchYSupportsMultiArg;
	curveDat.matS = matS;
	%
	curveDat.indexOfMin = indexOfMin;
	%
	curveDat.rvecNu = rvecNu;
	curveDat.matY = matY;
	curveDat.matDelta = matDelta;
	curveDat.matX = matX;
	curveDat.matF = matF;
	curveDat.rvecOmega = rvecOmega;
	%
	curveDat.rvecDeltaNorm = rvecDeltaNorm;
	curveDat.matFLin = matFLin;
	curveDat.rvecOmegaLin = rvecOmegaLin;
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt_genCurveDat;
