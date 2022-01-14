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
	commoninit;
	thisFile = "studyPt_genCurveDat";
	%
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
		assert( ( (matW'*matW)-matH ) <= (eps^0.75)*wSqScale );
		assert( ( (matW'*vecF0) + vecG ) <= (eps^0.75)*wScale*fScale );
	end
	%
	numNuValsDesired = mygetfield( prm, "numNuValsDesired", 100 );
	assert( isposintscalar(numNuValsDesired-2) );
	funchFSupportsMultiArg = mygetfield( prm, "funchFSupportsMultiArg", true );
	nuDiffSqThresh = mygetfield( prm, "nuDiffSqThresh", 1.0/(10.0*numNuValsDesired) );
	assert( isrealscalar(nuDiffSqThresh) );
	assert( nuDiffSqThresh > 0.0 );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PRE-WORK.
	%
	vecDiagH = diag(matH);
	assert( 0.0 <= min(vecDiagH) );
	if ( 0.0 == max(vecDiagH) )
		msg_warn( verbLev, thisFile, __LINE__, "Warning: matH has all zero diagonals." );
	elseif ( 0.0 == min(vecDiagH) )
		msg_warn( verbLev, thisFile, __LINE__, "Warning: matH has a zero diagonal." );
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
		%
		%
	case {STEPTYPE__LEVCURVE_SCALED}
		matA = matH - matD;
		%
		funchYOfNu = @(nuDummy)( nuDummy*( (matD+(nuDummy*matA)) \ vecG )  );
		funchYIsLinear = false;
		funchYSupportsMultiArg = false;
		%
		%
	case {STEPTYPE__GRADCURVE}
		[ matPsi, matLambda ] = eig( matH );
		if (VALLEV__HIGH<= valLev)
			assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
			assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
		end
		vecPsiTN = matPsi'*(matH\vecG);
		lambdaMin = min(diag(matLambda));
		matSigma = matLambda / lambdaMin;
		%
		funchYOfNu = @(nu)( matPsi*(vecPsiTN - (diag((1.0-nu).^diag(matSigma))*vecPsiTN)) );
		funchYIsLinear = false;
		funchYSupportsMultiArg = false;
		%
		clear matSigma;
		clear lambdaMin;
		clear vecPsiTN;
		clear matPsi;
		clear matLambda;
		%
		%
	case {STEPTYPE__GRADCURVE_SCALED}
		matDInvSqrt = diag(1./sqrt(diag(matD)));
		matHScl = matDInvSqrt * matH * matDInvSqrt;
		vecGScl = matDInvSqrt * vecG;
		[ matPsi, matLambda ] = eig( matHScl );
		if (VALLEV__HIGH<= valLev)
			assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
			assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
		end
		vecPsiTN = matPsi'*(matHScl\vecGScl);
		lambdaMin = min(diag(matLambda));
		matSigma = matLambda / lambdaMin;
		%
		funchYOfNu = @(nu)( matDInvSqrt*matPsi*(vecPsiTN - (diag((1.0-nu).^diag(matSigma))*vecPsiTN)) );
		funchYIsLinear = false;
		funchYSupportsMultiArg = false;
		%
		clear matSigma;
		clear lambdaMin;
		clear vecPsiTN;
		clear matLambda;
		clear vecGScl;
		clear matHScl;
		clear matDInvSqrt;
		%
		%
	otherwise
		% Handle all linear cases here...
		switch (stepType)
		case {STEPTYPE__NEWTON}
			vecTemp = matH \ vecG;
		case {STEPTYPE__GRADDIR}
			vecTemp = vecG;
		case {STEPTYPE__GRADDIR_SCALED}
			vecTemp = matD \ vecG;
		case {STEPTYPE__PICARD}
			vecTemp = matV' * (eye(sizeX,sizeF) * vecF0);
		case {STEPTYPE__PICARD_SCALED}
			vecTemp = matD \ (matV' * (eye(sizeX,sizeF) * vecF0));
		case {STEPTYPE__SPECIFIED_VECTOR}
			% This is a linear case, but, don't follow normal limit....
			vecTemp = getfield( prm, "vecY" );
			assert( isrealarray(vecTemp,[sizeK,1]) );
		otherwise
			error(sprintf( "Invalid value of stepType (%d).", stepType ));
		end
		%
		fTemp = vecTemp' * matH * vecTemp;
		if (STEPTYPE__NEWTON==stepType)
			%echo__fTemp = fTemp
			%echo__vecTempPTVecG = vecTemp'*vecG
			assert( abs(fTemp-vecTemp'*vecG) < sqrt(eps)*(abs(fTemp)+abs(vecTemp'*vecG)) );
			vecY = vecTemp;
		elseif (STEPTYPE__SPECIFIED_VECTOR == stepType)
			vecY = vecTemp;
		elseif ( 0.0 < fTemp )
			vecY = vecTemp * ((vecTemp'*vecG)/fTemp);
		else
			vecY = vecTemp * 0.0;
		end
		%
		funchYOfNu = @(nuDummy)( vecY * nuDummy );
		funchYIsLinear = true;
		funchYSupportsMultiArg = true;
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
	% Get other ofMin quantities and overwrite.
	%
	vecYOfMin = funchYOfNu(nuOfMin);
	vecDeltaOfMin = matV * vecYOfMin;
	vecXOfMin = vecX0 + vecDeltaOfMin;
	vecFOfMin = funchF(vecXOfMin);
	omegaOfMin = 0.5*sum(vecFOfMin.^2);
	%
	[ nuDiffSqTemp, indexOfMin ] = min( (rvecNu-nuOfMin).^2 );
	if ( nuDiffSqTemp < nuDiffSqThresh )
		if ( 1 == indexOfMin )
			% Do nothing more.
			minResultIsIncluded = false;
		elseif ( numNuVals == indexOfMin )
			% Do nothing more.
			minResultIsIncluded = false;
		else
			% Overwrite the point.
			rvecNu(indexOfMin) = nuOfMin;
			matY(:,indexOfMin) = vecYOfMin(:);
			matDelta(:,indexOfMin) = vecDeltaOfMin(:);
			matX(:,indexOfMin) = vecXOfMin(:);
			matF(:,indexOfMin) = vecFOfMin(:);
			rvecOmega(:,indexOfMin) = omegaOfMin(:);
			minResultIsIncluded = true;
		end
	else
		if ( nuOfMin > rvecNu(indexOfMin) )
			indexOfMin++;
		end
		rvecNu(indexOfMin:numNuVals+1) = [ nuOfMin, rvecNu(indexOfMin:numNuVals) ];
		matY(:,indexOfMin:numNuVals+1) = [ vecYOfMin, matY(:,indexOfMin:numNuVals) ];
		matDelta(:,indexOfMin:numNuVals+1) = [ vecDeltaOfMin, matDelta(:,indexOfMin:numNuVals) ];
		matX(:,indexOfMin:numNuVals+1) = [ vecXOfMin, matX(:,indexOfMin:numNuVals) ];
		matF(:,indexOfMin:numNuVals+1) = [ vecFOfMin, matF(:,indexOfMin:numNuVals) ];
		rvecOmega(indexOfMin:numNuVals+1) = [ omegaOfMin, rvecOmega(indexOfMin:numNuVals) ];
		numNuVals++;
		minResultIsIncluded = true;
	end
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate other quantiies.
	%
	rvecDeltaNorm = sqrt(sum((matDelta).^2,1));
	matFLin = repmat(vecF0,[1,numNuVals]) + (matW * matY);
	rvecOmegaLin = 0.5*sum(matFLin.^2,1);
	rvecDAC = [ 0.0, cumsum(sqrt(sum((matDelta(:,2:end)-matDelta(:,1:end-1)).^2,1))) ];
	% DRaburn 2020.03.24: This measure of distance-along-curve is crude.
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
	%
	curveDat.indexOfMin = indexOfMin;
	curveDat.minResultIsIncluded = minResultIsIncluded;
	curveDat.nuOfMin = nuOfMin;
	curveDat.vecYOfMin = vecYOfMin;
	curveDat.vecDeltaOfMin = vecDeltaOfMin;
	curveDat.vecXOfMin = vecXOfMin;
	curveDat.vecFOfMin = vecFOfMin;
	curveDat.omegaOfMin = omegaOfMin;
	%
	curveDat.rvecIndex = (1:numNuVals);
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
	curveDat.rvecDAC = rvecDAC;
	%
	% Distance between points, distance to nearest neighbor.
	rvecDBP = sqrt(sum((matDelta(:,2:end)-matDelta(:,1:end-1)).^2,1));
	curveDat.rvecDTNN(1) = rvecDBP(1);
	curveDat.rvecDTNN(2:numNuVals-1) = min([ rvecDBP(1:end-1); rvecDBP(2:end) ]);
	curveDat.rvecDTNN(numNuVals) = rvecDBP(end);
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt_genCurveDat;
