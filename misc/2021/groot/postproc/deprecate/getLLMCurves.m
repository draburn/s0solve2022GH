function [ curveDat, retCode, datOut ] = getLLMCurves( vecF, matV, matW, numPts, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "getLLMCurves";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% VALIDATE INPUT
	%
	%
	sizeF = size(vecF,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF,[sizeF,1]) );
	sizeX = size(matV,1);
	sizeK = size(matV,2);
	assert( 1 <= sizeX );
	assert( 1 <= sizeK );
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( sum(abs(matV'*matV-eye(sizeK,sizeK))) < 10.0*(sizeK^3)*(eps^0.75) );
	assert( isrealarray(matW,[sizeF,sizeK]) );
	assert( 2 <= numPts );
	%
	curveTypes_default = [ ...
	  STEPTYPE__NEWTON, ...
	  STEPTYPE__GRADDIR, ...
	  STEPTYPE__LEVCURVE, ...
	  STEPTYPE__GRADCURVE ];
	curveTypes = mygetfield( prm, "curveTypes", curveTypes_default );
	numCurves = max(size(curveTypes));
	assert( 1 <= numCurves );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	matI = eye(sizeK,sizeK);
	vecG = matW' * vecF;
	matH = matW' * matW;
	vecD = diag(matH);
	assert( min(vecD) > 0.0 );
	matD = diag(diag(matH));
	sVals = (0:numPts-1)/(numPts-1.0);
	%
	for curveIndex = 1 : numCurves
		switch (curveTypes(curveIndex))
		case {STEPTYPE__NEWTON}
			vecTemp = -matH\vecG;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			assert( abs(s0-1.0) < (sizeK^3)*10.0*(eps^0.75) );
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			clear denom;
			clear s0;
			clear vecTemp;
		case {STEPTYPE__PICARD}
			vecTemp = -matV' * vecF;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			clear denom;
			clear s0;
			clear vecTemp;
		case {STEPTYPE__PICARD_SCALED}
			vecTemp = -matD\(matV'*vecF);
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			clear denom;
			clear s0;
			clear vecTemp;
		case {STEPTYPE__GRADDIR}
			vecTemp = -vecG;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			clear denom;
			clear s0;
			clear vecTemp;
		case {STEPTYPE__GRADDIR_SCALED}
			vecTemp = -matD\vecG;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			clear denom;
			clear s0;
			clear vecTemp;
		case {STEPTYPE__LEVCURVE}
			matL = max(vecD) * matI;
			matA = matH-matL;
			funchY = @(nu)( -nu*((matL+(nu*matA))\vecG) );
			funchF = @(nu)( sqrt(sum((funchY(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, numPts, funchF );
			for n=1:max(size(nuVals))
				curveDat(curveIndex).matY(:,n) = funchY(nuVals(n));
			end
			clear nuVals;
			clear funchF;
			clear funchY;
			clear matA;
			clear matL;
		case {STEPTYPE__LEVCURVE_SCALED}
			matDInvSqrt = diag(1./sqrt(diag(matD)));
			matHScl = matDInvSqrt * matH * matDInvSqrt;
			vecGScl = matDInvSqrt * vecG;
			matAScl = matHScl-matI;
			funchY = @(nu)( -nu*((matI+(nu*matAScl))\vecGScl) );
			funchF = @(nu)( sqrt(sum((funchY(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, numPts, funchF );
			for n=1:max(size(nuVals))
				curveDat(curveIndex).matY(:,n) = matDInvSqrt * funchY(nuVals(n));
			end
			clear nuVals;
			clear funchF;
			clear funchY;
			clear matAScl;
			clear vecGScl;
			clear matHScl;
			clear matDInvSqrt;
		case {STEPTYPE__GRADCURVE}
			[ matPsi, matLambda ] = eig( matH );
			assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
			assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
			vecPsiTN = matPsi'*(matH\vecG);
			lambdaMin = min(diag(matLambda));
			matSigma = matLambda / lambdaMin;
			funchY = @(nu)( ...
			  vecPsiTN - (diag(nu.^diag(matSigma))*vecPsiTN) );
			funchF = @(nu)( sqrt(sum((funchY(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, numPts, funchF );
			for n=1:max(size(nuVals))
				curveDat(curveIndex).matY(:,n) = -matPsi * funchY(nuVals(n));
			end
			clear nuVals;
			clear funcF;
			clear funcY;
			clear matSigma;
			clear lambdaMin;
			clear vecPsiTN;
			clear matPsi;
			clear matLambda;
		case {STEPTYPE__GRADCURVE_SCALED}
			matDInvSqrt = diag(1./sqrt(diag(matD)));
			matHScl = matDInvSqrt * matH * matDInvSqrt;
			vecGScl = matDInvSqrt * vecG;
			[ matPsi, matLambda ] = eig( matHScl );
			assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
			assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < 10.0*(sizeK^3)*(eps^0.75) );
			vecPsiTN = matPsi'*(matHScl\vecGScl);
			lambdaMin = min(diag(matLambda));
			matSigma = matLambda / lambdaMin;
			funchY = @(nu)( ...
			  vecPsiTN - (diag(nu.^diag(matSigma))*vecPsiTN) );
			funchF = @(nu)( sqrt(sum((funchY(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, numPts, funchF );
			for n=1:max(size(nuVals))
				curveDat(curveIndex).matY(:,n) = -matDInvSqrt * matPsi * funchY(nuVals(n));
			end
			clear nuVals;
			clear funcF;
			clear funcY;
			clear matSigma;
			clear lambdaMin;
			clear vecPsiTN;
			clear matPsi;
			clear matLambda;
			clear vecGScl;
			clear matHScl;
			clear matDInvSqrt;
		otherwise
			error(sprintf( "Value of curveTypes(%g) is invalid (%g).", ...
			  curveIndex, curveTypes(curveIndex) ));
		end
		curveDat(curveIndex).curveType = curveTypes(curveIndex);
	end
	%
	retCode = RETCODE__SUCCESS;
	datOut = [];
return;
end

%!test
%!	test_getLLMCurves
