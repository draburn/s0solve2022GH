function [ curveDat, retCode, datOut ] = getLocLinModelCurves( ...
  vecF, matV, matW, numPts, prm=[], datIn=[] );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "getLocLinModelCurves";
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
	assert(isrealarray(vecF,[sizeF,1]));
	sizeX = size(matV,1);
	sizeK = size(matV,2);
	assert(isrealarray(matV,[sizeX,sizeK]));
	assert( sum(abs(matV'*matV-eye(sizeK,sizeK))) < sizeK*(eps^0.75) );
	assert(isrealarray(matW,[sizeF,sizeK]));
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	if ( (sizeF ~= sizeX) || (sizeK ~= sizeX) )
		msg_warn( verbLev, thisFile, __LINE__, ...
		  "WARNING: Subspace curves and non-square functions may not be fully supported." );
	end
	%
	matI = eye(sizeK,sizeK);
	vecG = matW' * vecF;
	matH = matW' * matW;
	matD = diag(diag(matH));
	vecF = matD \ vecG;
	vecN = matH \ vecG;
	hScale = max(max(matD));
	sVals = (0:numPts-1)/(numPts-1.0);
	%
	temp0 = vecG' * (matH * vecG);
	assert( temp0 > 0.0 );
	temp1 = (vecG' * vecG) / temp0;
	assert( temp1 > 0.0 );
	vecDeltaG = -temp1 * matV * vecG;
	curveDat.matDeltaG = vecDeltaG * sVals;
	%
	temp0 = vecF' * (matH * vecF);
	assert( temp0 > 0.0 );
	temp1 = (vecF' * vecG) / temp0;
	assert( temp1 > 0.0 );
	vecDeltaF = -temp1 * matV * vecF;
	curveDat.matDeltaF = vecDeltaF * sVals;
	%
	temp0 = vecN' * (matH * vecN);
	assert( temp0 > 0.0 );
	temp1 = (vecN' * vecG) / temp0;
	assert( abs(temp1-1.0) < (sizeK^2) * (eps^0.75) );
	vecDeltaN = -matV * vecN;
	curveDat.matDeltaN = vecDeltaN * sVals;
	%
	%curveDat.matPicard = ...
	%curveDat.matSPicard = ...
	%
	%curveDat.matLevenberg = ...
	matHScaleI = hScale*matI;
	funchYLev = @(mu_dummy)( (matH+(mu_dummy*matHScaleI))\vecG );
	funchFLev = @(mu_dummy)( sqrt(sum((funchYLev(mu_dummy)).^2)) );
	muValsLev = flinspace( 0.0, numPts^2, numPts, funchFLev );
	for n=1:max(size(muValsLev))
		curveDat.matDeltaL(:,n) = -matV * funchYLev(muValsLev(n));
	end
	%
	matDInvSqrt = diag(1./sqrt(diag(matD)));
	matVScl = matV * matDInvSqrt;
	matHScl = matDInvSqrt * matH * matDInvSqrt;
	vecGScl = matDInvSqrt * vecG;
	%
	%curveDat.matLevMarq = ...
	funchYLevScl = @(mu_dummy)( (matHScl+(mu_dummy*matI))\vecGScl );
	funchFLevScl = @(mu_dummy)( sqrt(sum((funchYLevScl(mu_dummy)).^2)) );
	muValsLevScl = flinspace( 0.0, numPts^2, numPts, funchFLevScl );
	for n=1:max(size(muValsLevScl))
		curveDat.matDeltaLM(:,n) = -matVScl * funchYLevScl(muValsLevScl(n));
	end
	%
	%curveDat.matGradCurve = ...
	[ matPsi, matLambda ] = eig( matH );
	assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < (sizeK^3)*(eps^0.75) );
	assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < (sizeK^3)*(eps^0.75) );
	vecPsiTN = matPsi'*vecN;
	lambdaMin = min(diag(matLambda));
	matSigma = matLambda / lambdaMin;
	funchYGradCurve = @(nu_dummy)( ...
	  vecPsiTN - (diag(nu_dummy.^diag(matSigma))*vecPsiTN) );
	funchFGradCurve = @(nu_dummy)( sqrt(sum((funchYGradCurve(nu_dummy)).^2)) );
	nuValsGradCurve = flinspace( 0.0, 1.0, numPts, funchFGradCurve );
	%%%nuValsGradCurve = ((0:numPts-1)/(numPts-1.0)).^2;
	for n=1:max(size(nuValsGradCurve))
		curveDat.matDeltaGC(:,n) = -matV * (matPsi * funchYGradCurve(nuValsGradCurve(n)));
	end
	%
	%curveDat.matGradScurve = ...
	%
retCode = RETCODE__SUCCESS;
datOut = [];
return;
end

%!test
%!	clear;
%!	commondefs;
%!	thisFile = "test getLocLinModelCurves";
%!	tic();
%!	%
%!	%randnSeed = mod(round(1E6*time()),1E6)
%!	randnSeed = 835470
%!	randn("seed",randnSeed);
%!	sizeF = 2;
%!	sizeX = 2;
%!	sizeK = 2;
%!	vecF = randn(sizeF,1);
%!	matV = eye(sizeX,sizeK);
%!	matW = randn(sizeF,sizeK);
%!	numPts = 50;
%!	%
%!	[ curveDat, retCode, datOut ] = getLocLinModelCurves( vecF, matV, matW, numPts );
%!	assert( RETCODE__SUCCESS == retCode );
%!	matDeltaG = curveDat.matDeltaG;
%!	matDeltaF = curveDat.matDeltaF;
%!	matDeltaN = curveDat.matDeltaN;
%!	matDeltaL = curveDat.matDeltaL;
%!	matDeltaLM = curveDat.matDeltaLM;
%!	matDeltaGC = curveDat.matDeltaGC;
%!	%
%!	matF = repmat(vecF,[1,numPts]);
%!	matFPG = matF + (matW * (matV'*matDeltaG));
%!	matFPF = matF + (matW * (matV'*matDeltaF));
%!	matFPN = matF + (matW * (matV'*matDeltaN));
%!	matFPL = matF + (matW * (matV'*matDeltaL));
%!	matFPLM = matF + (matW * (matV'*matDeltaLM));
%!	matFPGC = matF + (matW * (matV'*matDeltaGC));
%!	deltaGNormVals = sqrt(sum(matDeltaG.^2,1));
%!	deltaFNormVals = sqrt(sum(matDeltaF.^2,1));
%!	deltaNNormVals = sqrt(sum(matDeltaN.^2,1));
%!	deltaLNormVals = sqrt(sum(matDeltaL.^2,1));
%!	deltaLMNormVals = sqrt(sum(matDeltaLM.^2,1));
%!	deltaGCNormVals = sqrt(sum(matDeltaGC.^2,1));
%!	omegaGVals = 0.5*(sum(matFPG.^2,1));
%!	omegaFVals = 0.5*(sum(matFPF.^2,1));
%!	omegaNVals = 0.5*(sum(matFPN.^2,1));
%!	omegaLVals = 0.5*(sum(matFPL.^2,1));
%!	omegaLMVals = 0.5*(sum(matFPLM.^2,1));
%!	omegaGCVals = 0.5*(sum(matFPGC.^2,1));
%!	%
%!	numFigs = 0;
%!	xn = 1;
%!	yn = 2;
%!	colorG = [ 1.0, 0.0, 0.0 ];
%!	colorF = [ 0.0, 0.0, 1.0 ];
%!	colorN = [ 0.0, 0.8, 0.0 ];
%!	colorL = [ 0.8, 0.8, 0.0 ];
%!	colorLM = [ 0.8, 0.4, 1.0 ];
%!	colorGC = [ 0.4, 0.6, 0.6 ];
%!	%
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  deltaGNormVals, omegaGVals, 'o-', 'color', colorG, ...
%!	  deltaFNormVals, omegaFVals, 'x-', 'color', colorF, ...
%!	  deltaNNormVals, omegaNVals, 's-', 'color', colorN, ...
%!	  deltaLNormVals, omegaLVals, '^-', 'color', colorL, ...
%!	  deltaLMNormVals, omegaLMVals, 'v-', 'color', colorLM, ...
%!	  deltaGCNormVals, omegaGCVals, '*-', 'color', colorGC, 'markersize', 15  );
%!	grid on;
%!	%
%!	xMin = max(min([ ...
%!	  matDeltaG(xn,:), ...
%!	  matDeltaF(xn,:), ...
%!	  matDeltaN(xn,:) ]));
%!	xMax = max(max([ ...
%!	  matDeltaG(xn,:), ...
%!	  matDeltaF(xn,:), ...
%!	  matDeltaN(xn,:) ]));
%!	yMin = max(min([ ...
%!	  matDeltaG(yn,:), ...
%!	  matDeltaF(yn,:), ...
%!	  matDeltaN(yn,:) ]));
%!	yMax = max(max([ ...
%!	  matDeltaG(yn,:), ...
%!	  matDeltaF(yn,:), ...
%!	  matDeltaN(yn,:) ]));
%!	if (1)
%!	aMin = min([xMin,yMin]);
%!	aMax = max([xMax,yMax]);
%!	xLo = aMin-0.3*abs(aMax-aMin);
%!	xHi = aMax+0.3*abs(aMax-aMin);
%!	yLo = xLo;
%!	yHi = xHi;
%!	else
%!	xLo = xMin-0.3*abs(xMax-xMin);
%!	xHi = xMax+0.3*abs(xMax-xMin);
%!	yLo = yMin-0.3*abs(yMax-yMin);
%!	yHi = yMax+0.3*abs(yMax-yMin);
%!	end
%!	%
%!	numXVals = 50;
%!	numYVals = 51;
%!	xVals = xLo + ((xHi-xLo)*(0:numXVals-1)/(numXVals-1.0));
%!	yVals = yLo + ((yHi-yLo)*(0:numYVals-1)/(numYVals-1.0));
%!	[ xGrid, yGrid ] = meshgrid( xVals, yVals );
%!	%
%!	xv = 1;
%!	yv = 2;
%!	omegaGrid = 0.5*( ...
%!	   (( (vecF(xv) + (matW(xv,xn)*xGrid) + (matW(xv,yn)*yGrid)) ).^2) ...
%!	 + (( (vecF(yv) + (matW(yv,xn)*xGrid) + (matW(yv,yn)*yGrid)) ).^2) );
%!	%
%!	numFigs++; figure(numFigs);
%!	contour( xGrid, yGrid, sqrt(omegaGrid), 50 );
%!	hold on;
%!	plot( ...
%!	  matDeltaG(xn,:), matDeltaG(yn,:), 'o-', 'color', colorG, ...
%!	  matDeltaF(xn,:), matDeltaF(yn,:), 'x-', 'color', colorF, ...
%!	  matDeltaN(xn,:), matDeltaN(yn,:), 's-', 'color', colorN, ...
%!	  matDeltaL(xn,:), matDeltaL(yn,:), '^-', 'color', colorL, ...
%!	  matDeltaLM(xn,:), matDeltaLM(yn,:), '^-', 'color', colorLM, ...
%!	  matDeltaGC(xn,:), matDeltaGC(yn,:), '*-', 'color', colorGC, 'markersize', 15 );
%!	hold off;
%!	axis equal;
%!	grid on;
%!	%
%!	toc();
