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
	vecG = matW' * vecF;
	matH = matW' * matW;
	matD = diag(diag(matH));
	vecF = matD \ vecG;
	vecN = matH \ vecG;
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
	assert( abs(temp1-1.0) < (sizeK^2) * 10.0 * eps );
	vecDeltaN = -matV * vecN;
	curveDat.matDeltaN = vecDeltaN * sVals;
	%
	%curveDat.matPicard = ...
	%curveDat.matRiker = ...
	%
	%curveDat.matLevenberg = ...
	%curveDat.matLevMarq = ...
	%curveDat.matGradCurve = ...
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
%!	randn("seed",0);
%!	sizeF = 2;
%!	sizeX = 2;
%!	sizeK = 2;
%!	vecF = randn(sizeF,1);
%!	matV = eye(sizeX,sizeK);
%!	matW = randn(sizeF,sizeK);
%!	numPts = 100;
%!	%
%!	[ curveDat, retCode, datOut ] = getLocLinModelCurves( vecF, matV, matW, numPts );
%!	assert( RETCODE__SUCCESS == retCode );
%!	matDeltaG = curveDat.matDeltaG;
%!	matDeltaF = curveDat.matDeltaF;
%!	matDeltaN = curveDat.matDeltaN;
%!	%
%!	matF = repmat(vecF,[1,numPts]);
%!	matFPG = matF + (matW * (matV'*matDeltaG));
%!	matFPF = matF + (matW * (matV'*matDeltaF));
%!	matFPN = matF + (matW * (matV'*matDeltaN));
%!	deltaGNormVals = sqrt(sum(matDeltaG.^2,1));
%!	deltaFNormVals = sqrt(sum(matDeltaF.^2,1));
%!	deltaNNormVals = sqrt(sum(matDeltaN.^2,1));
%!	omegaGVals = 0.5*(sum(matFPG.^2,1));
%!	omegaFVals = 0.5*(sum(matFPF.^2,1));
%!	omegaNVals = 0.5*(sum(matFPN.^2,1));
%!	%
%!	numFigs = 0;
%!	xn = 1;
%!	yn = 2;
%!	%
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  deltaGNormVals, omegaGVals, 'ro-', ...
%!	  deltaFNormVals, omegaFVals, 'gx-', ...
%!	  deltaNNormVals, omegaNVals, 'bs-' );
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
%!	if (0)
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
%!	  matDeltaG(xn,:), matDeltaG(yn,:), 'ro-', ...
%!	  matDeltaF(xn,:), matDeltaF(yn,:), 'gx-', ...
%!	  matDeltaN(xn,:), matDeltaN(yn,:), 'bs-' );
%!	hold off;
%!	%axis equal;
%!	grid on;
%!	%
%!	toc();
