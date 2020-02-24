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
	assert( abs(temp1-1.0) < eps^0.75 );
	vecDeltaN = -matV * vecN;
	curveDat.matDeltaN = vecDeltaN * sVals;
	%
retCode = RETCODE__SUCCESS;
datOut = [];
return;
end

%!test
%!	commondefs;
%!	thisFile = "test getLocLinModelCurves";
%!	randn("seed",0);
%!	sizeF = 5;
%!	sizeX = 5;
%!	sizeK = 3;
%!	vecF = randn(sizeF,1);
%!	matV = eye(sizeX,sizeK);
%!	matW = randn(sizeF,sizeK);
%!	numPts = 100;
%!	[ curveDat, retCode, datOut ] = getLocLinModelCurves( vecF, matV, matW, numPts );
%!	assert( RETCODE__SUCCESS == retCode );
%!	matF = repmat(vecF,[1,numPts]);
%!	matFPG = matF + (matW * (matV'*curveDat.matDeltaG));
%!	matFPF = matF + (matW * (matV'*curveDat.matDeltaF));
%!	matFPN = matF + (matW * (matV'*curveDat.matDeltaN));
%!	deltaGNormVals = sqrt(sum(curveDat.matDeltaG.^2,1));
%!	deltaFNormVals = sqrt(sum(curveDat.matDeltaF.^2,1));
%!	deltaNNormVals = sqrt(sum(curveDat.matDeltaN.^2,1));
%!	omegaGVals = 0.5*(sum(matFPG.^2,1));
%!	omegaFVals = 0.5*(sum(matFPF.^2,1));
%!	omegaNVals = 0.5*(sum(matFPN.^2,1));
%!	xn = 1;
%!	yn = 2;
%!	numFigs = 0;
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  curveDat.matDeltaG(xn,:), curveDat.matDeltaG(yn,:), 'o-', ...
%!	  curveDat.matDeltaF(xn,:), curveDat.matDeltaF(yn,:), 'x-', ...
%!	  curveDat.matDeltaN(xn,:), curveDat.matDeltaN(yn,:), 's-' );
%!	grid on;
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  deltaGNormVals, omegaGVals, 'o-', ...
%!	  deltaFNormVals, omegaFVals, 'x-', ...
%!	  deltaNNormVals, omegaNVals, 's-' );
%!	grid on;
