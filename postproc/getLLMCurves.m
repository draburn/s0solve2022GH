%

% todo...
%   o Change "delta" to "y"?

function [ curveDat, retCode, datOut ] = getLLMCurves( vecF, matJ, numPts, prm=[], datIn=[] )	
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	getLLMCurves_setCnsts;
	thisFile = "getCurves";
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
	sizeX = size(matJ,2);
	assert( 1 <= sizeX );
	assert( isrealarray(matJ,[sizeF,sizeX]) );
	assert( 2 <= numPts );
	%
	curveTypes_default = [ ...
	  GETCURVES_CURVETYPE__NEWTON, ...
	  GETCURVES_CURVETYPE__PICARD, ...
	  GETCURVES_CURVETYPE__PICARD_SCALED, ...
	  GETCURVES_CURVETYPE__GRADDIR, ...
	  GETCURVES_CURVETYPE__GRADDIR_SCALED, ...
	  GETCURVES_CURVETYPE__LEVCURVE, ...
	  GETCURVES_CURVETYPE__LEVCURVE_SCALED, ...
	  GETCURVES_CURVETYPE__GRADCURVE, ...
	  GETCURVES_CURVETYPE__GRADCURVE_SCALED ];
	curveTypes = mygetfield( prm, "curveTypes", curveTypes_default );
	numCurves = max(size(curveTypes));
	assert( 1 <= numCurves );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	matI = eye(sizeX,sizeX);
	vecG = matJ' * vecF;
	matH = matJ' * matJ;
	matD = diag(diag(matH));
	sVals = (0:numPts-1)/(numPts-1.0);
	%
	for curveIndex = 1 : numCurves
		clear vecTemp;
		clear s0;
		clear funchDelta;
		clear funchF;
		clear muVals;
		clear nuVals;
		switch (curveTypes(curveIndex))
		case {GETCURVES_CURVETYPE__NEWTON}
			vecTemp = -matH\vecG;
			s0 = -(vecTemp'*vecG)/(vecTemp'*matH*vecTemp);
			assert( abs(s0-1.0) < (sizeX^3)*10.0*(eps^0.75) );
			vecTemp *= s0;
			curveDat(curveIndex).matDelta = vecTemp * sVals;
			curveDat(curveIndex).strType = "Newton";
		case {GETCURVES_CURVETYPE__PICARD}
			vecTemp = -vecF;
			s0 = -(vecTemp'*vecG)/(vecTemp'*matH*vecTemp);
			vecTemp *= s0;
			curveDat(curveIndex).matDelta = vecTemp * sVals;
			curveDat(curveIndex).strType = "Picard";
		case {GETCURVES_CURVETYPE__PICARD_SCALED}
			vecTemp = -matD\vecF;
			s0 = -(vecTemp'*vecG)/(vecTemp'*matH*vecTemp);
			vecTemp *= s0;
			curveDat(curveIndex).matDelta = vecTemp * sVals;
			curveDat(curveIndex).strType = "PicardScl";
		case {GETCURVES_CURVETYPE__GRADDIR}
			vecTemp = -vecG;
			s0 = -(vecTemp'*vecG)/(vecTemp'*matH*vecTemp);
			vecTemp *= s0;
			curveDat(curveIndex).matDelta = vecTemp * sVals;
			curveDat(curveIndex).strType = "GradDir";
		case {GETCURVES_CURVETYPE__GRADDIR_SCALED}
			vecTemp = -matD\vecG;
			s0 = -(vecTemp'*vecG)/(vecTemp'*matH*vecTemp);
			vecTemp *= s0;
			curveDat(curveIndex).matDelta = vecTemp * sVals;
			curveDat(curveIndex).strType = "GradDirScl";
		case {GETCURVES_CURVETYPE__LEVCURVE}
			matL = max(diag(matD)) * matI;
			matA = matH-matL;
			funchY = @(nu)( -nu*((matL+(nu*matA))\vecG) );
			funchF = @(nu)( sqrt(sum((funchY(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, numPts, funchF );
			for n=1:max(size(nuVals))
				curveDat(curveIndex).matDelta(:,n) = funchY(nuVals(n));
			end
			curveDat(curveIndex).strType = "Leveneberg";
			clear funchF;
			clear funchY;
			clear nuVals;
			clear matA;
			clear matL;
		otherwise
			error(sprintf( "Value of curveTypes(%g) is invalid (%g).", ...
			  curveIndex, curveTypes(curveIndex) ));
		end
	end
	%
	retCode = RETCODE__SUCCESS;
	datOut = [];
return;
end

%!test
%!	commondefs;
%!	getLLMCurves_setCnsts;
%!	thisFile = "test getLLMCurvs";
%!	vecF = [1;2]
%!	matJ = [1,2;3,5]
%!	numPts = 20;
%!	curveTypes = [ ...
%!	  GETCURVES_CURVETYPE__NEWTON, ...
%!	  GETCURVES_CURVETYPE__PICARD, ...
%!	  GETCURVES_CURVETYPE__PICARD_SCALED, ...
%!	  GETCURVES_CURVETYPE__GRADDIR, ...
%!	  GETCURVES_CURVETYPE__GRADDIR_SCALED ];
%!	prm.curveTypes = curveTypes;
%!	[ curveDat, retCode, datOut ] = getLLMCurves( vecF, matJ, numPts, prm );
