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
	  GETCURVES_CURVETYPE__GRADIENT, ...
	  GETCURVES_CURVETYPE__GRADIENT_SCALED, ...
	  GETCURVES_CURVETYPE__LEVCURVE, ...
	  GETCURVES_CURVETYPE__LEVCURVE_SCALED, ...
	  GETCURVES_CURVETYPE__GRADCURVE, ...
	  GETCURVES_CURVETYPE__GRADCURVE_SCALED ];
	%%%curveTypes = mygetfield( prm, "curveTypes", curveTypes_default );
	curveTypes = mygetfield( prm, "curveTypes", [GETCURVES_CURVETYPE__NEWTON] );
	numCurves = max(size(curveTypes));
	assert( 1 <= numCurves );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	matI = eye(sizeX,sizeX);
	vecG = matJ' * vecF;
	matH = matJ' * matJ;
	sVals = (0:numPts-1)/(numPts-1.0);
	%
	for curveIndex = 1 : numCurves
		switch (curveTypes(curveIndex))
		case {GETCURVES_CURVETYPE__NEWTON}
			vecDeltaN = -matH\vecG
			curveDat(curveIndex).matDelta = vecDeltaN * sVals;
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
%!	vecF = [1;2];
%!	matJ = [1,2;3,5];
%!	numPts = 20;
%!	[ curveDat, retCode, datOut ] = getLLMCurves( vecF, matJ, numPts );
%!	plot( curveDat(1).matDelta(1,:), curveDat(1).matDelta(2,:), 'o-' );
%!	grid on;
