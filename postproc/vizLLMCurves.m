function [ numFigs, retCode, datOut ] = vizLLMCurves( ...
  funchF, vecX, matV, matW, numPts, prm=[], datIn=[] );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	getLLMCurves_setCnsts;
	thisFile = "vizLLMCurve";
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
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	assert( 1 <= sizeK );
	assert( sizeK <= sizeX );
	assert( sizeK <= sizeF );
	assert( isrealarray(vecX,[sizeX,1]) );
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( isrealarray(matW,[sizeF,sizeK]) );
	%
	vecF = funchF(vecX);
	assert( isrealarray(vecF,[sizeF,1]) );
	%
	assert( isrealscalar(numPts) );
	assert( 1<=numPts );
	%
	figsIndex0 = 0;
	%
	curveTypes = [ ...
	  GETCURVES_CURVETYPE__NEWTON, ...
	  GETCURVES_CURVETYPE__LEVCURVE, ...
	  GETCURVES_CURVETYPE__GRADCURVE, ...
	  GETCURVES_CURVETYPE__GRADDIR, ...
	  GETCURVES_CURVETYPE__LEVCURVE_SCALED, ...
	  GETCURVES_CURVETYPE__GRADCURVE_SCALED, ...
	  GETCURVES_CURVETYPE__GRADDIR_SCALED ];
	%
	numCurves = max(size(curveTypes));
	colMap = 0.8*jet(numCurves);
	mrkList0 = ['+x^v<>sdpho'];
	mrkList = mrkList0(1);
	mszList = [ 5 ];
	for n=2:numCurves
		mrkList = [ mrkList; mrkList0(mod(n,max(size(mrkList0)))) ];
		mszList = [ mszList; mszList(n-1)+1 ];
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO "WORK".
	%
	%
	prm.curveTypes = curveTypes;
	[ curveDat, retCode, datOut ] = getLLMCurves( vecF, matV, matW, numPts, prm );
	vecXStarLLN = vecX - (matV*(matW'*matW)\(matW*vecF));
	for n=1:numCurves
		matY = curveDat(n).matY;
		thisNumPts = size(matY,2);
		matX = repmat(vecX,[1,thisNumPts]) + (matV*matY);
		matResLLM = repmat(vecF,[1,thisNumPts]) + (matW*matY);
		matRes = funchF( repmat(vecX,[1,thisNumPts]) + (matV*matY) );
		matDTRLLN = repmat(vecXStarLLN,[1,thisNumPts]) - matX;
		myDat(n).deltaNorm = sqrt(sum(matY.^2,1));
		myDat(n).dtrLLN = sqrt(sum(matDTRLLN.^2,1));
		myDat(n).matResLLM = matResLLM;
		myDat(n).omegaLLM = 0.5*sum(matResLLM.^2,1);
		myDat(n).matRes = matRes;
		myDat(n).omega = 0.5*sum(matRes.^2,1);
		clear matDTRLLN;
		clear matRes;
		clear matResLLM;
		clear thisNumPts;
		clear matX;
		clear matY;
	end
	%
	numFigs = 0;
	strLegend = [ curveDat(1).strType ];
	for n=2:numCurves
		strLegend = [ strLegend; curveDat(n).strType ];
	end
	%
	numFigs++; figure(numFigs+figsIndex0);
	plot( ...
	  myDat(1).deltaNorm, myDat(1).omegaLLM, ...
	  'o-', 'color', colMap(1,:), ...
	  'marker', mrkList(1), ...
	  'markerSize', mszList(1) );
	hold on;
	for n=2:numCurves
		plot( ...
		  myDat(n).deltaNorm, myDat(n).omegaLLM, ...
		  'o-', 'color', colMap(n,:), ...
		  'marker', mrkList(n), ...
		  'markerSize', mszList(n) );
	end
	hold off;
	grid on;
	legend( strLegend );
	%
retCode = RETCODE__SUCCESS;
datOut = [];
return;
end

%!test
%!	clear;
%!	commondefs;
%!	getLLMCurves_setCnsts;
%!	thisFile = "test vizLLMCurvs";
%!	%
%!	vecF = [-2;-3];
%!	matJ = [-1,1;1,5];
%!	funchF = @(v)( repmat(vecF,[1,size(v,2)]) + (matJ*v) );
%!	matV = eye(2,2);
%!	vecX = [0;0];
%!	numPts = 20;
%!	%
%!	matW = matJ*matV;
%!	vizLLMCurves( funchF, vecX, matV, matW, numPts );
