function [ numFigs, retCode, datOut ] = vizLLMCurves( ...
  funchF, vecX, matV, matW, numPts, vecXSecret, prm=[], datIn=[] );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
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
	assert( isrealarray(vecXSecret,[sizeX,1]) );
	%
	vecF = funchF(vecX);
	assert( isrealarray(vecF,[sizeF,1]) );
	%
	assert( isrealscalar(numPts) );
	assert( 1<=numPts );
	%
	numFigsOffset = mygetfield( prm, "numFigsOffset", 0 );
	%
	mrkList0 = ['+x^v<>sdpho'];
	mrkList = mrkList0(1);
	mszList = [ 5 ];
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO CALCULATIONS.
	%
	[ curveDat, retCode, datOut ] = getLLMCurves( vecF, matV, matW, numPts );
	numCurves = max(size(curveDat));
	colMap = 0.7*jet(numCurves);
	for n=2:numCurves
		mrkList = [ mrkList; mrkList0(mod(n,max(size(mrkList0)))) ];
		mszList = [ mszList; mszList(n-1)+1 ];
	end
	vecXStarLLN = vecX - (matV*((matW'*matW)\(matW'*vecF)));
	for n=1:numCurves
		matY = curveDat(n).matY;
		thisNumPts = size(matY,2);
		matX = repmat(vecX,[1,thisNumPts]) + (matV*matY);
		matResLLM = repmat(vecF,[1,thisNumPts]) + (matW*matY);
		%matRes = funchF( repmat(vecX,[1,thisNumPts]) + (matV*matY) );
		for j=1:thisNumPts
			matRes(:,j) = funchF( vecX + (matV*matY(:,j)) );
		end
		matDistLLMR = repmat(vecXStarLLN,[1,thisNumPts]) - matX;
		matDistSecret = repmat(vecXSecret,[1,thisNumPts]) - matX;
		myDat(n).matX = matX;
		myDat(n).deltaNorm = sqrt(sum(matY.^2,1));
		myDat(n).distLLMR = sqrt(sum(matDistLLMR.^2,1));
		myDat(n).distSecret = sqrt(sum(matDistSecret.^2,1));
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
	strLegend = [ steptype2str(curveDat(1).curveType) ];
	for n=2:numCurves
		strLegend = [ strLegend; steptype2str(curveDat(n).curveType) ];
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAKE PLOTS
	%
	numFigs++; figure(numFigs+numFigsOffset);
	plot( ...
	  myDat(1).deltaNorm, myDat(1).omega, ...
	  'o-', 'color', colMap(1,:), ...
	  'marker', mrkList(1), ...
	  'markerSize', mszList(1) );
	hold on;
	for n=2:numCurves
		plot( ...
		  myDat(n).deltaNorm, myDat(n).omega, ...
		  'o-', 'color', colMap(n,:), ...
		  'marker', mrkList(n), ...
		  'markerSize', mszList(n) );
	end
	hold off;
	xlabel( "Step Size" );
	ylabel( "Omega" );
	title( "Omega vs Step Size" );
	grid on;
	legend( strLegend );
	%
	%
	numFigs++; figure(numFigs+numFigsOffset);
	plot( ...
	  myDat(1).deltaNorm, myDat(1).distSecret, ...
	  'o-', 'color', colMap(1,:), ...
	  'marker', mrkList(1), ...
	  'markerSize', mszList(1) );
	hold on;
	for n=2:numCurves
		plot( ...
		  myDat(n).deltaNorm, myDat(n).distSecret, ...
		  'o-', 'color', colMap(n,:), ...
		  'marker', mrkList(n), ...
		  'markerSize', mszList(n) );
	end
	hold off;
	xlabel( "Step Size" );
	ylabel( "Distance to True Root" );
	title( "Step Size vs Dist TR" );
	grid on;
	legend( strLegend );
	%
	return;
	%
	%
	numFigs++; figure(numFigs+numFigsOffset);
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
	xlabel( "Step Size" );
	ylabel( "OmegaLLM" );
	title( "OmegaLLM vs Step Size" );
	grid on;
	legend( strLegend );
	%
	%
	numFigs++; figure(numFigs+numFigsOffset);
	plot( ...
	  myDat(1).deltaNorm, myDat(1).distLLMR, ...
	  'o-', 'color', colMap(1,:), ...
	  'marker', mrkList(1), ...
	  'markerSize', mszList(1) );
	hold on;
	for n=2:numCurves
		plot( ...
		  myDat(n).deltaNorm, myDat(n).distLLMR, ...
		  'o-', 'color', colMap(n,:), ...
		  'marker', mrkList(n), ...
		  'markerSize', mszList(n) );
	end
	hold off;
	xlabel( "Step Size" );
	ylabel( "Distance to LLM Root" );
	title( "Step Size vs Dist LLMR" );
	grid on;
	legend( strLegend );
	%
	%
	%
	if ( 2==sizeK )
	numXVals = 50;
	numYVals = 51;
	forceSquare = true;
	cx = 1;
	cy = 2;
	%
	funchFLLM = @(x)( repmat(vecF,[1,size(x,2)]) + (matW*(matV'*x)) );
	%
	xMin = min(myDat(1).matX(cx,:));
	xMax = max(myDat(1).matX(cx,:));
	yMin = min(myDat(1).matX(cy,:));
	yMax = max(myDat(1).matX(cy,:));
	for n=2:numCurves
		xMin = min([ xMin, min(myDat(n).matX(cx,:)) ]);
		xMax = max([ xMax, max(myDat(n).matX(cx,:)) ]);
		yMin = min([ yMin, min(myDat(n).matX(cy,:)) ]);
		yMax = max([ yMax, max(myDat(n).matX(cy,:)) ]);
	end
	if (forceSquare)
		xVar = max([ xMax-xMin, yMax-yMin ])/2.0;
		yVar = max([ xMax-xMin, yMax-yMin ])/2.0;
	else
		xVar = (xMax-xMin)/2.0;
		yVar = (yMax-yMin)/2.0;
	end
	xAvg = (xMin+xMax)/2.0;
	yAvg = (yMin+yMax)/2.0;
	xLo = xAvg - (1.5*xVar);
	xHi = xAvg + (1.5*xVar);
	yLo = yAvg - (1.5*yVar);
	yHi = yAvg + (1.5*yVar);
	%
	matH = matW' * matW;
	vecG = matW' * vecF;
	vecXNStep = vecX - (matV*(matH\vecG));
	clear vecG;
	clear math;
	%
	xVals = xLo + ((xHi-xLo)*(0:numXVals-1)/(numXVals-1.0));
	yVals = yLo + ((yHi-yLo)*(0:numYVals-1)/(numYVals-1.0));
	[ xGrid, yGrid ] = meshgrid( xVals, yVals );
	zVals = repmat(vecXNStep,[1,numXVals*numYVals]);
	zVals(cx,:) = reshape(xGrid,1,[]);
	zVals(cy,:) = reshape(yGrid,1,[]);
	%
	for n=1:size(zVals,2)
		fVals(:,n) = funchF( zVals(:,n) );
		fLLMVals(:,n) = funchFLLM( zVals(:,n) );
	end
	omegaLLMVals = 0.5*sum(fLLMVals.^2,1);
	omegaLLMGrid = reshape( omegaLLMVals, numYVals, numXVals );
	%
	numFigs++; figure(numFigs+numFigsOffset);
	contour(xGrid,yGrid,sqrt(omegaLLMGrid),50);
	title( sprintf("OmegaLLM vs %g, %g", cx, cy) );
	if (forceSquare)
		axis square;
	end
	hold on;
	for n=1:numCurves
		plot( ...
		  myDat(n).matX(cx,:), myDat(n).matX(cy,:), ...
		  'o-', 'color', colMap(n,:), ...
		  'marker', mrkList(n), ...
		  'markerSize', mszList(n) );
	end
	hold off;
	grid on;
	%
	%
	omegaVals = 0.5*sum(fVals.^2,1);
	omegaGrid = reshape( omegaVals, numYVals, numXVals );
	%
	numFigs++; figure(numFigs+numFigsOffset);
	contour(xGrid,yGrid,sqrt(omegaGrid),50);
	title( sprintf("Omega vs %g, %g", cx, cy) );
	if (forceSquare)
		axis square;
	end
	hold on;
	for n=1:numCurves
		plot( ...
		  myDat(n).matX(cx,:), myDat(n).matX(cy,:), ...
		  'o-', 'color', colMap(n,:), ...
		  'marker', mrkList(n), ...
		  'markerSize', mszList(n) );
	end
	hold off;
	grid on;
	else
	numFigs++; figure(numFigs+numFigsOffset);
	clf();
	numFigs++; figure(numFigs+numFigsOffset);
	clf();
	end
	%
	%
	%
retCode = RETCODE__SUCCESS;
datOut = [];
return;
end

%!test
%!	test_vizLLMCurves
