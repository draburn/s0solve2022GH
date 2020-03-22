function [ retCode, datOut ] = vizPlane( ...
  rvecS1Vals, vecV1, rvecS2Vals, vecV2, ...
  matC, colmap, ...
  curveDat, ...
  prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "vizPlane";
	startTime = time();
	retCode = RETCODE__NOT_SET;
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	sizeX = size(vecV1,1);
	sizeC = size(colmap,1);
	numS1Vals = size(rvecS1Vals,2);
	numS2Vals = size(rvecS2Vals,2);
	assert( 2 <= sizeX );
	assert( 3 <= numS1Vals );
	assert( 3 <= numS2Vals );
	assert( isrealarray(rvecS1Vals,[1,numS1Vals]) );
	assert( isrealarray(rvecS2Vals,[1,numS2Vals]) );
	assert( isrealarray(vecV1,[sizeX,1]) );
	assert( isrealarray(vecV2,[sizeX,1]) );
	assert( isrealarray(matC,[numS1Vals,numS2Vals]) );
	assert( isrealarray(colmap,[sizeC,3]) );
	assert( 0.0 <= colmap(:,:) );
	assert( 1.0 >= colmap(:,:) );
	%
	numCurves = size(curveDat,2);
	if ( 1 <= numCurves )
		assert( issize(curveDat,[1,numCurves]) );
	end
	%
	minNumCLevs = mygetfield( prm, "minNumCLevs", 4 )
	maxNumCLevs = mygetfield( prm, "maxNumCLevs", min([32,round(sizeC/2)]) )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	[ matS1, matS2 ] = ndgrid( rvecS1Vals, rvecS2Vals );
	%
	% Values of matC should ideally run from 2 to sizeC-1.
	matCMod = round(cap( matC, 1, sizeC ));
	cModMax = max(max(matCMod))
	cModMin = min(min(matCMod))
	numCLevs = median([ minNumCLevs, maxNumCLevs, cModMax-cModMin ])
	size(matS1)
	size(matS2)
	size(matCMod)
	contourf( matS1, matS2, matCMod, numCLevs );
	%image( rvecS2Vals, rvecS1Vals, matCMod );
	%image( matCMod );
	colormap( colmap )
	grid on;
	%
return;
end
