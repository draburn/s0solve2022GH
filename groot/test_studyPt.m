clear;
commondefs;
thisFile = "test_studyPt";
tic();
%
%sizeX = 250;
%sizeX = 100;
sizeX = 10;
sizeF = sizeX;
%seedPrm = demoFunc0101_genSeedPrm("lin-easy");
%seedPrm = demoFunc0101_genSeedPrm("easy");
seedPrm = demoFunc0101_genSeedPrm("moderate");
%randState = mod(round(time),1E6)
%randState = 448932 % Newton better then Lev.
randState = 448959 % All bad, but gradient least bad.
%randState = 454940 % A less bad example with sizeX = 100.
seedPrm.randState = randState;
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
vecX0 = zeros(seedPrm.sizeX,1);
%%%vecX0 = funcPrm.x0+1e-7;
vecF0 = demoFunc0101_eval( vecX0, funcPrm );
studyPtPrm = [];
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
studyPtPrm.funchJ = funchJ;
[ vecDelta_suggested, retCode, datOut ] = studyPt( funchF, vecX0, studyPtPrm );
%
numFigs = 0;
%
if (1)
	numFigs++; figure(numFigs);
	hold off;
	plot( [0,1], [0,0], 'ko-', 'linewidth', 7, 'markersize', 20 );
	hold on;
	%
	vecB = funcPrm.x0 - vecX0;
	[ c, normR ] = decompose( datOut.vecDelta_newton, vecB );
	plot( [0,c], [0,normR], 'bx-', 'linewidth', 6, 'markersize', 18 );
	[ c, normR ] = decompose( datOut.vecDelta_newton_omegaMin, vecB );
	plot( [0,c], [0,normR], 'cx-', 'linewidth', 5, 'markersize', 16 );
	[ c, normR ] = decompose( datOut.vecDelta_gradDir, vecB );
	plot( [0,c], [0,normR], 'g+-', 'linewidth', 4, 'markersize', 14 );
	[ c, normR ] = decompose( datOut.vecDelta_gradDir_omegaMin, vecB );
	plot( [0,c], [0,normR], 'ys-', 'linewidth', 3, 'markersize', 12 );
	[ c, normR ] = decompose( datOut.vecDelta_levenberg_omegaMin, vecB );
	plot( [0,c], [0,normR], 'rp-', 'linewidth', 2, 'markersize', 10 );
	%
	axis equal;
	grid on;
	hold off;
	axpand;
end
%
if (1)
	vecU1 = funcPrm.x0 - vecX0;
	%vecU2 = datOut.vecDelta_levenberg_omegaMin;
	vecU2 = datOut.vecDelta_newton;
	[ vecV1, vecV2, matX, rvecD1, rvecD2, numD1Vals, numD2Vals ] ...
	  = spanspace( vecX0, vecU1, vecU2 );
	numPts = size(matX,2);
	assert( numD1Vals*numD2Vals == numPts );
	matF = funchF(matX);
	rvecOmega = 0.5*sum( matF.^2, 1 );
	%
	gridZ = reshape( rvecOmega, numD2Vals, numD1Vals );
	gridX = reshape( rvecD1, numD2Vals, numD1Vals );
	gridY = reshape( rvecD2, numD2Vals, numD1Vals );
	%
	numFigs++; figure(numFigs);
	hold off;
	contour( gridX, gridY, gridZ.^0.2, 31 );
	hold on;
	plot( [ 0.0, vecV1'*vecU1 ], [ 0.0, vecV2'*vecU1 ], 'ko-', 'linewidth', 4 );
	plot( [ 0.0, vecV1'*vecU2 ], [ 0.0, vecV2'*vecU2 ], 'bx-', 'linewidth', 4 );
	[ c, normR ] = decompose( datOut.vecDelta_newton_omegaMin, [vecV1, vecV2] )
	plot( [ 0.0, c(1) ], [ 0.0, c(2) ], 'k-' );
	[ c, normR ] = decompose( datOut.vecDelta_gradDir_omegaMin, [vecV1, vecV2] )
	plot( [ 0.0, c(1) ], [ 0.0, c(2) ], 'k-' );
	grid on;
	axis equal;
	%
	%
	%
	matF0 = funchF( vecX0 );
	matJ0 = funchJ( vecX0 );
	matFLM = repmat(matF0,[1,numPts]) + ( matJ0*matX );
	rvecOLM = 0.5*sum( matFLM.^2, 1 );
	%
	gridZ = reshape( rvecOLM, numD2Vals, numD1Vals );
	gridX = reshape( rvecD1, numD2Vals, numD1Vals );
	gridY = reshape( rvecD2, numD2Vals, numD1Vals );
	%
	numFigs++; figure(numFigs);
	hold off;
	contour( gridX, gridY, gridZ.^0.2, 31 );
	hold on;
	plot( [ 0.0, vecV1'*vecU1 ], [ 0.0, vecV2'*vecU1 ], 'ko-', 'linewidth', 4 );
	plot( [ 0.0, vecV1'*vecU2 ], [ 0.0, vecV2'*vecU2 ], 'bx-', 'linewidth', 4 );
	[ c, normR ] = decompose( datOut.vecDelta_newton_omegaMin, [vecV1, vecV2] )
	plot( [ 0.0, c(1) ], [ 0.0, c(2) ], 'k-' );
	[ c, normR ] = decompose( datOut.vecDelta_gradDir_omegaMin, [vecV1, vecV2] )
	plot( [ 0.0, c(1) ], [ 0.0, c(2) ], 'k-' );
	grid on;
	axis equal;
end

%
toc();
