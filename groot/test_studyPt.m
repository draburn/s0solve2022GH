clear;
commondefs;
thisFile = "test_studyPt";
tic();
%
%sizeX = 250;
%sizeX = 100;
sizeX = 2;
sizeF = sizeX;
%seedPrm = demoFunc0101_genSeedPrm("lin-easy");
%seedPrm = demoFunc0101_genSeedPrm("easy");
seedPrm = demoFunc0101_genSeedPrm("moderate");
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
[ vecXSuggested, retCode, datOut ] = studyPt( funchF, vecX0, studyPtPrm );
%
if (1)
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
toc();
