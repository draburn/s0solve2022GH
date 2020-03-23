
clear;
commondefs;
thisFile = "test_studyPt_genCurveDat";
tic();
numFigs = 0;
%
sizeK = 100;
sizeX = 100;
sizeF = 100;
%
%seedPrm = demoFunc0101_genSeedPrm("lin-easy");
%seedPrm = demoFunc0101_genSeedPrm("easy");
seedPrm = demoFunc0101_genSeedPrm("moderate");
%seedPrm = demoFunc0101_genSeedPrm("xhard");
%randState = mod(round(time),1E6);
%randState = 0;
%randState = 677832;
%randState = 952523; % Try with 5x100x100xSTEPTYPE__GRADDIR.
%randState = 88; % Try with 100x100x100x"moderate" STEPTYPE__LEVCURVE_SCALED.
%randState = 953150; % Try with 100x100x100x"moderate" STEPTYPE__LEVCURVE_SCALED.
%randState = 955286;
randState = 955618;
echo__randStat = randState
seedPrm.randState = randState;
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%
%
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
%
vecX0 = zeros(seedPrm.sizeX,1);
vecF0 = funchF(vecX0);
matJ = funchJ(vecX0);
matV = eye( sizeX, sizeK );
matW = matJ * matV;
matH = matW' * matW;
vecG = -matW' * vecF0;
vecXSecret = funcPrm.x0;
%
prm = [];
prm.vecXSecret = funcPrm.x0;
[ studyPtDat, retCode, datOut ] = studyPt( ...
  funchF, vecX0, matV, matW, prm );
toc;
