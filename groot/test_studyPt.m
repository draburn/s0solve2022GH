clear;
commondefs;
thisFile = "test_studyPt";
tic();
%
sizeX = 250;
%sizeX = 100;
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
[ vecXSuggested, retCode, studyPtDatOut ] = studyPt( funchF, vecX0, studyPtPrm );
%
toc();
