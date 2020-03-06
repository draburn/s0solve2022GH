clear;
commondefs;
thisFile = "test_minScan";
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
prm = [];
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
vecF0 = funchF(vecX0);
matJ0 = funchJ(vecX0);
vecDeltaNewton = -matJ0\vecF0;
funchDeltaOfS = @(s)( s*vecDeltaNewton );
[ sOfMin, retCode, datOut ] = minScan( funchF, vecX0, funchDeltaOfS, prm );
%
toc();
