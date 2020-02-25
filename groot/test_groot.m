clear;
commondefs;
thisFile = "test_groot";
sizeX = 250;
sizeF = sizeX;
seedPrm = demoFunc0101_genSeedPrm("lin-easy");
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
vecX0 = zeros(seedPrm.sizeX,1);
%%%vecX0 = funcPrm.x0+1e-7;
vecF0 = demoFunc0101_eval( vecX0, funcPrm );
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
prm.funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
[ vecXEnd, retCode, datOut ] = groot( funchF, vecX0, prm );