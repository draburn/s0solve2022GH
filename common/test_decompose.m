clear;
tic();
sizeX = 10;
vecX0 = zeros(sizeX,1);
vecU1 = randn(sizeX,1);
vecU2 = randn(sizeX,1);
[ vecV1, vecV2, matX, rvecD1, rvecD2, numD1Vals, numD2Vals ] = spanspace( vecX0, vecU1, vecU2 );
vecU3 = randn(sizeX,1);
[ rvecC1, normR1 ] = decompose( vecU1, [vecV1,vecV2] )
[ rvecC2, normR2 ] = decompose( vecU2, [vecV1,vecV2] )
[ rvecC3, normR3 ] = decompose( vecU3, [vecV1,vecV2] )
toc();
