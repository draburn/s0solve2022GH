clear;
rand( "seed", 0 );
randn( "seed", 0 );
%
sizeX = 7
numPts = 10
%
matX = randn( sizeX, numPts )
%
tic();
matV = genbasis( matX )
vecX0 = matX(:,1);
matV'*(matX-vecX0)
vecX0 + matV*( matV'*(matX-vecX0) )
toc();
