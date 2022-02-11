clear;
setprngstates();
numFigs = 0;
tic();
%
numPts = 10000;
rPts = randn([1,numPts]);
rAvg = sum(rPts)/numPts
rSqAvg = sumsq(rPts)/numPts
rVarSq = rSqAvg - (rAvg^2)
rVar = sqrt(rVarSq)
%
toc();
