clear;
setprngstates();
numFigs = 0;
%
numTrials = 1000000;
rVals = randn(numTrials,1);
rAvg = sum(rVals)/numTrials,
rSqAvg = sumsq(rVals)/numTrials;
rVarSq = rSqAvg - (rAvg^2);
if ( rVarSq > 0.0 )
	rVar = sqrt(rVarSq);
else
	rVar = 0.0;
endif
[ rAvg, rVar ]
