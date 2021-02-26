clear;
tic();
msg("newtoptim20210225_calc2",__LINE__,"Performing calculations...");
setprngstates(0);
%
if (1)
	xLo = 0.0;
	xHi = 2.0;
	numXVals = 51;
	numTrials = 1E5;
	%
	a0 =  1.0;
	a1 =  0.3;
	b0 = -1.0;
	b1 =  0.3;
	absC0 = 0.0;
	absC1 = 0.0;
else
	%xLo = 0.0;
	%xHi = 4.0;
	xLo = 1.3;
	xHi = 2.0;
	numXVals = 5;
	%numTrials = 1E4;
	numTrials = 2E7;
	%
	a0 = 3.0;
	a1 = 1.0;
	b0 = -1.0;
	b1 = 1.0;
	absC0 = 1.0;
	absC1 = 1.0;
endif
%
XVALS_DIMENSION = 1;
TRIALS_DIMENSION = 2;
%
aVals = a0 + a1*randn(1,numTrials);
bVals = b0 + b1*randn(1,numTrials);
cVals = abs( absC0 + absC1*randn(1,numTrials) );
matR = randn(numXVals,numTrials);
%
vecX = linspace( xLo, xHi, numXVals )';
%
matF = repmat(aVals,[numXVals,1]) + vecX * bVals + repmat(cVals,[numXVals,1]) .* matR;
% This misses the fact that we do know f exactly at several values of x.
%
toc();
