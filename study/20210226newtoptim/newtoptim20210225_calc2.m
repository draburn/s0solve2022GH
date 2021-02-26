clear;
setprngstates(0);
tic();
%
if (1)
	xLo = 0.0;
	xHi = 2.0;
	numXVals = 51;
	numTrials = 1E6;
	%
	a0 =  1.0;
	a1 =  0.1;
	b0 = -1.0;
	b1 =  0.1;
	absC0 = 0.1;
	absC1 = 0.1;
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
vecA = a0 + a1*randn(numTrials,1);
vecB = b0 + b1*randn(numTrials,1);
vecC = abs( absC0 + absC1*randn(numTrials,1) );
matR = randn(numTrials,numXVals);
%
xVals = linspace( xLo, xHi, numXVals );
%
matF = repmat(vecA,[1,numXVals]) + vecB * xVals + repmat(vecC,[1,numXVals]).*matR;
% This misses the fact that we do know f exactly at several values of x.
%
toc();
