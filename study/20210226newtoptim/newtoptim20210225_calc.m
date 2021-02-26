clear;
tic();
msg("newtoptim20210225_calc2",__LINE__,"Performing calculations...");
setprngstates(0);
%
if (0)
	xLo = 0.0;
	xHi = 2.0;
	numXVals = 21;
	numTrials = 1E6;
	%
	a0 =  1.0;
	a1 =  1.0;
	b0 = -1.0;
	b1 =  1.0;
	absC0 = 0.3;
	absC1 = 0.3;
else
	xLo = 0.0;
	xHi = 2.0;
	numXVals = 21;
	numTrials = 1E6;
	%
	a0 =  1.0;
	a1 =  0.0;
	b0 = -1.0;
	b1 =  1.0;
	absC0 = 0.0;
	absC1 = 0.0;
endif
%
XVALS_DIMENSION = 1;
TRIALS_DIMENSION = 2;
%
aVals = a0 + a1*randn(1,numTrials);
%!
%%%bVals = b0 + b1*randn(1,numTrials);
%bVals = sign(b0)*exp( log(abs(b0)) + log(b1)*randn(1,numTrials) );
bVals = b0 * ( (1.0+b1).^randn(1,numTrials) );
%!
cVals = abs( absC0 + absC1*randn(1,numTrials) );
matR = randn(numXVals,numTrials);
%
vecX = linspace( xLo, xHi, numXVals )';
%
matF = repmat(aVals,[numXVals,1]) + vecX * bVals + repmat(cVals,[numXVals,1]) .* matR;
% This misses the fact that we do know f exactly at several values of x.
%
% This stuff is technically post, but indp post params and time consuming.
matF_abs = abs(matF);
matF_sign = sign(matF);
matF_sort = sort( matF, TRIALS_DIMENSION );
matF_absSort = sort( matF_abs, TRIALS_DIMENSION );
%
%
toc();
