clear;
numFigs = 0;
setprngstates();
tic();
%
sizeX = 3;
sizeF = sizeX;
cEye = 1.0;
cDelta = 0.1;
%
matR = randn(sizeF,sizeX);
%
numTrials = 1000000;
%
frobVals = zeros(1,numTrials);
condVals = zeros(1,numTrials);
for n=1:numTrials
	%matR = ( matR + cDelta*randn(sizeF,sizeX) ) / sqrt( 1.0 + cDelta^2 );
	matR = ( matR + cDelta*randn(sizeF,sizeX) ) / sqrt( 1.0 + cDelta^2 );
	frobVals(n) = sqrt( sum(sumsq(matR))/(sizeX*sizeF) );
	matJ = cEye*eye(sizeF,sizeX) + matR;
	condVals(n) = cond(matJ);
endfor
%
numFigs++; figure(numFigs);
plot( frobVals, "o-" );
grid on;
%
numFigs++; figure(numFigs);
semilogy( condVals, "o-" );
grid on;
%
numFigs++; fig0 = figure(numFigs);
hist( frobVals(1:ceil(numTrials/3.0)), floor(sqrt(numTrials)/3.0) );
grid on;
ax0 = axis();
%
numFigs++; figure(numFigs);
hist( frobVals(numTrials-floor(numTrials/3.0):end), floor(sqrt(numTrials)/3.0) );
grid on;
ax1 = axis();
axis([ min([ax0(1),ax1(1)]), max([ax0(2),ax1(2)]), min([ax0(3),ax1(3)]), max([ax0(4),ax1(4)]) ]);
figure(fig0);
axis([ min([ax0(1),ax1(1)]), max([ax0(2),ax1(2)]), min([ax0(3),ax1(3)]), max([ax0(4),ax1(4)]) ]);
%
%
%
numFigs++; fig0 = figure(numFigs);
hist( log(condVals(1:ceil(numTrials/3.0))), floor(sqrt(numTrials)/3.0));
grid on;
ax0 = axis();
%
numFigs++; figure(numFigs);
hist( log(condVals(numTrials-floor(numTrials/3.0):end)), floor(sqrt(numTrials)/3.0) );
grid on;
ax1 = axis();
axis([ min([ax0(1),ax1(1)]), max([ax0(2),ax1(2)]), min([ax0(3),ax1(3)]), max([ax0(4),ax1(4)]) ]);
figure(fig0);
axis([ min([ax0(1),ax1(1)]), max([ax0(2),ax1(2)]), min([ax0(3),ax1(3)]), max([ax0(4),ax1(4)]) ]);

toc();
