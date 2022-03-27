clear;
%setprngstates(44953376); % 44953376 is an unusual case.
%setprngstates(67018368); % 67018368 is another unusual case.
setprngstates();
numFigs = 0;
%
%sizeX = 2; sizeF = 2;
%sizeX = 20; sizeF = 20;
sizeX = 100; sizeF = 100;
vecF = randn(sizeF,1);
matJ0 = randn(sizeF,sizeX);
matJ1 = 0.01*abs(randn(sizeF,sizeX));
% matJ = matJ0 + matJ1.*(randn(sizeF,sizeX));
% This "should" have the property that < matJ' * matJ > = matJ0' * matJ0 + diag(daig( matJ1' * matJ1 )), I think?
%
vecDelta1 = -( matJ0'*matJ0 ) \ ( matJ0'*vecF );
vecDelta2 = -( matJ0'*matJ0 + diag(diag(matJ1'*matJ1)) ) \ ( matJ0'*vecF );
%
numTrials = 10000;
matJTrials = repmat(matJ0,[1,1,numTrials]) + repmat(matJ1,[1,1,numTrials]).*randn(sizeF,sizeX,numTrials);
%
vecF1Trials = zeros(sizeF,numTrials);
vecF2Trials = zeros(sizeF,numTrials);
for n=1:numTrials
	vecF1Trials(:,n) = matJTrials(:,:,n)*vecDelta1;
	vecF2Trials(:,n) = matJTrials(:,:,n)*vecDelta2;
end
vecF1Trials += repmat(vecF,[1,numTrials]);
vecF2Trials += repmat(vecF,[1,numTrials]);
%
vecF1Avg = sum(vecF1Trials,2)/numTrials;
vecF2Avg = sum(vecF2Trials,2)/numTrials;
vecF1SqAvg = sum(vecF1Trials.^2,2)/numTrials;
vecF2SqAvg = sum(vecF2Trials.^2,2)/numTrials;
%
vecF1AvgNormSq = sumsq(vecF1Avg);
vecF2AvgNormSq = sumsq(vecF2Avg);
vecF1NormSqAvg = sum(vecF1SqAvg);
vecF2NormSqAvg = sum(vecF2SqAvg);
%
msg( __FILE__, __LINE__, "Results..." );
msg( __FILE__, __LINE__, sprintf( " Option 1:  %10.3e, %10.3e.", sqrt(vecF1AvgNormSq), sqrt(vecF1NormSqAvg) ) );
msg( __FILE__, __LINE__, sprintf( " Option 2:  %10.3e, %10.3e.", sqrt(vecF2AvgNormSq), sqrt(vecF2NormSqAvg) ) );
%
f1NormTrials = sqrt(sum(vecF1Trials.^2,1));
f2NormTrials = sqrt(sum(vecF2Trials.^2,1));
%
numFigs++; figure(numFigs);
hist( f1NormTrials, ceil(sqrt(numTrials)) );
grid on;
title( "f1 norm histogram" );
ax1 = axis();
%
numFigs++; figure(numFigs);
hist( f2NormTrials, ceil(sqrt(numTrials)) );
grid on;
title( "f2 norm histogram" );
ax2 = axis();
%
%figure(1); axis([ 0.0, min([ ax1(2) ax2(2) ]), ax1(3), ax1(4) ]);
%figure(2); axis([ 0.0, min([ ax1(2) ax2(2) ]), ax1(3), ax1(4) ]);
