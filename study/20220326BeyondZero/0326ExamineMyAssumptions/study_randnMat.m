clear;
tic();
setprngstates();
numFigs = 0;
%
sizeX = 2;
sizeF = 2;
%matJ0 = randn(sizeF,sizeX);
matJ1 = abs(randn(sizeF,sizeX));
matJ0 = zeros(sizeF,sizeX);
%matJ1 = ones(sizeF,sizeX);
%
numTrials = 1000000;
matJTrials = repmat(matJ0,[1,1,numTrials]) + repmat(matJ1,[1,1,numTrials]).*randn(sizeF,sizeX,numTrials);
matJAvg = sum(matJTrials,3)/numTrials
matJTJAvg = zeros(sizeF,sizeX);
for n=1:numTrials
	matJTJAvg += matJTrials(:,:,n)' * matJTrials(:,:,n);
end
matJTJAvg /= numTrials
diag(diag(matJ1'*matJ1))
%
toc();
