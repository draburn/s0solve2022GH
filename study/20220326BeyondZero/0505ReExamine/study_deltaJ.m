clear;
setprngstates();
numFigs = 0;
%
sizeX = 10; sizeF = sizeX;
vecF = randn(sizeF,1);
matJ0 = randn(sizeF,sizeX);
matDJ = 0.1*abs(randn(sizeF,sizeX));
% matJ = matJ0 + matDJ.*(randn(sizeF,sizeX));
%
numNuVals = 101;
nuVals = linspace(0.0,1.0,numNuVals);
matH0 = matJ0'*matJ0;
matDH = diag(diag(matDJ'*matDJ));
%%%matDH = matDJ'*matDJ;
vecMG = -( matJ0'*vecF );
vecXVals = zeros(sizeX,numNuVals);
for n=1:numNuVals
	matR = chol( matH0 + nuVals(n)*matDH );
	vecXVals(:,n) = matR \ ( matR'\vecMG );
endfor
%
numTrials = 10000;
matRTrials = randn(sizeX,sizeX,numTrials);
vecFValsTrials = zeros(sizeF,numNuVals,numTrials);
for n=1:numTrials
	vecFValsTrials(:,:,n) = vecF + ( (matJ0+(matRTrials(:,:,n).*matDJ)) * vecXVals ); % Automatic broadcasting over nuVals.
endfor
omegaValsTrials = 0.5*sum(vecFValsTrials.^2,1);
omegaValsAvg = sum(omegaValsTrials,3) / numTrials;
omegaValsSqAvg = sum(omegaValsTrials.^2,3) / numTrials;
omegaValsVar = sqrt( omegaValsSqAvg - (omegaValsAvg.^2) );
%
numFigs++; figure(numFigs);
semilogy( ...
  nuVals, omegaValsAvg - omegaValsVar, '^-', ...
  nuVals, omegaValsAvg, '-', ...
  nuVals, omegaValsAvg + omegaValsVar, 'v-' );
grid on;
xlabel( "nu" );
ylabel( "omega" );
