clear;
tic();
%setprngstates(14920800);
%setprngstates(70883440);
%setprngstates(95901904);
setprngstates(51041264);
numFigs = 0;
%
sizeX = 2; sizeF = sizeX;
vecF = randn(sizeF,1);
omega = 0.5*sumsq(vecF)
matJ0 = randn(sizeF,sizeX);
%%%matDJ = 0.2*abs(randn(sizeF,sizeX));
matDJ = 10.0*abs(randn(sizeF,sizeX));
% matJ = matJ0 + matDJ.*(randn(sizeF,sizeX));
%
numThetaVals = 101;
%%%thetaVals = 2.0*pi*linspace(0.0,1.0,numThetaVals);
thetaVals = linspace(2.0,4.0,numThetaVals);
vecVVals = [ cos(thetaVals); sin(thetaVals) ];
matVVTVals = zeros(sizeX,sizeX,numThetaVals);
matIMVVTVals = zeros(sizeX,sizeX,numThetaVals);
for m=1:numThetaVals
	matVVTVals(:,:,m) = vecVVals(:,m) * (vecVVals(:,m)');
	matIMVVTVals(:,:,m) = eye(sizeX,sizeX) - matVVTVals(:,:,m);
endfor
vecMF = -vecF;
%
vecX0I = (matJ0'*matJ0)\(matJ0*vecMF);
vecX0P = (matJ0'*matJ0 + diag(diag(matDJ'*matDJ)) )\(matJ0*vecMF);
(vecX0P'*vecX0I)/(norm(vecX0P)*norm(vecX0I))
%
%
numTrials = 100;
matRTrials = randn(sizeX,sizeX,numTrials);
omegaValsTrials = zeros(numThetaVals,numTrials);
vecFValsTrials = zeros(sizeF,numThetaVals,numTrials);
for n=1:numTrials
	matDeltaJ = matRTrials(:,:,n).*matDJ;
	matJ1 = matJ0 + matDeltaJ;
	vecF0ITrials(:,n) = vecF + matJ1*vecX0I;
	vecF0PTrials(:,n) = vecF + matJ1*vecX0P;
	%
	for m=1:numThetaVals;
		matJ = matJ0 + matDeltaJ*matVVTVals(:,:,m);
		vecFValsTrials(:,m,n) = vecF + matJ1*((matJ'*matJ)\(matJ'*vecMF));
		%
		matDJ1 = matDeltaJ*matIMVVTVals(:,:,m);
		vecF2ValsTrials(:,m,n) = vecF + matJ1*((matJ'*matJ+diag(diag(matDJ1'*matDJ1)))\(matJ'*vecMF));
	endfor
endfor
omega0ITrials = 0.5*sum(vecF0ITrials,1);
omega0IAvg = sum(omega0ITrials,2) / numTrials;
omega0ISqAvg = sum(omega0ITrials.^2,2) / numTrials;
omega0IVar = sqrt( omega0ISqAvg - (omega0IAvg.^2) );
%
omega0PTrials = 0.5*sum(vecF0PTrials,1);
omega0PAvg = sum(omega0PTrials,2) / numTrials;
omega0PSqAvg = sum(omega0PTrials.^2,2) / numTrials;
omega0PVar = sqrt( omega0PSqAvg - (omega0PAvg.^2) );
%
%
omegaValsTrials = 0.5*sum(vecFValsTrials.^2,1);
omegaAvgVals = sum(omegaValsTrials,3) / numTrials;
omegaSqAvgVals = sum(omegaValsTrials.^2,3) / numTrials;
omegaVarVals = sqrt( omegaSqAvgVals - (omegaAvgVals.^2) );
%
omega2ValsTrials = 0.5*sum(vecF2ValsTrials.^2,1);
omega2AvgVals = sum(omega2ValsTrials,3) / numTrials;
omega2SqAvgVals = sum(omega2ValsTrials.^2,3) / numTrials;
omega2VarVals = sqrt( omega2SqAvgVals - (omega2AvgVals.^2) );
%
%
%
numFigs++; figure(numFigs);
plot( ...
  thetaVals, omegaAvgVals + omegaVarVals, 'o-', ...
  thetaVals, omegaAvgVals, 'o-', ...
  thetaVals, omega2AvgVals + omega2VarVals, '*-', ...
  thetaVals, omega2AvgVals, '*-' );
grid on;
xlabel( "theta" );
ylabel( "omega" );
%
%
%
numFigs++; figure(numFigs);
plot( ...
  thetaVals, vecVVals'*vecX0I./norm(vecX0I), 'o-', ...
  thetaVals, vecVVals'*vecX0P./norm(vecX0P), 'p-' );
grid on;
xlabel( "theta" );
ylabel( "v' * vecX" );
%
%
%
numFigs++; figure(numFigs);
plot( .....
  thetaVals, (0*thetaVals)+omega0IAvg + omega0IVar, 'o-', ...
  thetaVals, (0*thetaVals)+omega0IAvg, 'o-', ...
  thetaVals, (0*thetaVals)+omega0PAvg + omega0PVar, 'x-', ...
  thetaVals, (0*thetaVals)+omega0PAvg, 'x-' );
grid on;
xlabel( "theta" );
ylabel( "omega" );
%
toc();
