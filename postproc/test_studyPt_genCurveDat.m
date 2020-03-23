clear;
commondefs;
thisFile = "test_studyPt_genCurveDat";
tic();
numFigs = 0;
%
sizeK = 100;
sizeX = 100;
sizeF = 100;
%
seedPrm = demoFunc0101_genSeedPrm("moderate");
%seedPrm = demoFunc0101_genSeedPrm("easy");
%randState = mod(round(time),1E6);
randState = 0;
%randState = 677832;
echo__randStat = randState
seedPrm.randState = randState;
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%
%
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
%
vecX0 = zeros(seedPrm.sizeX,1);
vecF0 = funchF(vecX0);
matJ = funchJ(vecX0);
matV = eye( sizeX, sizeK );
matW = matJ * matV;
matH = matW' * matW;
vecG = -matW' * vecF0;
vecXSecret = funcPrm.x0;
%
prm = [];
[ curveDat, retCode, datOut ] = studyPt_genCurveDat( ...
  funchF, vecX0, matV, matW, matH, vecG, STEPTYPE__LEVCURVE, prm );
indexOfMin = curveDat.indexOfMin;
numNuVals = size(curveDat.rvecNu,2);
toc;
%
rvecDAC = sqrt(sum((curveDat.matY(:,2:end)-curveDat.matY(:,1:end-1)).^2,1));
numFigs++; figure(numFigs);
plot( ...
  (1:numNuVals-1)-0.5, rvecDAC, 'o-', ...
  [0], '.', ...
  (indexOfMin-1)*[1,1], [min(rvecDAC),max(rvecDAC)], 's-', 'linewidth', 3, 'markersize', 10 );
xlabel( "point index" );
ylabel( "distance between points" );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  curveDat.rvecNu, curveDat.rvecDeltaNorm, 'o-', ...
  [0], '.', ...
  curveDat.rvecNu(indexOfMin), curveDat.rvecDeltaNorm(indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
xlabel( "nu" );
ylabel( "deltaNorm" );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  curveDat.rvecDeltaNorm, curveDat.rvecOmegaLin, 'o-', ...
  curveDat.rvecDeltaNorm, curveDat.rvecOmega, 'x-', ...
  curveDat.rvecDeltaNorm(indexOfMin), curveDat.rvecOmega(indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
xlabel( "deltaNorm" );
ylabel( "omega" );
legend( ...
  "lin", ...
  "actual", ...
  "location", "northeast" );
grid on;
