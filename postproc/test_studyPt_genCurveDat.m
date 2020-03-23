
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
%seedPrm = demoFunc0101_genSeedPrm("lin-easy");
%seedPrm = demoFunc0101_genSeedPrm("easy");
%seedPrm = demoFunc0101_genSeedPrm("moderate");
seedPrm = demoFunc0101_genSeedPrm("xhard");
%randState = mod(round(time),1E6);
%randState = 0;
%randState = 677832;
%randState = 952523; % Try with 5x100x100xSTEPTYPE__GRADDIR.
%randState = 88; % Try with 100x100x100x"moderate" STEPTYPE__LEVCURVE_SCALED.
%randState = 953150; % Try with 100x100x100x"moderate" STEPTYPE__LEVCURVE_SCALED.
%randState = 955286;
randState = 955618;
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
indexOfMin = curveDat.indexOfMin
minResultIsIncluded = curveDat.minResultIsIncluded
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
%
% More stuff...
%
numFigs++; figure(numFigs);
viz_vecV = myorth(curveDat.matY(:,end));
viz_rvecX1 = viz_vecV'*curveDat.matY;
viz_matTemp = curveDat.matY - (viz_vecV*viz_rvecX1);
viz_rvecX2 = sqrt(sum(viz_matTemp.^2,1));
plot( ...
  viz_rvecX1, viz_rvecX2, 'o-', ...
  [0], [0], 'k+', 'linewidth', 3, 'markersize', 20, ...
  viz_rvecX1(indexOfMin), viz_rvecX2(indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
%
grid on;
%
%
viz_matV = eye(sizeK,sizeK);
viz_matY = viz_matV' * curveDat.matY;

for n=1:floor(sizeK/2)
numFigs++; figure(numFigs);
i1 = (2*n)-1; i2 = 2*n;
plot( ...
  viz_matY(i1,:), viz_matY(i2,:), 'o-', ...
  [0], [0], 'k+', 'linewidth', 3, 'markersize', 20, ...
  viz_matY(i1,indexOfMin), viz_matY(i2,indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
grid on;
end
return


%
numFigs++; figure(numFigs);
i1 = 1; i2 = 2;
plot( ...
  viz_matY(i1,:), viz_matY(i2,:), 'o-', ...
  [0], [0], 'k+', 'linewidth', 3, 'markersize', 20, ...
  viz_matY(i1,indexOfMin), viz_matY(i2,indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
grid on;
%
numFigs++; figure(numFigs);
i1 = 3; i2 = 4;
plot( ...
  viz_matY(i1,:), viz_matY(i2,:), 'o-', ...
  [0], [0], 'k+', 'linewidth', 3, 'markersize', 20, ...
  viz_matY(i1,indexOfMin), viz_matY(i2,indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
grid on;
%
numFigs++; figure(numFigs);
i1 = 1; i2 = 4;
plot( ...
  viz_matY(i1,:), viz_matY(i2,:), 'o-', ...
  [0], [0], 'k+', 'linewidth', 3, 'markersize', 20, ...
  viz_matY(i1,indexOfMin), viz_matY(i2,indexOfMin), 's', 'linewidth', 3, 'markersize', 10 );
grid on;
