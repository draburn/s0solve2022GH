numFigs = 0;
%
vecXE = (1:sizeX)';
assert( norm(funchFJ(vecXE)) < eps );
%
vecXVals = datOut_lesquj_scan.vecXVals;
vecFVals = datOut_lesquj_scan.vecFVals;
matJVals = datOut_lesquj_scan.matJVals;
numVals = size(matJVals,3);
%
foo1Vals = zeros(1,numVals);
foo2Vals = zeros(1,numVals);
matJVals_secret = zeros(sizeF,sizeX,numVals);
foo1Vals_secret = zeros(1,numVals);
foo2Vals_secret = zeros(2,numVals);
for n=1:numVals
	vecFHat = vecFVals(:,n)/norm(vecFVals(:,n));
	vecJDX = matJVals(:,:,n)*( vecXVals(:,n) - vecXE );
	foo1Vals(n) = vecFHat'*vecJDX;
	foo2Vals(n) = norm( vecJDX - vecFHat*foo1Vals(n) );
	%
	[ blarg, matJVals_secret(:,:,n) ] = funchFJ( vecXVals(:,n) );
	vecJDX_secret = matJVals_secret(:,:,n)*( vecXVals(:,n) - vecXE );
	foo1Vals_secret(n) = vecFHat'*vecJDX_secret;
	foo2Vals_secret(n) = norm( vecJDX_secret - vecFHat*foo1Vals_secret(n) );
endfor
%
numFigs++; figure(numFigs);
plot( ...
  datOut_lesquj_scan.fevalCountVals, foo1Vals, '^-', 'linewidth', 2, 'markersize', 15, ...
  datOut_lesquj_scan.fevalCountVals, foo2Vals, 'v-', 'linewidth', 2, 'markersize', 15, ...
  datOut_lesquj_scan.fevalCountVals, foo1Vals_secret, 's-', 'linewidth', 2, 'markersize', 15, ...
  datOut_lesquj_scan.fevalCountVals, foo2Vals_secret, 'p-', 'linewidth', 2, 'markersize', 15, 'color', [0.3,0.7,0.0] );
grid on;
legend( ...
  "foo1", ...
  "foo2", ...
  "foo1 secret", ...
  "foo2 secret", ...
  "location", "northeast" );
%
%
numFigs++; figure(numFigs);
loglog( datOut_lesquj_scan.sVals, 'o-', 'linewidth', 2, 'markersize', 15 );
grid on;
%
myIter = 150;
%myIter = 60;
vecX = datOut_lesquj_scan.vecXVals(:,myIter);
vecF = datOut_lesquj_scan.vecFVals(:,myIter);
matJ = datOut_lesquj_scan.matJVals(:,:,myIter);
%
collected_vecXVals = datOut_lesquj_scan.iterDat(myIter).collected_vecXVals;
collected_vecFVals = datOut_lesquj_scan.iterDat(myIter).collected_vecFVals;
yVals = (vecXE-vecX)'*(collected_vecXVals-vecX)/sumsq(vecXE-vecX);
omegaVals = sumsq(collected_vecFVals,1)/2.0;
%
vecXPts = datOut_lesquj_scan.iterDat(myIter).lesquj_datOut.vecXPts;
vecFPts = datOut_lesquj_scan.iterDat(myIter).lesquj_datOut.vecFPts;
wPts = datOut_lesquj_scan.iterDat(myIter).lesquj_datOut.wPts;
yPts = (vecXE-vecX)'*(vecXPts-vecX)/sumsq(vecXE-vecX);
omegaPts = sumsq(vecFPts,1)/2.0;
%
vVals = linspace( min(yVals), max(yVals), 101 );
vecXModelVals = vecX + (vecXE-vecX)*vVals;
vecFModelVals = vecF + matJ*(vecXModelVals-vecX);
omegaModelVals = sumsq(vecFModelVals,1)/2.0;
yModelVals = (vecXE-vecX)'*(vecXModelVals-vecX)/sumsq(vecXE-vecX);
assert( reldiff(vVals,yModelVals) < sqrt(eps) );
%
numFigs++; figure(numFigs);
plot( ...
  yVals, omegaVals, 'o', 'linewidth', 2, 'markersize', 15, ...
  yPts, omegaPts, '*', 'linewidth', 2', 'markersize', 10, ...
  yPts(myIter+1:end), omegaPts(myIter+1:end), 'p', 'linewidth', 2', 'markersize', 20, ...
  vVals, omegaModelVals, '-', 'linewidth', 2, ...
  0.0, sumsq(vecF)/2.0, 'x', 'linewidth', 3, 'markersize', 15 );
grid on;
numFigs++; figure(numFigs);
semilogy( ...
  yPts, wPts, 'o', 'linewidth', 2, 'markersize', 15 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  omegaVals, 'o', 'linewidth', 2, 'markersize', 15, ...
  omegaPts, '*', 'linewidth', 2', 'markersize', 10 );
grid on;
numFigs++; figure(numFigs);
plot( ...
  wPts, 'o', 'linewidth', 2, 'markersize', 15 );
grid on;

error( "TODO: Call lesquj from here with generated data,and tweak lesquj to make results good." );
