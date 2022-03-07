numFigs = 0;
vecXE = (1:sizeX)';
vecFE = funchFJ(vecXE);
assert( norm(vecFE) < eps );
%
myIter = 60;
msg( __FILE__, __LINE__, sprintf( "Analyzing iter %d.", myIter ) );
datOut = datOut_lesquj_scan; haveLesquj = true;
%datOut = datOut_broyd_scan; haveLesquj = false;
% Handle step/jupdate dependent stuff here.
if (haveLesquj)
	vecXPts = datOut.iterDat(myIter).lesquj_datOut.vecXPts;
	vecFPts = datOut.iterDat(myIter).lesquj_datOut.vecFPts;
	wPts = datOut.iterDat(myIter).lesquj_datOut.wPts;
else
	vecXPts = datOut.vecXVals(:,1);
	vecFPts = datOut.vecFVals(:,1);
	wPts = 1.0;
endif
%
%
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut.sVals, 'o-', 'linewidth', 2, 'markersize', 10, ...
  sqrt(sumsq(datOut.vecXVals(:,2:end)-datOut.vecXVals(:,1:end-1),1)), 'x-', 'linewidth', 2, 'markersize', 10, ...
  reshape(sqrt(sum(sumsq(datOut.matJVals(:,:,2:end)-datOut.matJVals(:,:,1:end-1),1),2)),1,[]), '^-', 'linewidth', 2, 'markersize', 10 );
grid on;
%
%
%
if (0)
nStride = 2;
numFigs++; figure(numFigs);
semilogy( reshape(sqrt(sum(sumsq( ...
  datOut.matJVals(:,:,1+nStride:end)-datOut.matJVals(:,:,1:end-nStride) ...
  ,1),2)),1,[]), '^-', 'linewidth', 2, 'markersize', 10 );
grid on;
endif
%
%
%
%
%"C" is for "current".
vecXC = datOut.vecXVals(:,myIter);
vecFC = datOut.vecFVals(:,myIter);
matJC = datOut.matJVals(:,:,myIter);
vecDeltaC = datOut.vecDeltaVals(:,myIter);
vecXVals = datOut.iterDat(myIter).collected_vecXVals;
vecFVals = datOut.iterDat(myIter).collected_vecFVals;
% Move step/jupdate dependent stuff earlier.
%vecXPts = datOut.iterDat(myIter).lesquj_datOut.vecXPts;
%vecFPts = datOut.iterDat(myIter).lesquj_datOut.vecFPts;
%wPts = datOut.iterDat(myIter).lesquj_datOut.wPts;
numVals = size(vecXVals,2);
numPts = size(vecXPts,2);
% We also have vecX0, vecF0, and matJ0.
%
%
%
vecDE = vecXE - vecXC;
assert( norm(vecDE) > eps );
vecVEHat = vecDE/norm(vecDE);
assert( norm(vecFC) > eps );
vecFCHat = vecFC/norm(vecFC);
%
dC = sqrt(sumsq(vecXC-vecXC,1)); % == 0.
deC = sqrt(sumsq(vecXC-vecXE,1));
omegaC = sumsq(vecFC)/2.0;
pC = vecVEHat'*(vecXC-vecXC); % = 0.
fC = vecFCHat'*vecFC;
%
dE = sqrt(sumsq(vecXE-vecXC,1));
deE = sqrt(sumsq(vecXE-vecXE,1)); % = 0.
omegaE = sumsq(vecFE,1)/2.0;
pE = vecVEHat'*(vecXE-vecXC);
fE = vecFCHat'*(vecFE);
%
dVals = sqrt(sumsq(vecXVals-vecXC,1));
deVals = sqrt(sumsq(vecXVals-vecXE,1));
omegaVals = sumsq(vecFVals,1)/2.0;
pVals = vecVEHat'*(vecXVals-vecXC);
fVals = vecFCHat'*(vecFVals);
%
dPts = sqrt(sumsq(vecXPts-vecXC,1));
dePts = sqrt(sumsq(vecXPts-vecXE,1));
omegaPts = sumsq(vecFPts,1)/2.0;
pPts = vecVEHat'*(vecXPts-vecXC);
fPts = vecFCHat'*vecFPts;
%
numPost = 101;
qLo = min(pVals)/norm(vecXE-vecXC);
qHi = max(pVals)/norm(vecXE-vecXC);
qPost = linspace( qLo - 0.3*(qHi-qLo), qHi + 0.3*(qHi-qLo), numPost );
vecXPost = vecXC + ( vecXE - vecXC ) * qPost;
for n=1:numPost
	vecFPost(:,n) = funchFJ( vecXPost(:,n) );
endfor
dPost = sqrt(sumsq(vecXPost-vecXC,1));
dePost = sqrt(sumsq(vecXPost-vecXE,1));
omegaPost = sumsq(vecFPost,1)/2.0;
pPost = vecVEHat'*(vecXPost-vecXC);
fPost = vecFCHat'*vecFPost;
%
numFigs++; figure(numFigs);
plot( ...
  pC, fC, 'p', 'linewidth', 3, 'markersize', 30, ...
  pE, fE, 'x', 'linewidth', 3, 'markersize', 30, ...
  pVals, fVals, 'o', 'linewidth', 1, 'markersize', 10, ...
  pPts, fPts, '+', 'linewidth', 1, 'markersize', 10, ...
  pPost, fPost, 'x-', 'linewidth', 2, 'markersize', 5 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  pC, numVals, 'p', 'linewidth', 3, 'markersize',30, ...
  pE*[1,1], [1,numVals], '-', 'linewidth', 2, ...
  pVals, (1:numVals), 'o-', 'linewidth', 1, 'markersize', 10 );
grid on;
%
%
%
assert( norm(vecDeltaC) > eps );
vecDeltaCHat = vecDeltaC/norm(vecDeltaC);
%
gC = vecDeltaCHat'*(vecXC-vecXC); % = 0;
gE = vecDeltaCHat'*(vecXE-vecXC);
gVals = vecDeltaCHat'*(vecXVals-vecXC);
gPts = vecDeltaCHat'*(vecXPts-vecXC);
gPost = vecDeltaCHat'*(vecXPost-vecXC);
%
numToast = 201;
tLo = min(gVals)/norm(vecDeltaC);
tHi = max(gVals)/norm(vecDeltaC);
tToast = linspace( tLo - 0.3*(tHi-tLo), tHi + 0.3*(tHi-tLo), numToast );
vecXToast = vecXC + vecDeltaC * tToast;
for n=1:numToast
	vecFToast(:,n) = funchFJ( vecXToast(:,n) );
endfor
gToast = vecDeltaCHat'*(vecXToast-vecXC);
dToast = sqrt(sumsq(vecXToast-vecXC,1));
deToast = sqrt(sumsq(vecXToast-vecXE,1));
omegaToast = sumsq(vecFToast,1)/2.0;
pToast = vecVEHat'*(vecXToast-vecXC);
fToast = vecFCHat'*vecFToast;
%
mToast = vecFCHat'*( vecFC + matJC*(vecXToast-vecXC) );
mC = vecFCHat'*vecFC;
%fog = (gVals*(gVals'))\( gVals*((fVals-fC)') );
wModPts = double(gPts>-0.2);
%wModPts = wPts;
fog = ( (gPts.*wModPts)*(gPts') )\( (gPts.*wModPts)*((fPts-fC)') );
nToast = fC + fog*gToast;
%
numFigs++; figure(numFigs);
plot( ...
  gC, omegaC, 'p', 'linewidth', 3, 'markersize', 30, ...
  gE, omegaE, 'x', 'linewidth', 3, 'markersize', 30, ...
  gVals, omegaVals, 'o', 'linewidth', 1, 'markersize', 10, ...
  gPts, omegaPts, '+', 'linewidth', 1, 'markersize', 10, ...
  gToast, omegaToast, 'x-', 'linewidth', 2, 'markersize', 5 );
grid on;
title( "omega" );
%
numFigs++; figure(numFigs);
plot( ...
  gC, fC, 'p', 'linewidth', 3, 'markersize', 30, ...
  gE, fE, 'x', 'linewidth', 3, 'markersize', 30, ...
  gVals, fVals, 'o', 'linewidth', 1, 'markersize', 10, ...
  gPts, fPts, '+', 'linewidth', 1, 'markersize', 10, ...
  gToast, fToast, 'x-', 'linewidth', 2, 'markersize', 5, ...
  gToast, mToast, 's-', 'linewidth', 2, 'markersize', 5, ...
  gToast, nToast, 'v-', 'linewidth', 2, 'markersize', 5 );
grid on;
title( "f" );
%
numFigs++; figure(numFigs);
semilogy( ...
  gPts, wPts, '+', 'linewidth', 1, 'markersize', 10, ...
  gPts, wModPts, '+', 'linewidth', 1, 'markersize', 10 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  gC, numVals, 'p', 'linewidth', 3, 'markersize',30, ...
  gE*[1,1], [1,numVals], '-', 'linewidth', 2, ...
  gVals, (1:numVals), 'o-', 'linewidth', 1, 'markersize', 10 );
grid on;
%
return;
