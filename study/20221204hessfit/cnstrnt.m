clear;
tic();
setprngstates(0);
numFigs = 0;
numXVals = 101; numYVals = 101;
xLo = -2.0; xHi = 2.0; yLo = -2.0; yHi = 2.0;
%matV = [ 1, 1; 1, 1 ];
matV = [ 1, 1; 0, 1 ];
%matV = [ 1, 0; 0, 1 ];
alpha = (matV(:,1)'*matV(:,2))^2 / ( sumsq(matV(:,1)) * sumsq(matV(:,2)) )
%%%alpha = sqrt(alpha)

matC = [ 1.0, alpha; alpha, 1.0 ];
matCSqrt = [ 1.0, sqrt(alpha); sqrt(alpha), 1.0 ];
vecB = max(abs(matV),[],2);

valsX = linspace( xLo, xHi, numXVals );
valsY = linspace( yLo, yHi, numYVals );
[ meshX, meshY ] = meshgrid( valsX, valsY );
numPts = numXVals * numYVals;
ptsVecZ = [ reshape( meshX, 1, numPts ); reshape( meshY, 1, numPts ) ];
ptsVecDelta = matV*ptsVecZ;
ptsXIndexOfZ = 1 + round( (numXVals-1) * ( ptsVecDelta(1,:) - xLo ) / ( xHi - xLo ) );
ptsYIndexOfZ = 1 + round( (numYVals-1) * ( ptsVecDelta(2,:) - yLo ) / ( yHi - yLo ) );

ptsSatsTrue = ( max( abs(ptsVecDelta)./(eps+vecB), [], 1 ) < 1.0 );
ptsSatsZL1 = ( sum( abs(ptsVecZ), 1 ) < 1.0 );
ptsSatsZL2 = ( sum( ptsVecZ.^2, 1 ) < 1.0 );
ptsSatsCZL2 = sum( ptsVecZ .* ( matC * ptsVecZ ), 1 ) < 1.0;
ptsSatsCZAbsL2 = sum( abs(ptsVecZ) .* ( matC * abs(ptsVecZ) ), 1 ) < 1.0;
ptsSatsCSqrtZL2 = sum( ptsVecZ .* ( matCSqrt * ptsVecZ ), 1 ) < 1.0;

meshSatsTrue = reshape( ptsSatsTrue, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshSatsTrue) );
grid on;
axis equal; axis equal;

meshSatsZL1 = reshape( ptsSatsZL1, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshSatsZL1) );
grid on;
axis equal; axis equal;

meshSatsZL2 = reshape( ptsSatsZL2, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshSatsZL2) );
grid on;
axis equal; axis equal;

meshSatsCZL2 = reshape( ptsSatsCZL2, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshSatsCZL2) );
grid on;
axis equal; axis equal;

meshSatsCZAbsL2 = reshape( ptsSatsCZAbsL2, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshSatsCZAbsL2) );
grid on;
axis equal; axis equal;

meshSatsCSqrtZL2 = reshape( ptsSatsCSqrtZL2, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshSatsCSqrtZL2) );
grid on;
axis equal; axis equal;








return;


clear;
tic();
setprngstates(0);
numFigs = 0;

numXVals = 51; numYVals = 51;
%xLo = -2.0; xHi = 2.0; yLo = -2.0; yHi = 2.0;
xLo = -3.0; xHi = 3.0; yLo = -3.0; yHi = 3.0;
matV = [ 1, 1; 1, 1 ];
%matV = [ 1, 0; 0, 1 ];
%matV = [ 1, 1; 1, 2 ];
alpha = (matV(:,1)'*matV(:,2))^2 / ( sumsq(matV(:,1)) * sumsq(matV(:,2)) );
alpha = sqrt(alpha)

matC = [ 1.0, alpha; alpha, 1.0 ];
vecB = max(abs(matV),[],2);

valsX = linspace( xLo, xHi, numXVals );
valsY = linspace( yLo, yHi, numYVals );
[ meshX, meshY ] = meshgrid( valsX, valsY );
numPts = numXVals * numYVals;
ptsVecZ = [ reshape( meshX, 1, numPts ); reshape( meshY, 1, numPts ) ];
ptsVecDelta = matV*ptsVecZ;
ptsXIndexOfDelta = 1 + round( (numXVals-1) * ( ptsVecDelta(1,:) - xLo ) / ( xHi - xLo ) );
ptsYIndexOfDelta = 1 + round( (numYVals-1) * ( ptsVecDelta(2,:) - yLo ) / ( yHi - yLo ) );

ptsSatsTrue = ( max( abs(ptsVecDelta)./(eps+vecB), [], 1 ) < 1.0 );
meshCountsTrue = zeros( size(meshX) );
for p = (1 : numPts)(ptsSatsTrue)
	meshCountsTrue( ptsXIndexOfDelta(p), ptsYIndexOfDelta(p) )++;
endfor
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshCountsTrue>0.5) );
grid on;
axis equal; axis equal;

ptsSatsZL1 = ( sum( abs(ptsVecZ), 1 ) < 1.0 );
meshCountsZL1 = zeros( size(meshX) );
for p = (1 : numPts)(ptsSatsZL1)
	meshCountsZL1( ptsXIndexOfDelta(p), ptsYIndexOfDelta(p) )++;
endfor
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshCountsZL1>0.5) );
grid on;
axis equal; axis equal;

ptsSatsZL2 = ( sum( ptsVecZ.^2, 1 ) < 1.0 );
meshCountsZL2 = zeros( size(meshX) );
for p = (1 : numPts)(ptsSatsZL2)
	meshCountsZL2( ptsXIndexOfDelta(p), ptsYIndexOfDelta(p) )++;
endfor
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshCountsZL2>0.5) );
grid on;
axis equal; axis equal;

ptsSatsCZL2 = ( sum( ptsVecZ .* ( matC * ptsVecZ ), 1 ) < 1.0 );
meshCountsCZL2 = zeros( size(meshX) );
for p = (1 : numPts)(ptsSatsCZL2)
	meshCountsCZL2( ptsXIndexOfDelta(p), ptsYIndexOfDelta(p) )++;
endfor
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshCountsZL2>0.5) );
grid on;
axis equal; axis equal;

ptsSatsCZAbsL2 = ( sum( abs(ptsVecZ) .* ( matC * abs(ptsVecZ) ), 1 ) < 1.0 );
meshCountsCZAbsL2 = zeros( size(meshX) );
for p = (1 : numPts)(ptsSatsCZAbsL2)
	meshCountsCZAbsL2( ptsXIndexOfDelta(p), ptsYIndexOfDelta(p) )++;
endfor
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshCountsZL2>0.5) );
grid on;
axis equal; axis equal;

toc();
return



numPts = numXVals * numYVals;

matV = [ 1, 1; 1, 1 ];
alpha = (matV(:,1)'*matV(:,2))^2 / ( sumsq(matV(:,1)) * sumsq(matV(:,2)) );
matC = [ 1.0, alpha; alpha, 1.0 ];
vecB = max(abs(matV),[],2);
ptsVecDelta = matV*ptsVecZ;
ptsXIndex = 1 + round( (numXVals-1) * ( ptsVecDelta(1,:) - xLo ) / ( xHi - xLo ) );
ptsYIndex = 1 + round( (numYVals-1) * ( ptsVecDelta(2,:) - yLo ) / ( yHi - yLo ) );

ptsSatsTrue = max( abs(ptsVecDelta)./(eps+vecB), 1 ) < 1.0;
meshHitsTrue = zeros( size(meshX) );


return;

ptsNumHitsTrue = zeros( 1, numPts );




meshHitTrue = zeros(size(meshX));

return;


ptsSatsTrue = max( abs(ptsVecDelta)./(eps+vecB), 1 ) < 1.0;

meshHitTrue = zeros(size(meshX));



return;
ptsSatsZL1 = sum( abs(ptsVecZ), 1 ) < 1.0;
ptsSatsZL2 = sum( ptsVecZ.^2, 1 ) < 1.0;
ptsSatsCZL2 = sum( ptsVecZ .* ( matC * ptsVecZ ), 1 ) < 1.0;

ptsHitZL1 = zeros(1,numPts);
ptsHitZL2 = zeros(1,numPts);
ptsHitCZL2 = zeros(1,numPts);




return;

matV = [ 1, 1; 1, 1 ];
%matV = [ 1, 0; 0, 1 ];
matVPInv = pinv(matV);
ptsVecZ = matVPInv * ptsVecDelta;

ptsBetaTrue = sum(abs(ptsVecZ),1).^2;
meshBetaTrue = reshape( ptsBetaTrue, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshBetaTrue<=1.0) );
grid on;
axis equal; axis equal;

ptsBetaCrude = sum( ptsVecZ.^2, 1 );
meshBetCrude = reshape( ptsBetaCrude, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshBetCrude<=1.0) );
grid on;
axis equal; axis equal;

alpha = (matV(:,1)'*matV(:,2))^2 / ( sumsq(matV(:,1)) * sumsq(matV(:,2)) );
alpha = 0.0;
matC = [ 1.0, alpha; alpha, 1.0 ]
ptsBetaModel = sum( ptsVecZ .* (matC*ptsVecZ), 1 );
meshBetaModel = reshape( ptsBetaModel, size(meshX) );
numFigs++; figure(numFigs);
contourf( meshX, meshY, double(meshBetaModel<=1.0) );
grid on;
axis equal; axis equal;

return;

alpha = cos( matV(:,1)'*matV(:,2) )^2 / ( sumsq(matV(:,1)) * sumsq(matV(:,2)) );
matC = [ 1.0, alpha; alpha, 1.0 ];
vecB = max(abs(matV),[],2);
ptsVecDelta = matV*ptsVecR;


ptsGood = ptsVec


return;

ptsBetaTrue = max( ptsVecDelta ./ vecB ).^2;
ptsBetaWant = sum( abs(ptsR), 1 ).^2;
ptsBetaCrud = sum( ptsR.*ptsR, 1 );
%%%ptsBetaPred = sum( ptsR .* (matC*ptsR), 1 );
ptsBetaPred = sum( abs(ptsR) .* (matC*abs(ptsR)), 1 );

meshBetaTrue = reshape( ptsBetaTrue, size(meshX) );
meshBetaWant = reshape( ptsBetaWant, size(meshX) );
meshBetaCrud = reshape( ptsBetaCrud, size(meshX) );
meshBetaPred = reshape( ptsBetaPred, size(meshX) );

numFigs++; figure(numFigs);
contourf( meshX, meshY, meshBetaTrue );
grid on;
axis equal; axis equal;

numFigs++; figure(numFigs);
contourf( meshX, meshY, meshBetaWant );
grid on;
axis equal; axis equal;

numFigs++; figure(numFigs);
contourf( meshX, meshY, meshBetaCrud );
grid on;
axis equal; axis equal;

numFigs++; figure(numFigs);
contourf( meshX, meshY, meshBetaPred );
grid on;
axis equal; axis equal;




return;
%
%
fp = zeros(1,size(vp,2));
gp = zeros(1,size(vp,2));
hp = zeros(1,size(vp,2));
%
%
numTrials = 1;
for t = 1 : numTrials
	%V = randn( 2, 2 );
	V = [ 1, 1; 2, 3 ];
	alpha = cos( V(:,1)'*V(:,2) ) / ( norm(V(:,1)) * norm(V(:,2)) );
	alpha = sqrt(abs(alpha))
	C = [ 1.0, alpha; alpha, 1.0 ];
	sqrtC = sqrt(C);
	fp += double( sum( vp .* (C*vp), 1 ) < 1.0 );
	gp += double( sum( abs(sqrtC*vp), 1 ) < 1.0 );
endfor
fp /= numTrials;
gp /= numTrials;
ff = reshape( fp, size(xx) );
gg = reshape( gp, size(xx) );
%
%
numFigs++; figure(numFigs);
contourf( xx, yy, ff );
grid on;
axis equal;
%
%
numFigs++; figure(numFigs);
contourf( xx, yy, gg );
grid on;
axis equal;
