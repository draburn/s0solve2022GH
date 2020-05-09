myclear;
%
probSize = 5
%
%setprngstates(95746016);
vecX1 = randn(probSize,1)
vecX2 = randn(probSize,1)
%
%setprngstates(0);
vecV = randn(probSize,1);
vecVHat = vecV/norm(vecV)
vecW = randn(probSize,1);
%%vecW -= (vecVHat'*vecW)*vecVHat;
vecWHat = vecW/norm(vecW)
c = abs(randn);
%
alpha = vecVHat'*vecWHat
funchJ = @(x)( c*( vecWHat*(vecVHat'*x) + vecVHat*(vecWHat'*x) ) );
%
vecJ1 = funchJ(vecX1);
vecJ2 = funchJ(vecX2);
%
vecY1 = vecJ1;
vecY1Hat = vecY1/norm(vecY1);
vecY2 = vecJ2;
vecY2 -= (vecJ2'*vecY1Hat)*vecY1Hat;
vecY2Hat = vecY2/norm(vecY2);
rvecPhi = linspace(-pi,pi,10001);
%matX = vecJ1*cos(rvecPhi) + vecJ2*sin(rvecPhi);
matX = vecY1Hat*cos(rvecPhi) + vecY2Hat*sin(rvecPhi);
matJ = funchJ(matX);
%
rvecF0 = sum(matJ.^2,1);
rvecF1 = sum((vecVHat'*matJ).^2,1);
rvecF2 = sum((vecWHat'*matJ).^2,1);
rvecF3 = sum((vecY1Hat'*matJ).^2,1);
rvecF4 = sum((vecY2Hat'*matJ).^2,1);
rvecF5 = sum(matX.*matJ,1);
%
%%[ dummy, indexOfMax ] = max(rvecF1);
%%vecZ = matJ(:,indexOfMax);
%%vecZHat = vecZ/norm(vecZ)
%
plot( ...
  rvecPhi, rvecF0, 'o-', ...
  rvecPhi, rvecF1, 'v-', ...
  rvecPhi, rvecF2, '^-', ...
  rvecPhi, rvecF3, 'v-', ...
  rvecPhi, rvecF4, '^-', ...
  rvecPhi, rvecF5, 's-' );
grid on;
