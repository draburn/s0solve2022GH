clear;
tic();
numFigs = 0;
%
numS1Vals = 51;
numS2Vals = 81;
%
rvecS1Vals = linspace(0.0,1.0,numS1Vals);
rvecS2Vals = linspace(0.0,2.0,numS2Vals);
matV = eye(5,2);
vecV1 = matV(:,1);
vecV2 = matV(:,2);
%
%matC = round(2+(61*rand(numS1Vals,numS2Vals)));
sizeC = 64;
[ matS1, matS2 ] = ndgrid( rvecS1Vals, rvecS2Vals );
matC = cos(2*pi*matS1)+sin(2*pi*matS2);
cMax = max(max(matC));
cMin = min(min(matC));
matC = 2 + ((sizeC-3)*(matC-cMin)/(cMax-cMin));
colmap = jet(sizeC);
%
curveDat = [];
%
numFigs++; figure(numFigs);
vizPlane( ...
  rvecS1Vals, vecV1, rvecS2Vals, vecV2, ...
  matC, colmap, ...
  curveDat );
%
toc();
return;


clear;
tic();
numFigs = 0;
%randnState = mod(round(1E6*time),1E6);
%randnState = 0;
randnState = 399011;
echo_randnState = randnState
randn( "state", randnState );
rand( "state", randnState );
%
%
sizeX = 3;
sizeK = 2;
sizeF = sizeX;
vecX0 = zeros(sizeX,1);
matU = randn(sizeX,sizeK);
%
matJ = randn(sizeF,sizeX);
vecXSecret = matU*rand(sizeK,1)*0.5;
funchOmega = @(matXDummy)( ...
  0.5 * sum( matJ * (matXDummy-repmat(vecXSecret,[1,size(matXDummy,2)])),1 ) ...
  );
%
%
toc;
