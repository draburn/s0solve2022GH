clear;
commondefs;
thisFile = "test_getLocLinModelCurves";
tic();
%
%%%randnSeed = mod(round(1E6*time()),1E6);
randnSeed = 1365; % Non-monotonic LM.
%%%randnSeed = 849761; % Causes GC flinspace to fail.
echo__randnSeed = randnSeed
randn("seed",randnSeed);
sizeF = 2;
sizeX = 2;
sizeK = 2;
matS = [1,0;0,10];
vecF = matS*randn(sizeF,1);
matV = eye(sizeX,sizeK);
matW = matS*randn(sizeF,sizeK);
numPts = 50;
%
[ curveDat, retCode, datOut ] = getLocLinModelCurves( vecF, matV, matW, numPts );
assert( RETCODE__SUCCESS == retCode );
matDeltaG = curveDat.matDeltaG;
matDeltaF = curveDat.matDeltaF;
matDeltaN = curveDat.matDeltaN;
matDeltaL = curveDat.matDeltaL;
matDeltaLM = curveDat.matDeltaLM;
matDeltaGC = curveDat.matDeltaGC;
%
matF0 = repmat(vecF,[1,numPts]);
matFPG = matF0 + (matW * (matV'*matDeltaG));
matFPF = matF0 + (matW * (matV'*matDeltaF));
matFPN = matF0 + (matW * (matV'*matDeltaN));
matFPL  = repmat(vecF,[1,size(matDeltaL, 2)]) + (matW * (matV'*matDeltaL ));
matFPLM = repmat(vecF,[1,size(matDeltaLM,2)]) + (matW * (matV'*matDeltaLM));
matFPGC = repmat(vecF,[1,size(matDeltaGC,2)]) + (matW * (matV'*matDeltaGC));
deltaGNormVals = sqrt(sum(matDeltaG.^2,1));
deltaFNormVals = sqrt(sum(matDeltaF.^2,1));
deltaNNormVals = sqrt(sum(matDeltaN.^2,1));
deltaLNormVals = sqrt(sum(matDeltaL.^2,1));
deltaLMNormVals = sqrt(sum(matDeltaLM.^2,1));
deltaGCNormVals = sqrt(sum(matDeltaGC.^2,1));
omegaGVals = 0.5*(sum(matFPG.^2,1));
omegaFVals = 0.5*(sum(matFPF.^2,1));
omegaNVals = 0.5*(sum(matFPN.^2,1));
omegaLVals = 0.5*(sum(matFPL.^2,1));
omegaLMVals = 0.5*(sum(matFPLM.^2,1));
omegaGCVals = 0.5*(sum(matFPGC.^2,1));
%
numFigs = 0;
xn = 1;
yn = 2;
colorG = [ 1.0, 0.0, 0.0 ];
colorF = [ 0.0, 0.0, 1.0 ];
colorN = [ 0.0, 0.8, 0.0 ];
colorL = [ 0.8, 0.8, 0.0 ];
colorLM = [ 0.8, 0.4, 1.0 ];
colorGC = [ 0.4, 0.6, 0.6 ];
%
numFigs++; figure(numFigs);
plot( ...
  deltaGNormVals, omegaGVals, 'o-', 'color', colorG, ...
  deltaFNormVals, omegaFVals, 'x-', 'color', colorF, ...
  deltaNNormVals, omegaNVals, 's-', 'color', colorN, ...
  deltaLNormVals, omegaLVals, '^-', 'color', colorL, ...
  deltaLMNormVals, omegaLMVals, 'v-', 'color', colorLM, ...
  deltaGCNormVals, omegaGCVals, '*-', 'color', colorGC, 'markersize', 15  );
grid on;
%
xMin = max(min([ ...
  matDeltaG(xn,:), ...
  matDeltaF(xn,:), ...
  matDeltaN(xn,:) ]));
xMax = max(max([ ...
  matDeltaG(xn,:), ...
  matDeltaF(xn,:), ...
  matDeltaN(xn,:) ]));
yMin = max(min([ ...
  matDeltaG(yn,:), ...
  matDeltaF(yn,:), ...
  matDeltaN(yn,:) ]));
yMax = max(max([ ...
  matDeltaG(yn,:), ...
  matDeltaF(yn,:), ...
  matDeltaN(yn,:) ]));
if (1)
aMin = min([xMin,yMin]);
aMax = max([xMax,yMax]);
xLo = aMin-0.3*abs(aMax-aMin);
xHi = aMax+0.3*abs(aMax-aMin);
yLo = xLo;
yHi = xHi;
else
xLo = xMin-0.3*abs(xMax-xMin);
xHi = xMax+0.3*abs(xMax-xMin);
yLo = yMin-0.3*abs(yMax-yMin);
yHi = yMax+0.3*abs(yMax-yMin);
end
%
numXVals = 50;
numYVals = 51;
xVals = xLo + ((xHi-xLo)*(0:numXVals-1)/(numXVals-1.0));
yVals = yLo + ((yHi-yLo)*(0:numYVals-1)/(numYVals-1.0));
[ xGrid, yGrid ] = meshgrid( xVals, yVals );
%
xv = 1;
yv = 2;
omegaGrid = 0.5*( ...
   (( (vecF(xv) + (matW(xv,xn)*xGrid) + (matW(xv,yn)*yGrid)) ).^2) ...
 + (( (vecF(yv) + (matW(yv,xn)*xGrid) + (matW(yv,yn)*yGrid)) ).^2) );
%
numFigs++; figure(numFigs);
contour( xGrid, yGrid, sqrt(omegaGrid), 50 );
hold on;
plot( ...
  matDeltaG(xn,:), matDeltaG(yn,:), 'o-', 'color', colorG, ...
  matDeltaF(xn,:), matDeltaF(yn,:), 'x-', 'color', colorF, ...
  matDeltaN(xn,:), matDeltaN(yn,:), 's-', 'color', colorN, ...
  matDeltaL(xn,:), matDeltaL(yn,:), '^-', 'color', colorL, ...
  matDeltaLM(xn,:), matDeltaLM(yn,:), '^-', 'color', colorLM, ...
  matDeltaGC(xn,:), matDeltaGC(yn,:), '*-', 'color', colorGC, 'markersize', 15 );
hold off;
axis equal;
grid on;
%
toc();
