clear;
commondefs;
thisFile = "test_vizLLMCurvs";
%
randnSeed = mod(round(1E6*time),1E6);
%randnSeed = 0;
%randnSeed = 792958; % Newton is best.
%randnSeed = 44936; % Grad Dir is best.
%randnSeed = 952130; % Overshoot.
echo_randnSeed = randnSeed
wNoiseLevel = 0.5;
sizeX = 10;
sizeF = 10;
sizeK = 2;
numPts = 50;

randn("seed",randnSeed);
vecF = randn(sizeF,1);
matJ = randn(sizeF,sizeX);
funchF = @(v)( repmat(vecF,[1,size(v,2)]) + (matJ*v) );
matV = eye(sizeX,sizeK);
vecX = zeros(sizeX,1);
vecXSecret = -matJ \ vecF;
%
matW = matJ*matV;
matW += wNoiseLevel * randn(size(matW));
prm.curveTypes = [ ...
  STEPTYPE__NEWTON, ...
  STEPTYPE__PICARD, ...
  STEPTYPE__PICARD_SCALED, ...
  STEPTYPE__GRADDIR, ...
  STEPTYPE__GRADDIR_SCALED, ...
  STEPTYPE__LEVCURVE, ...
  STEPTYPE__LEVCURVE_SCALED, ...
  STEPTYPE__GRADCURVE, ...
  STEPTYPE__GRADCURVE_SCALED ];
vizLLMCurves( funchF, vecX, matV, matW, numPts, vecXSecret );
